#!/usr/bin/env python3

# Author: Austin Marshall
# Date: 26/Oct/25

"""
Pigeon - post-hoc branchwater/sourmash comparison (unitigs from GFA)

Compares:
  * Unitigs (from: gfatools asm -u <assembly.gfa>)
  * Assembly contigs (FASTA or .sig.zip)
  * Bins (FASTA directory or .sig.zip of bins)

Outputs an HTML with:
  * 1−Jaccard heatmap + MDS (unitigs/assembly/bins)
  * Unitigs partition: A-only / B-only / both / unexplained (+ Sankey)
  * Top-bins by explained unitig hashes
  * Greedy cumulative curve: fraction of unitig hashes explained vs #bins
  * Summary metrics table incl. AUC + PAM score

Inputs:
  --gfa assembly.gfa   (required unless --unitigs-fa/--unitigs-sig given)
  --assembly contigs.fasta(.gz) or .sig.zip
  --bins_dir dir of bin FASTAs or a .sig.zip of bins

Dependencies:
    * sourmash
    * sourmash-minimal
    * biopython
    * gfatools (for unitig extraction)
    * plotly
    * python-rocksdb
    * numpy
    * tqdm
    * pandas
    * pigz

"""

import os, sys, json, argparse, subprocess, gzip, struct, shutil
from pathlib import Path
from typing import Dict, List, Tuple, Set, Iterable, Optional

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sourmash import load_file_as_signatures, MinHash
from sourmash.signature import SourmashSignature

# Optional RocksDB
try:
    import rocksdb  # type: ignore
    HAVE_ROCKS = True
except Exception:
    HAVE_ROCKS = False


# ---------------- I/O helpers ----------------
def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)

def safe_name(s: str) -> str:
    allowed = set("-_.")
    return "".join(ch if ch.isalnum() or ch in allowed else "_" for ch in s)[:200]

def short_name(s: str) -> str:
    return Path(s).name


# ---------------- GFA -> unitigs ----------------
def extract_unitigs_from_gfa(gfa_path: str, out_fa: Path):
    """
    Ensure FASTA output for unitigs.
    `gfatools asm -u` produces a GFA then convert to FASTA.
    """
    # Run asm -u and capture stdout
    p1 = subprocess.run(
        ["gfatools", "asm", "-u", gfa_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
    )
    data = p1.stdout

    # FASTA should have '>' at start or soon after newlines
    def looks_like_fasta(b: bytes) -> bool:
        return b.startswith(b">") or b"\n>" in b[:1_000_000]

    if not looks_like_fasta(data):
        # Convert GFA -> FASTA
        p2 = subprocess.run(
            ["gfatools", "gfa2fa", "-"],
            input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
        data = p2.stdout

    # Basic sanity check
    if not data or not looks_like_fasta(data):
        raise RuntimeError(
            "Unitig conversion did not yield FASTA. "
            "Check that your GFA has segment sequences and try upgrading gfatools."
        )

    with open(out_fa, "wb") as w:
        w.write(data)


# ---------------- Sketching ----------------
def singlesketch(path: str, name: str, ks: List[int], scaled: int, seed: int, out_sig: Path):
    cmd = ["sourmash", "scripts", "singlesketch", path, "-o", str(out_sig), "--name", name]
    for ki in ks:
        cmd += ["-p", f"k={ki},scaled={scaled},seed={seed}"]
    subprocess.run(cmd, check=True)

def sketch_many(files: List[str], ks: List[int], scaled: int, seed: int, out_sig: Path):
    lst = out_sig.with_suffix(".list.txt")
    with open(lst, "w") as f:
        for p in files:
            f.write(str(p) + "\n")
    cmd = [
        "sourmash", "sketch", "dna",
        "--from-file", str(lst),
        "-o", str(out_sig),
        "--name-from-first",
    ]
    for ki in ks:
        cmd += ["-p", f"k={ki},scaled={scaled},seed={seed}"]
    subprocess.run(cmd, check=True)
    try:
        os.remove(lst)
    except Exception:
        pass


# ---------------- Signatures & math ----------------
def load_sig_zip(sig_zip: str, ksize: int) -> List[SourmashSignature]:
    return list(load_file_as_signatures(sig_zip, ksize=ksize))

def get_hash_set(sig: SourmashSignature) -> Set[int]:
    mh = sig.minhash
    return set(mh.hashes.keys()) if hasattr(mh, "hashes") else set(mh.get_mins())

def jaccard(A: Set[int], B: Set[int]) -> float:
    if not A and not B:
        return 1.0
    if not A or not B:
        return 0.0
    I = len(A & B); U = len(A | B)
    return I / U if U else 0.0


# ---------------- RocksDB helpers ----------------
def to_key(h: int) -> bytes:
    # store as unsigned 64-bit big-endian
    return struct.pack(">Q", h & 0xFFFFFFFFFFFFFFFF)

class HashDB:
    """Single-CF RocksDB with namespaced keys:
       b'U' (unitigs), b'A' (assembly), b'B' (bins_union).
    """
    def __init__(self, dbpath: Path):
        self.ok = HAVE_ROCKS
        self.db = None
        if not self.ok:
            print("[DB] python-rocksdb not available; skipping DB creation.")
            return
        ensure_dir(dbpath)
        opts = rocksdb.Options()
        opts.create_if_missing = True
        opts.max_open_files = 64
        self.db = rocksdb.DB(str(dbpath), opts)

    def put_set(self, prefix: bytes, S: Set[int]):
        if not self.ok or self.db is None:
            return
        wb = rocksdb.WriteBatch()
        for h in S:
            wb.put(prefix + to_key(h), b"1")
        self.db.write(wb)


# ---------------- Metric compute ----------------
def compute_novel_metrics(sets: Dict[str, Set[int]], bin_names: List[str], primary_key: str):
    Uset = sets.get(primary_key, set())          # unitigs
    A     = sets.get("assembly", set())          # assembly contigs
    Bunion = set().union(*(sets[n] for n in bin_names)) if bin_names else set()

    Uall = A | Bunion
    UA   = Uset & A
    UB   = Uset & Bunion
    U_only    = Uset - Uall
    U_A_not_B = UA - Bunion
    U_B_not_A = UB - A
    U_both    = Uset & A & Bunion

    # Greedy cumulative coverage of primary (unitigs) by bins
    covered = set()
    curve = []
    ordered_bins = sorted(bin_names, key=lambda b: len((sets[b] & Uset) - covered), reverse=True)
    for i, b in enumerate(ordered_bins, 1):
        new = (sets[b] & Uset) - covered
        covered |= new
        curve.append((i, b, len(covered)))

    per_bin_explained = [(b, len(sets[b] & Uset)) for b in bin_names]
    per_bin_explained.sort(key=lambda x: x[1], reverse=True)

    Un = len(Uset) if Uset else 1

    # AUC of cumulative explained fraction vs bins (normalized by mean(y))
    if curve:
        ys = [v / float(Un) for _, _, v in curve]
        auc_norm = float(sum(ys)) / float(len(ys))
    else:
        auc_norm = 0.0

    # Pigeon Appraisal Metric (PAM) in [0,1]
    frac_U_only  = (len(U_only) / Un)
    frac_U_in_A  = (len(UA)     / Un)
    frac_U_in_B  = (len(UB)     / Un)
    frac_U_Aonly = (len(U_A_not_B) / Un)
    frac_U_Bonly = (len(U_B_not_A) / Un)
    frac_U_both  = (len(U_both) / Un)

    explained = 1.0 - frac_U_only
    balance   = 1.0 - abs(frac_U_Aonly - frac_U_Bonly)
    pam = 0.6 * explained + 0.3 * auc_norm + 0.1 * balance

    return {
        f"|{primary_key}|": len(Uset), "|A|": len(A), "|B|": len(Bunion), "|U|": len(Uall),
        f"|{primary_key}∩A|": len(UA), f"|{primary_key}∩B|": len(UB),
        f"frac_{primary_key}_in_A": frac_U_in_A,
        f"frac_{primary_key}_in_B": frac_U_in_B,
        f"frac_{primary_key}_only": frac_U_only,
        f"frac_{primary_key}_A_not_B": frac_U_Aonly,
        f"frac_{primary_key}_B_not_A": frac_U_Bonly,
        f"frac_{primary_key}_both": frac_U_both,
        "counts": {
            "A_only": len(U_A_not_B),
            "B_only": len(U_B_not_A),
            "A_not_B": len(U_both),
            "unexplained": len(U_only),
        },
        "cum_bins_curve": curve,
        "per_bin_explained": per_bin_explained[:50],
        "ordered_bins": ordered_bins,
        "auc_norm": auc_norm,
        "pam": pam,
    }


# ---------------- Report ----------------
def make_report(outdir: Path,
                names: List[str],
                sets: Dict[str, Set[int]],
                sizes: Dict[str, int],
                ksize: int, scaled: int, seed: int,
                novel_metrics: Dict,
                primary_key: str,
                top_bins: int = 20):
    ensure_dir(outdir)
    labels = [short_name(n) for n in names]

    # Pairwise Jaccard & MDS
    N = len(names)
    J = np.ones((N, N), dtype=float)
    for i in range(N):
        for j in range(i + 1, N):
            v = jaccard(sets[names[i]], sets[names[j]])
            J[i, j] = v; J[j, i] = v
    D = 1.0 - J
    J2 = D ** 2
    H = np.eye(N) - np.ones((N, N)) / N
    B = -0.5 * H @ J2 @ H
    w, V = np.linalg.eigh(B)
    idx = np.argsort(-w)
    vals = np.maximum(w[idx[:2]], 0)
    coords = V[:, idx[:2]] * np.sqrt(vals + 1e-12)

    title_primary = "Unitigs" if primary_key == "unitigs" else primary_key

    fig = make_subplots(
        rows=3, cols=3,
        specs=[
            [{"type": "heatmap"}, {"type": "scatter"}, {"type": "bar"}],
            [{"type": "bar"}, {"type": "scatter"}, {"type": "table"}],
            [{"type": "table"}, {"type": "bar"}, {"type": "table"}],
        ],
        column_widths=[0.38, 0.32, 0.30],
        subplot_titles=(
            f"1 - Jaccard ({title_primary}/assembly/bins)",
            "MDS on 1 - Jaccard",
            "Top bins by explained hashes",
            f"{title_primary} partition (A only / B only / both / unexplained)",
            "Cumulative explained vs #bins",
            "Key metrics",
            f"{title_primary} → {{A only, A not B, B only, Unexplained}} (Sankey)",
            f"Per-bin explained (top bins, % of {title_primary})",
            "Parameters",
        ),
    )

    # Heatmap & MDS
    fig.add_trace(go.Heatmap(z=D, x=labels, y=labels, colorscale="Viridis",
                             colorbar=dict(title="1-J")), row=1, col=1)
    fig.add_trace(go.Scatter(x=coords[:, 0], y=coords[:, 1], mode="markers+text",
                             text=[short_name(n) for n in names],
                             textposition="top center"), row=1, col=2)

    # Top bins bar
    top = novel_metrics.get("per_bin_explained", [])[:top_bins]
    fig.add_trace(go.Bar(x=[b for b, _ in top], y=[v for _, v in top], showlegend=False), row=1, col=3)
    fig.update_yaxes(row=1, col=3, title="Hashes explained (count)")

    # Primary partition bar (percent)
    parts = [
        ("A only", novel_metrics.get(f"frac_{primary_key}_A_not_B", 0.0)),
        ("B only", novel_metrics.get(f"frac_{primary_key}_B_not_A", 0.0)),
        ("A not B", novel_metrics.get(f"frac_{primary_key}_both", 0.0)),
        ("Unexplained", novel_metrics.get(f"frac_{primary_key}_only", 0.0)),
    ]
    fig.add_trace(go.Bar(x=[p[0] for p in parts],
                         y=[p[1] for p in parts],
                         text=[f"{p[1]*100:.1f}%" for p in parts],
                         textposition="auto",
                         showlegend=False), row=2, col=1)
    fig.update_yaxes(row=2, col=1, tickformat=".0%")

    # Cumulative explained curve (percent)
    cum = novel_metrics.get("cum_bins_curve", [])
    total_primary = float(novel_metrics.get(f"|{primary_key}|", 1))
    if cum:
        xs = [i for i, _, _ in cum]
        ys = [v / total_primary for _, _, v in cum]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines+markers", name="explained"), row=2, col=2)
        fig.update_xaxes(row=2, col=2, title="# bins")
        fig.update_yaxes(row=2, col=2, title=f"Explained {title_primary}", tickformat=".0%")
    else:
        fig.add_trace(go.Scatter(x=[], y=[]), row=2, col=2)

    # Metrics (AUC + PAM)
    rows = [
        f"|{title_primary}| (hashes)", "|A|", "|B| (union)",
        f"{title_primary} in A (frac)", f"{title_primary} in B (frac)", "Unexplained (frac)",
        "A only (frac)", "B only (frac)", "A not B (frac)",
        "AUC (cum curve)", "PAM score"
    ]
    vals = [
        str(novel_metrics.get(f"|{primary_key}|", 0)),
        str(novel_metrics.get("|A|", 0)),
        str(novel_metrics.get("|B|", 0)),
        f"{novel_metrics.get(f'frac_{primary_key}_in_A', 0.0):.3f}",
        f"{novel_metrics.get(f'frac_{primary_key}_in_B', 0.0):.3f}",
        f"{novel_metrics.get(f'frac_{primary_key}_only', 0.0):.3f}",
        f"{novel_metrics.get(f'frac_{primary_key}_A_not_B', 0.0):.3f}",
        f"{novel_metrics.get(f'frac_{primary_key}_B_not_A', 0.0):.3f}",
        f"{novel_metrics.get(f'frac_{primary_key}_both', 0.0):.3f}",
        f"{novel_metrics.get('auc_norm', 0.0):.3f}",
        f"{novel_metrics.get('pam', 0.0):.3f}",
    ]
    fig.add_trace(go.Table(header=dict(values=["Metric", "Value"]),
                           cells=dict(values=[rows, vals])), row=2, col=3)

    # Sankey: primary -> {A only, A not B, B only, Unexplained}
    counts = novel_metrics.get("counts", {})
    s_labels  = [title_primary, "A only", "A not B", "B only", "Unexplained"]
    s_sources = [0, 0, 0, 0]
    s_targets = [1, 2, 3, 4]
    s_values  = [
        int(counts.get("A_only", 0)),
        int(counts.get("A_not_B", 0)),
        int(counts.get("B_only", 0)),
        int(counts.get("unexplained", 0)),
    ]
    fig.add_trace(go.Sankey(
        arrangement="snap",
        node=dict(label=s_labels, pad=10, thickness=12),
        link=dict(source=s_sources, target=s_targets, value=s_values)
    ), row=3, col=1)

    # Per-bin explained as % of primary
    if top:
        xs = [b for b, _ in top]
        ys = [v / total_primary for _, v in top]
        fig.add_trace(go.Bar(x=xs, y=ys, showlegend=False,
                             text=[f"{y*100:.1f}%" for y in ys],
                             textposition="auto"), row=3, col=2)
        fig.update_yaxes(row=3, col=2, tickformat=".0%", title=f"% of {title_primary} explained")

    # Params table
    p_rows = ["ksize", "scaled", "seed", "#bins", "RocksDB"]
    p_vals = [
        str(ksize), str(scaled), str(seed),
        str(len([n for n in names if n not in (primary_key, "assembly")])),
        "yes" if HAVE_ROCKS else "no",
    ]
    fig.add_trace(go.Table(header=dict(values=["Param", "Value"]),
                           cells=dict(values=[p_rows, p_vals])), row=3, col=3)

    fig.update_layout(
        height=1100, template="plotly_white",
        font=dict(size=12),
        margin=dict(l=60, r=20, t=60, b=60),
        title=f"Post-hoc comparison of {title_primary}, assembly, and bins"
    )
    fig.write_html(str(outdir / "report.html"), include_plotlyjs="cdn")


# ---------------- Main ----------------
def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Primary source is now unitigs extracted from a GFA (or provided directly)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--gfa", help="assembly .gfa to extract unitigs via 'gfatools asm -u'")
    src.add_argument("--unitigs-fa", help="unitigs FASTA(.gz) if already extracted")
    src.add_argument("--unitigs-sig", help="prebuilt sourmash .sig.zip of unitigs")

    ap.add_argument("--assembly", required=True,
                    help="assembly contigs FASTA(.gz) or sourmash .sig.zip")
    ap.add_argument("--bins_dir", required=True,
                    help="directory of bin FASTAs(.gz) or a signature zip of bins")
    ap.add_argument("--outdir", required=True)

    ap.add_argument("--ksize", type=int, default=17)
    ap.add_argument("--scaled", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cores", type=int, default=max(1, (os.cpu_count() or 2)//2))
    ap.add_argument("--skip-db", action="store_true",
                    help="don’t build RocksDB even if available")
    ap.add_argument("--top-bins", type=int, default=20,
                    help="number of bins to show in plots")

    args = ap.parse_args()
    outdir = Path(args.outdir); ensure_dir(outdir)
    sigdir = outdir / "sigs"; ensure_dir(sigdir)
    tmpdir = outdir / "tmp"; ensure_dir(tmpdir)
    dbdir  = outdir / "hashdb"

    def is_sig(p: str) -> bool:
        return str(p).endswith(('.sig', 'sig.zip', '.zip'))

    # Build / gather signatures
    asm_sig    = sigdir / "assembly.sig.zip"
    unitig_sig = sigdir / "unitigs.sig.zip"
    bins_sig   = sigdir / "bins.sig.zip"

    # 1) Assembly (contigs)
    if is_sig(args.assembly):
        asm_sig = Path(args.assembly)
    else:
        print("[sketch] assembly (contigs)")
        singlesketch(args.assembly, "assembly", [args.ksize], args.scaled, args.seed, asm_sig)

    # 2) Unitigs
    if args.unitigs_sig:
        unitig_sig = Path(args.unitigs_sig)
    else:
        unitigs_fa: Optional[Path] = None
        if args.unitigs_fa:
            unitigs_fa = Path(args.unitigs_fa)
        elif args.gfa:
            unitigs_fa = tmpdir / "unitigs.from_gfa.fasta"
            print(f"[gfatools] extracting unitigs from GFA -> {unitigs_fa}")
            extract_unitigs_from_gfa(args.gfa, unitigs_fa)
        else:
            print("One of --gfa, --unitigs-fa, or --unitigs-sig is required."); sys.exit(1)

        print("[sketch] unitigs")
        singlesketch(str(unitigs_fa), "unitigs", [args.ksize], args.scaled, args.seed, unitig_sig)

    # 3) Bins
    if is_sig(args.bins_dir):
        bins_sig = Path(args.bins_dir)
    else:
        print("[sketch] bins")
        bin_files = []
        for f in sorted(Path(args.bins_dir).iterdir()):
            if f.is_file() and any(str(f).endswith(ext) for ext in (".fa",".fna",".fasta",".fa.gz",".fna.gz",".fasta.gz")):
                bin_files.append(str(f))
        if not bin_files:
            print("No bin FASTAs found in --bins_dir"); sys.exit(1)
        sketch_many(bin_files, [args.ksize], args.scaled, args.seed, bins_sig)

    # Load signatures -> sets
    primary_key = "unitigs"
    names = [primary_key, "assembly"]
    sets: Dict[str, Set[int]] = {}
    sizes: Dict[str, int] = {}

    sig_unitigs = load_sig_zip(str(unitig_sig), ksize=args.ksize)[0]
    sets[primary_key] = get_hash_set(sig_unitigs); sizes[primary_key] = len(sets[primary_key])

    sig_asm = load_sig_zip(str(asm_sig), ksize=args.ksize)[0]
    sets["assembly"] = get_hash_set(sig_asm); sizes["assembly"] = len(sets["assembly"])

    bin_sigs = load_sig_zip(str(bins_sig), ksize=args.ksize)
    bin_names: List[str] = []
    for sig in bin_sigs:
        nm = sig.name or sig.filename or "bin"
        nm = safe_name(short_name(nm))
        base = nm; i = 1
        while nm in sets:
            nm = f"{base}_{i}"; i += 1
        S = get_hash_set(sig)
        sets[nm] = S; sizes[nm] = len(S); bin_names.append(nm)
    names += bin_names

    # Build RocksDB unions (optional)
    if HAVE_ROCKS and not args.skip_db:
        print("[DB] building RocksDB unions")
        db = HashDB(dbdir)
        db.put_set(b'U', sets[primary_key])
        db.put_set(b'A', sets["assembly"])
        Bunion = set().union(*(sets[n] for n in bin_names)) if bin_names else set()
        db.put_set(b'B', Bunion)
        meta = {
            "ksize": args.ksize, "scaled": args.scaled, "seed": args.seed,
            "n_bins": len(bin_names), "counts": {k: sizes[k] for k in [primary_key,"assembly"]},
            "primary": primary_key,
        }
        with open(dbdir / "meta.json", "w") as w:
            json.dump(meta, w, indent=2)
    else:
        if not HAVE_ROCKS:
            print("[DB] WARNING: python-rocksdb not installed; skipping DB creation.")

    # Compute metrics & report
    M = compute_novel_metrics(sets, bin_names, primary_key=primary_key)
    with open(outdir / "novel_metrics.json", "w") as w:
        json.dump(M, w, indent=2)

    print("[report] building report.html")
    make_report(outdir, [primary_key,"assembly"] + bin_names, sets, sizes,
                args.ksize, args.scaled, args.seed, M, primary_key=primary_key, top_bins=args.top_bins)

    print("\nDone.")
    print(f" - {outdir/'report.html'}")
    print(f" - {outdir/'novel_metrics.json'}")
    if HAVE_ROCKS and not args.skip_db:
        print(f" - {outdir} / hashdb (RocksDB with unions: {primary_key}, assembly, bins_union)")

if __name__ == "__main__":
    main()
