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
    * markdown (optional, for rendering plot guide in HTML report)
"""

import re
import os, sys, json, argparse, subprocess, gzip, struct, shutil, math
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


PIGEON_LOSS_MODEL = "adaptive_focus_v1"
PIGEON_LOSS_COMPONENTS = [
    "missing_bin_recovery",
    "unexplained_unitigs",
    "assembly_bin_discordance",
    "fragmented_recovery",
    "support_imbalance",
    "novelty_partition_distance",
]


# ---------------- I/O helpers ----------------
def ensure_dir(p: Path):
    """Create directory and any missing parent directories."""
    p.mkdir(parents=True, exist_ok=True)

def open_maybe_gzip(path, mode="rt"):
    """Open file, automatically handling gzip compression based on extension."""
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)

def safe_name(s: str) -> str:
    """Sanitize string for use as filename (alphanumeric, dash, underscore, dot only)."""
    allowed = set("-_.")
    return "".join(ch if ch.isalnum() or ch in allowed else "_" for ch in s)[:200]

def format_name(s: str) -> str:
    """If the filename contains 'contig_' prefix it with 'bin_', otherwise keep the file basename (for assembly/unitig)."""
    p = Path(s)
    name = p.name
    if "contig_" in name:
        name = "bin_" + name
    return safe_name(name)

def display_name(s: str) -> str:
    """Shorten bin labels for display. E.g., 'bin_contig_1.fa.gz' -> '1'."""
    p = Path(s)

    # replace either "contig_" or "bin_contig_" with ""
    name = re.sub(r'(?:bin_)?contig_', '', p.name)
    return name


def require_path(path: str, label: str) -> Path:
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"{label} does not exist: {p}")
    return p


def fasta_inputs(path: Path) -> List[str]:
    return [
        str(f) for f in sorted(path.iterdir())
        if f.is_file() and any(str(f).endswith(ext) for ext in (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"))
    ]


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


def clamp01(value: float) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return 0.0
    if math.isnan(x) or math.isinf(x):
        return 0.0
    return max(0.0, min(1.0, x))


def adaptive_focus_loss(components: Dict[str, float]) -> Tuple[float, Dict[str, float], float]:
    """Return an adaptive bounded loss, focus weights, and temperature.

    This is intentionally not a static weighted sum. Each component is a
    normalized deficit where 0 is ideal and 1 is poor. The candidate's own
    deficit spread controls the softmax temperature, and the softmax-derived
    focus weights make the largest current failure mode lead the objective.
    """
    clean = {key: clamp01(value) for key, value in components.items()}
    if not clean:
        return 0.0, {}, 0.0

    values = list(clean.values())
    mean = sum(values) / float(len(values))
    variance = sum((value - mean) ** 2 for value in values) / float(len(values))
    temperature = 2.0 + (12.0 * math.sqrt(variance))

    max_value = max(values)
    scaled = {
        key: math.exp(temperature * (value - max_value))
        for key, value in clean.items()
    }
    total = sum(scaled.values()) or 1.0
    focus = {key: value / total for key, value in scaled.items()}
    loss = sum(focus[key] * clean[key] for key in clean)
    return clamp01(loss), focus, temperature


def dominant_loss_component(components: Dict[str, float], focus: Dict[str, float]) -> str:
    if not components:
        return ""
    return max(
        components,
        key=lambda key: (clamp01(components.get(key, 0.0)) * clamp01(focus.get(key, 0.0)), clamp01(components.get(key, 0.0))),
    )


def loss_profile_rows(metrics: Dict) -> List[Dict[str, float]]:
    components = metrics.get("pigeon_loss_components", {}) or {}
    focus = metrics.get("pigeon_loss_focus", {}) or {}
    rows = []
    for name in PIGEON_LOSS_COMPONENTS:
        value = clamp01(components.get(name, 0.0))
        weight = clamp01(focus.get(name, 0.0))
        rows.append({
            "component": name,
            "defect": value,
            "focus": weight,
            "contribution": value * weight,
        })
    rows.sort(key=lambda row: row["contribution"], reverse=True)
    return rows


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
    if not Uset:
        raise ValueError(f"Primary signature has no hashes: {primary_key}")
    if not A:
        raise ValueError("Assembly signature has no hashes")
    if not bin_names:
        raise ValueError("No bin signatures were provided")
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

    # Pigeon novelty loss: a bounded adaptive objective for iterative binning.
    # Lower is better. Components are normalized deficits; the final loss is an
    # adaptive focus mean, so the dominant current failure mode leads the search
    # without fixed component weights.
    novelty_partition_distance = math.sqrt(clamp01(1.0 - math.sqrt(clamp01(frac_U_both))))
    loss_components = {
        "missing_bin_recovery": clamp01(1.0 - frac_U_in_B),
        "unexplained_unitigs": clamp01(frac_U_only),
        "assembly_bin_discordance": clamp01(frac_U_Aonly + frac_U_Bonly),
        "fragmented_recovery": clamp01(1.0 - auc_norm),
        "support_imbalance": clamp01(abs(frac_U_Aonly - frac_U_Bonly)),
        "novelty_partition_distance": novelty_partition_distance,
    }
    pigeon_loss, loss_focus, loss_temperature = adaptive_focus_loss(loss_components)
    dominant_component = dominant_loss_component(loss_components, loss_focus)

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
        "pigeon_loss": pigeon_loss,
        "pigeon_loss_score": 1.0 - pigeon_loss,
        "pigeon_loss_components": loss_components,
        "pigeon_loss_focus": loss_focus,
        "pigeon_loss_temperature": loss_temperature,
        "pigeon_loss_model": PIGEON_LOSS_MODEL,
        "pigeon_loss_component_order": PIGEON_LOSS_COMPONENTS,
        "pigeon_loss_dominant_component": dominant_component,
        "pigeon_loss_profile": loss_profile_rows({
            "pigeon_loss_components": loss_components,
            "pigeon_loss_focus": loss_focus,
        }),
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
    labels = [display_name(n) for n in names]

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
                             text=[display_name(n) for n in names],
                             textposition="top center"), row=1, col=2)

    # Top bins bar
    top = novel_metrics.get("per_bin_explained", [])[:top_bins]
    fig.add_trace(go.Bar(x=[display_name(b) for b, _ in top], y=[v for _, v in top], showlegend=False), row=1, col=3)
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
        "AUC (cum curve)", "PAM score", "Pigeon loss", "Loss score"
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
        f"{novel_metrics.get('pigeon_loss', 0.0):.3f}",
        f"{novel_metrics.get('pigeon_loss_score', 1.0):.3f}",
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
        xs = [display_name(b) for b, _ in top]
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
    
    # Generate HTML with plot summaries appended
    html_content = fig.to_html(include_plotlyjs="cdn")
    
    # Load and convert markdown guide to HTML
    guide_md_path = Path(__file__).parent / "pigeon_report_guide.md"
    if guide_md_path.exists():
        try:
            import markdown
            with open(guide_md_path, "r") as f:
                md_content = f.read()
            guide_html = markdown.markdown(md_content, extensions=['extra', 'nl2br'])
            # Wrap in styled div
            summary_html = f"""
<div style="max-width: 1400px; margin: 40px auto; padding: 20px; font-family: Arial, sans-serif;">
    <style>
        .guide-content h1 {{ border-bottom: 2px solid #333; padding-bottom: 10px; }}
        .guide-content h2 {{ color: #2c5aa0; margin-top: 20px; }}
        .guide-content ul {{ line-height: 1.8; }}
        .guide-content hr + h2 {{ margin-top: 30px; }}
        .guide-content p {{ background-color: #f0f0f0; padding: 15px; border-left: 4px solid #2c5aa0; }}
    </style>
    <div class="guide-content">
{guide_html}
    </div>
</div>
"""
        except ImportError:
            # Fallback if markdown module not available
            summary_html = f"<!-- Markdown conversion requires 'markdown' package. Install with: pip install markdown -->"
    else:
        summary_html = "<!-- Plot guide markdown file not found -->"
    
    # Insert summary before closing </body> tag
    html_with_summary = html_content.replace("</body>", summary_html + "</body>")
    
    with open(outdir / "report.html", "w") as f:
        f.write(html_with_summary)


def write_loss_profile(outdir: Path, metrics: Dict):
    rows = loss_profile_rows(metrics)
    with open(outdir / "loss_profile.tsv", "w") as handle:
        handle.write("component\tdefect\tfocus\tcontribution\n")
        for row in rows:
            handle.write(
                f"{row['component']}\t{row['defect']:.8f}\t{row['focus']:.8f}\t{row['contribution']:.8f}\n"
            )
    with open(outdir / "loss_profile.json", "w") as handle:
        json.dump({
            "schema_version": "pigeon.loss_profile.v1",
            "model": metrics.get("pigeon_loss_model", PIGEON_LOSS_MODEL),
            "loss": metrics.get("pigeon_loss"),
            "loss_score": metrics.get("pigeon_loss_score"),
            "dominant_component": metrics.get("pigeon_loss_dominant_component"),
            "temperature": metrics.get("pigeon_loss_temperature"),
            "components": rows,
        }, handle, indent=2)


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
    if args.ksize <= 0:
        raise SystemExit("--ksize must be a positive integer")
    if args.scaled <= 0:
        raise SystemExit("--scaled must be a positive integer")
    if args.top_bins <= 0:
        raise SystemExit("--top-bins must be a positive integer")
    require_path(args.assembly, "--assembly")
    if args.gfa:
        require_path(args.gfa, "--gfa")
    if args.unitigs_fa:
        require_path(args.unitigs_fa, "--unitigs-fa")
    if args.unitigs_sig:
        require_path(args.unitigs_sig, "--unitigs-sig")
    require_path(args.bins_dir, "--bins_dir")

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
        bin_files = fasta_inputs(Path(args.bins_dir))
        if not bin_files:
            raise SystemExit(f"No bin FASTAs found in --bins_dir: {args.bins_dir}")
        sketch_many(bin_files, [args.ksize], args.scaled, args.seed, bins_sig)

    # Load signatures -> sets
    primary_key = "unitigs"
    names = [primary_key, "assembly"]
    sets: Dict[str, Set[int]] = {}
    sizes: Dict[str, int] = {}

    unitig_sigs = load_sig_zip(str(unitig_sig), ksize=args.ksize)
    if not unitig_sigs:
        raise SystemExit(f"No unitig signatures found for k={args.ksize}: {unitig_sig}")
    sig_unitigs = unitig_sigs[0]
    sets[primary_key] = get_hash_set(sig_unitigs); sizes[primary_key] = len(sets[primary_key])

    asm_sigs = load_sig_zip(str(asm_sig), ksize=args.ksize)
    if not asm_sigs:
        raise SystemExit(f"No assembly signatures found for k={args.ksize}: {asm_sig}")
    sig_asm = asm_sigs[0]
    sets["assembly"] = get_hash_set(sig_asm); sizes["assembly"] = len(sets["assembly"])

    # Load signatures for bins
    bin_sigs = load_sig_zip(str(bins_sig), ksize=args.ksize)
    if not bin_sigs:
        raise SystemExit(f"No bin signatures found for k={args.ksize}: {bins_sig}")
    bin_names: List[str] = []
    for sig in bin_sigs:
        nm = sig.name or sig.filename or "bin"
        nm = format_name(nm)
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
    write_loss_profile(outdir, M)

    print("[report] building report.html")
    make_report(outdir, [primary_key,"assembly"] + bin_names, sets, sizes,
                args.ksize, args.scaled, args.seed, M, primary_key=primary_key, top_bins=args.top_bins)

    print("\nDone.")
    print(f" - {outdir/'report.html'}")
    print(f" - {outdir/'novel_metrics.json'}")
    print(f" - {outdir/'loss_profile.tsv'}")
    if HAVE_ROCKS and not args.skip_db:
        print(f" - {outdir} / hashdb (RocksDB with unions: {primary_key}, assembly, bins_union)")

if __name__ == "__main__":
    main()
