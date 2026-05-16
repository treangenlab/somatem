#!/usr/bin/env python3

"""
Run bounded iterative MAG binning, score each candidate with Pigeon, and select
the best per-binner candidates for optional DAS Tool consensus.

This script is intended to be called by a single Nextflow module when users opt
into iterative binning. It keeps the loop logic in Python, where stopping rules
and per-binner trajectories are easier to express safely than in DSL2 dataflow.

Each iteration reruns each active binner with a seed offset. Pigeon then scores
the candidate bins against the same assembly graph unitig space, so the loop is
optimizing recovery of graph-supported sequence into bins rather than CheckM
marker estimates.
"""

import argparse
import csv
import gzip
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Set, Tuple


FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _as_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def script_dir() -> Path:
    return Path(__file__).resolve().parent


def run_command(cmd: List[str], cwd: Path, log_handle):
    log_handle.write(f"{cwd}\t{' '.join(shlex.quote(str(c)) for c in cmd)}\n")
    log_handle.flush()
    subprocess.run(cmd, cwd=str(cwd), check=True)


def split_args(value: str) -> List[str]:
    return shlex.split(value) if value else []


def symlink_or_copy(source: Path, target: Path, copy: bool = False):
    ensure_dir(target.parent)
    if target.exists() or target.is_symlink():
        target.unlink()
    if copy:
        shutil.copy2(source, target)
    else:
        os.symlink(source.resolve(), target)


def fasta_files(path: Path) -> List[Path]:
    if path.is_file():
        return [path]
    if not path.exists():
        return []
    return sorted(p for p in path.iterdir() if p.is_file() and str(p).endswith(FASTA_SUFFIXES))


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if name is not None:
        yield name, "".join(chunks)


def write_fasta_record(handle, name: str, sequence: str, width: int = 80):
    handle.write(f">{name}\n")
    for i in range(0, len(sequence), width):
        handle.write(sequence[i:i + width] + "\n")


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def canonical_kmers(seq: str, ksize: int) -> Set[str]:
    seq = seq.upper()
    kmers: Set[str] = set()
    for i in range(0, len(seq) - ksize + 1):
        kmer = seq[i:i + ksize]
        if set(kmer) <= {"A", "C", "G", "T"}:
            rc = revcomp(kmer)
            kmers.add(kmer if kmer <= rc else rc)
    return kmers


def fasta_kmer_union(paths: Iterable[Path], ksize: int) -> Set[str]:
    union: Set[str] = set()
    for path in paths:
        for _name, seq in iter_fasta(path):
            union.update(canonical_kmers(seq, ksize))
    return union


def read_json(path: Path) -> Dict[str, Any]:
    with open(path, "r") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return data


def telemetry_from_metrics(path: Path) -> Dict[str, Any]:
    metrics = read_json(path)
    if isinstance(metrics.get("telemetry"), dict):
        t = dict(metrics["telemetry"])
    else:
        t = dict(metrics)

    # Normalize the original pigeon.py novel_metrics.json schema into the
    # small telemetry vocabulary used by the iterative search.
    t.setdefault("pam", metrics.get("pam", 0.0))
    t.setdefault("auc_norm", metrics.get("auc_norm", 0.0))
    t.setdefault("bin_explained_fraction", metrics.get("frac_unitigs_in_B", 0.0))
    t.setdefault("assembly_explained_fraction", metrics.get("frac_unitigs_in_A", 0.0))
    t.setdefault("both_fraction", metrics.get("frac_unitigs_both", 0.0))
    t.setdefault("assembly_only_fraction", metrics.get("frac_unitigs_A_not_B", 0.0))
    t.setdefault("bin_only_fraction", metrics.get("frac_unitigs_B_not_A", 0.0))
    t.setdefault("unexplained_fraction", metrics.get("frac_unitigs_only", 0.0))
    t.setdefault("n_bins", len(metrics.get("ordered_bins", [])))
    t.setdefault("decision", "")
    return t


def score_from_telemetry(t: Dict[str, Any], mode: str) -> float:
    components = score_components_from_telemetry(t)
    if mode == "pam":
        return components["pam"]
    if mode == "bin_explained":
        return components["bin_explained"]
    if mode == "both":
        return components["both"]
    if mode == "unitig_recovery":
        return (
            0.30 * components["bin_explained"] +
            0.25 * components["both"] +
            0.15 * components["explained"] +
            0.10 * components["auc"] +
            0.10 * components["balance"] +
            0.10 * components["residual_resolution"]
        )
    return (
        0.60 * components["pam"] +
        0.25 * components["bin_explained"] +
        0.15 * components["both"]
    )


def score_components_from_telemetry(t: Dict[str, Any]) -> Dict[str, float]:
    balance = 1.0 - abs(
        _as_float(t.get("assembly_only_fraction")) -
        _as_float(t.get("bin_only_fraction"))
    )
    residual_resolution = 1.0 - _as_float(t.get("residual_kmer_fraction", 1.0))
    return {
        "pam": _as_float(t.get("pam")),
        "bin_explained": _as_float(t.get("bin_explained_fraction")),
        "both": _as_float(t.get("both_fraction")),
        "explained": 1.0 - _as_float(t.get("unexplained_fraction")),
        "auc": _as_float(t.get("auc_norm")),
        "balance": max(0.0, balance),
        "residual_resolution": max(0.0, residual_resolution),
    }


def iteration_seed(base_seed: int, iteration: int, seed_step: int) -> int:
    return base_seed + ((iteration - 1) * seed_step)


def run_semibin2(args: argparse.Namespace, iteration: int, target_fasta: Path, outdir: Path, log_handle) -> Path:
    root = outdir / "semibin2"
    cmd = [
        "SemiBin2",
        "single_easy_bin",
        "--input-fasta", str(target_fasta),
        "--input-bam", str(args.bam),
        "--output", str(root),
        "-t", str(args.threads),
        "--random-seed", str(iteration_seed(args.semibin2_seed, iteration, args.seed_step)),
    ]
    if args.sample_environment:
        cmd.extend(["--environment", args.sample_environment])
    else:
        cmd.append("--self-supervised")
    cmd.extend(split_args(args.semibin2_args))
    run_command(cmd, outdir, log_handle)
    return root / "output_bins"


def run_metabat2(args: argparse.Namespace, iteration: int, target_fasta: Path, outdir: Path, log_handle) -> Path:
    bins_dir = outdir / "metabat2_bins"
    ensure_dir(bins_dir)
    depth = outdir / "metabat2.depth.txt"
    run_command([
        "jgi_summarize_bam_contig_depths",
        "--outputDepth", str(depth),
        str(args.bam),
    ], outdir, log_handle)
    cmd = [
        "metabat2",
        "--inFile", str(target_fasta),
        "--abdFile", str(depth),
        "--outFile", str(bins_dir / "metabat2"),
        "--numThreads", str(args.threads),
        "--seed", str(iteration_seed(args.metabat2_seed, iteration, args.seed_step)),
    ]
    cmd.extend(split_args(args.metabat2_args))
    run_command(cmd, outdir, log_handle)
    return bins_dir


def run_vamb(args: argparse.Namespace, iteration: int, target_fasta: Path, outdir: Path, log_handle) -> Path:
    bamdir = outdir / "bam"
    ensure_dir(bamdir)
    symlink_or_copy(args.bam, bamdir / args.bam.name)
    root = outdir / "vamb"
    cmd = [
        "vamb",
        "bin",
        "default",
        "--outdir", str(root),
        "--fasta", str(target_fasta),
        "--bamdir", str(bamdir),
        "-p", str(args.threads),
        "--seed", str(iteration_seed(args.vamb_seed, iteration, args.seed_step)),
        "--minfasta", str(args.vamb_minfasta),
    ]
    cmd.extend(split_args(args.vamb_args))
    run_command(cmd, outdir, log_handle)
    return root / "bins"


def run_pigeon(args: argparse.Namespace,
               binner: str,
               iteration: int,
               bins_dir: Path,
               candidate_dir: Path,
               previous_metrics: List[Path],
               log_handle) -> Path:
    pigeon_out = candidate_dir / "pigeon"
    ensure_dir(pigeon_out)
    cmd = [
        sys.executable,
        str(args.pigeon_script),
        "--gfa", str(args.gfa),
        "--assembly", str(args.assembly),
        "--bins_dir", str(bins_dir),
        "--outdir", str(pigeon_out),
        "--ksize", str(args.pigeon_ksize),
        "--scaled", str(args.pigeon_scaled),
        "--seed", str(args.pigeon_seed),
        "--top-bins", str(args.pigeon_top_bins),
        "--skip-db",
    ]
    run_command(cmd, candidate_dir, log_handle)
    return pigeon_out / "novel_metrics.json"


def write_residual_target(args: argparse.Namespace,
                          source_fasta: Path,
                          bins_dir: Path,
                          out_fasta: Path,
                          summary_path: Path) -> Dict[str, Any]:
    """Write contigs from source_fasta whose k-mers remain poorly explained by bins."""
    ensure_dir(out_fasta.parent)
    bin_kmers = fasta_kmer_union(fasta_files(bins_dir), args.residual_ksize)

    total_contigs = 0
    retained_contigs = 0
    total_bases = 0
    retained_bases = 0
    total_kmers = 0
    retained_kmers = 0

    with open(out_fasta, "w") as writer:
        for name, seq in iter_fasta(source_fasta):
            total_contigs += 1
            total_bases += len(seq)
            kmers = canonical_kmers(seq, args.residual_ksize)
            if not kmers:
                continue
            explained = len(kmers & bin_kmers)
            residual_fraction = 1.0 - (explained / float(len(kmers)))
            total_kmers += len(kmers)
            if (
                len(seq) >= args.residual_min_contig_len and
                residual_fraction >= args.residual_min_unexplained_fraction
            ):
                retained_contigs += 1
                retained_bases += len(seq)
                retained_kmers += len(kmers)
                write_fasta_record(writer, name, seq)

    summary = {
        "source_fasta": str(source_fasta),
        "residual_fasta": str(out_fasta),
        "residual_ksize": args.residual_ksize,
        "total_contigs": total_contigs,
        "retained_contigs": retained_contigs,
        "total_bases": total_bases,
        "retained_bases": retained_bases,
        "total_kmers": total_kmers,
        "retained_kmers": retained_kmers,
        "residual_contig_fraction": retained_contigs / float(total_contigs or 1),
        "residual_base_fraction": retained_bases / float(total_bases or 1),
        "residual_kmer_fraction": retained_kmers / float(total_kmers or 1),
    }
    with open(summary_path, "w") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def write_manifest(path: Path, records: List[Dict[str, Any]]):
    fields = ["sample", "binner", "iteration", "seed", "target_fasta", "residual_fasta", "metrics", "bins"]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow({field: record[field] for field in fields})


def selected_dirs(selection_tsv: Path) -> List[Dict[str, str]]:
    with open(selection_tsv, "r", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_selection_outputs(args: argparse.Namespace, records: List[Dict[str, Any]]) -> Path:
    select_out = args.outdir / "selection"
    selected_root = select_out / "selected_bins" / args.sample_id
    ensure_dir(selected_root)

    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for record in records:
        grouped.setdefault(record["binner"], []).append(record)

    selected: List[Dict[str, Any]] = []
    for binner, rows in grouped.items():
        rows = sorted(rows, key=lambda r: _as_int(r["iteration"]))
        best = max(
            rows,
            key=lambda r: (
                _as_float(r["score"]),
                _as_float(r["telemetry"].get("bin_explained_fraction")),
                -_as_int(r["iteration"]),
            ),
        )
        dest = selected_root / binner
        ensure_dir(dest)
        for fasta in fasta_files(Path(best["bins"])):
            symlink_or_copy(fasta, dest / fasta.name)
        best["selected_bins_dir"] = str(dest)
        best["trajectory_stop_iteration"] = rows[-1]["iteration"]
        best["trajectory_stop_reason"] = rows[-1].get("stop_reason", "max_iterations")
        selected.append(best)

    selected_fields = [
        "sample", "binner", "iteration", "seed", "selection_score", "decision", "pam",
        "bin_explained_fraction", "both_fraction", "unexplained_fraction",
        "residual_contig_fraction", "residual_kmer_fraction",
        "trajectory_stop_iteration", "trajectory_stop_reason", "target_fasta",
        "residual_fasta", "metrics", "bins", "selected_bins_dir",
    ]
    with open(select_out / "selected_iterations.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=selected_fields, delimiter="\t")
        writer.writeheader()
        for row in selected:
            t = row["telemetry"]
            writer.writerow({
                "sample": row["sample"],
                "binner": row["binner"],
                "iteration": row["iteration"],
                "seed": row.get("seed"),
                "selection_score": f"{_as_float(row['score']):.6f}",
                "decision": t.get("decision", ""),
                "pam": t.get("pam"),
                "bin_explained_fraction": t.get("bin_explained_fraction"),
                "both_fraction": t.get("both_fraction"),
                "unexplained_fraction": t.get("unexplained_fraction"),
                "residual_contig_fraction": t.get("residual_contig_fraction"),
                "residual_kmer_fraction": t.get("residual_kmer_fraction"),
                "trajectory_stop_iteration": row.get("trajectory_stop_iteration"),
                "trajectory_stop_reason": row.get("trajectory_stop_reason"),
                "target_fasta": row.get("target_fasta"),
                "residual_fasta": row.get("residual_fasta"),
                "metrics": row["metrics"],
                "bins": row["bins"],
                "selected_bins_dir": row["selected_bins_dir"],
            })

    trajectory_fields = [
        "sample", "binner", "iteration", "seed", "score", "delta_score",
        "best_score_so_far", "no_gain_count", "gain_state", "decision", "pam",
        "bin_explained_fraction", "both_fraction", "unexplained_fraction",
        "residual_contig_fraction", "residual_kmer_fraction", "target_fasta",
        "residual_fasta",
    ]
    with open(select_out / "selection_trajectory.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=trajectory_fields, delimiter="\t")
        writer.writeheader()
        for row in records:
            t = row["telemetry"]
            writer.writerow({
                "sample": row["sample"],
                "binner": row["binner"],
                "iteration": row["iteration"],
                "seed": row.get("seed"),
                "score": row["score"],
                "delta_score": row.get("delta_score"),
                "best_score_so_far": row.get("best_score_so_far"),
                "no_gain_count": row.get("no_gain_count"),
                "gain_state": row.get("gain_state"),
                "decision": t.get("decision", ""),
                "pam": t.get("pam"),
                "bin_explained_fraction": t.get("bin_explained_fraction"),
                "both_fraction": t.get("both_fraction"),
                "unexplained_fraction": t.get("unexplained_fraction"),
                "residual_contig_fraction": t.get("residual_contig_fraction"),
                "residual_kmer_fraction": t.get("residual_kmer_fraction"),
                "target_fasta": row.get("target_fasta"),
                "residual_fasta": row.get("residual_fasta"),
            })

    with open(select_out / "selection_summary.json", "w") as handle:
        json.dump({
            "schema_version": "pigeon.selection.v1",
            "sample_id": args.sample_id,
            "score_mode": args.selection_score_mode,
            "max_iterations": args.max_iterations,
            "min_delta_score": args.selection_min_delta_score,
            "patience": args.selection_patience,
            "residual_rebinning": args.residual_rebinning,
            "residual_ksize": args.residual_ksize,
            "residual_min_unexplained_fraction": args.residual_min_unexplained_fraction,
            "residual_min_contig_len": args.residual_min_contig_len,
            "n_candidates": len(records),
            "n_selected": len(selected),
            "selected": selected,
        }, handle, indent=2)

    return select_out


def run_dastool(args: argparse.Namespace, selection_out: Path, log_handle) -> Path:
    consensus = args.outdir / "consensus"
    ensure_dir(consensus)
    table_paths: List[Path] = []
    labels: List[str] = []

    for row in selected_dirs(selection_out / "selected_iterations.tsv"):
        label = row["binner"]
        bins_dir = Path(row["selected_bins_dir"])
        normalized_bins = consensus / "normalized_bins" / label
        normalize_bins_for_dastool(bins_dir, normalized_bins)
        table = consensus / f"{label}.contigs2bin.tsv"
        cmd = [
            "Fasta_to_Contigs2Bin.sh",
            "-i", str(normalized_bins),
            "-e", "fa",
        ]
        log_handle.write(f"{consensus}\t{' '.join(shlex.quote(str(c)) for c in cmd)} > {table}\n")
        log_handle.flush()
        with open(table, "w") as handle:
            subprocess.run(cmd, cwd=str(consensus), stdout=handle, check=True)
        table_paths.append(table)
        labels.append(label)

    cmd = [
        "DAS_Tool",
        "-i", ",".join(str(p) for p in table_paths),
        "-l", ",".join(labels),
        "-c", str(args.assembly),
        "-o", str(consensus / args.sample_id),
        "--write_bins",
    ]
    cmd.extend(split_args(args.dastool_args))
    run_command(cmd, consensus, log_handle)
    bins = consensus / f"{args.sample_id}_DASTool_bins"
    return bins


def materialize_final_bins(source_dirs: Iterable[Path], final_dir: Path, copy: bool = False):
    ensure_dir(final_dir)
    for source_dir in source_dirs:
        label = source_dir.parent.name if source_dir.parent.name else source_dir.name
        for fasta in fasta_files(source_dir):
            target = final_dir / f"{label}_{fasta.name}"
            symlink_or_copy(fasta, target, copy=copy)


def normalize_bins_for_dastool(source_dir: Path, target_dir: Path):
    ensure_dir(target_dir)
    for fasta in fasta_files(source_dir):
        stem = fasta.name
        if stem.endswith(".gz"):
            stem = stem[:-3]
        stem = stem.rsplit(".", 1)[0]
        target = target_dir / f"{stem}.fa"
        if str(fasta).endswith(".gz"):
            with gzip.open(fasta, "rt") as reader, open(target, "w") as writer:
                shutil.copyfileobj(reader, writer)
        else:
            shutil.copy2(fasta, target)


def run_no_consensus(selection_out: Path, final_dir: Path):
    dirs = [Path(row["selected_bins_dir"]) for row in selected_dirs(selection_out / "selected_iterations.tsv")]
    materialize_final_bins(dirs, final_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run bounded iterative MAG binning and Pigeon-based selection."
    )
    parser.add_argument("--gfa", required=True, type=Path)
    parser.add_argument("--assembly", required=True, type=Path)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--binners", default="semibin2,metabat2,vamb",
                        help="Comma-separated binners to run: semibin2,metabat2,vamb")
    parser.add_argument("--max-iterations", type=int, default=5)
    parser.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1))

    parser.add_argument("--sample-environment", default="")
    parser.add_argument("--seed-step", type=int, default=1,
                        help="Amount added to each binner's base seed per iteration.")
    parser.add_argument("--semibin2-seed", type=int, default=42)
    parser.add_argument("--metabat2-seed", type=int, default=42)
    parser.add_argument("--vamb-seed", type=int, default=42)
    parser.add_argument("--vamb-minfasta", type=int, default=200000)
    parser.add_argument("--semibin2-args", default="")
    parser.add_argument("--metabat2-args", default="")
    parser.add_argument("--vamb-args", default="")

    parser.add_argument("--pigeon-script", type=Path, default=script_dir() / "pigeon.py")
    parser.add_argument("--pigeon-ksize", type=int, default=17)
    parser.add_argument("--pigeon-scaled", type=int, default=1000)
    parser.add_argument("--pigeon-seed", type=int, default=42)
    parser.add_argument("--pigeon-top-bins", type=int, default=20)
    parser.add_argument("--selection-score-mode", default="composite",
                        choices=["composite", "unitig_recovery", "pam", "bin_explained", "both"])
    parser.add_argument("--selection-min-delta-score", type=float, default=0.005)
    parser.add_argument("--selection-patience", type=int, default=2)
    parser.add_argument("--residual-rebinning", action="store_true",
                        help="Use poorly explained residual contigs as the next iteration's binner input.")
    parser.add_argument("--residual-ksize", type=int, default=21,
                        help="Exact canonical k-mer size used for residual contig selection.")
    parser.add_argument("--residual-min-unexplained-fraction", type=float, default=0.50,
                        help="Keep contigs whose k-mers are at least this unexplained by current bins.")
    parser.add_argument("--residual-min-contig-len", type=int, default=1500,
                        help="Minimum contig length retained in residual FASTA.")
    parser.add_argument("--residual-min-contigs", type=int, default=10,
                        help="Stop residual rebinning for a binner if fewer contigs remain.")

    parser.add_argument("--consensus-tool", default="dastool", choices=["dastool", "none"])
    parser.add_argument("--dastool-args", default="--search_engine diamond")
    args = parser.parse_args()

    args.gfa = args.gfa.resolve()
    args.assembly = args.assembly.resolve()
    args.bam = args.bam.resolve()
    args.outdir = args.outdir.resolve()
    args.pigeon_script = args.pigeon_script.resolve()
    args.max_iterations = max(1, min(args.max_iterations, 5))
    return args


def main():
    args = parse_args()
    ensure_dir(args.outdir)
    candidates_root = args.outdir / "candidates"
    ensure_dir(candidates_root)
    final_bins = args.outdir / "final_bins"

    binners = [b.strip() for b in args.binners.split(",") if b.strip()]
    runners = {
        "semibin2": run_semibin2,
        "metabat2": run_metabat2,
        "vamb": run_vamb,
    }
    unknown = sorted(set(binners) - set(runners))
    if unknown:
        raise SystemExit(f"Unknown binner(s): {', '.join(unknown)}")

    records: List[Dict[str, Any]] = []
    metrics_by_binner: Dict[str, List[Path]] = {binner: [] for binner in binners}
    best_score_by_binner: Dict[str, Optional[float]] = {binner: None for binner in binners}
    no_gain_by_binner: Dict[str, int] = {binner: 0 for binner in binners}
    target_fasta_by_binner: Dict[str, Path] = {binner: args.assembly for binner in binners}
    active = set(binners)

    command_log = args.outdir / "command_log.tsv"
    with open(command_log, "w") as log_handle:
        log_handle.write("cwd\tcommand\n")
        for iteration in range(1, args.max_iterations + 1):
            if not active:
                break
            for binner in list(binners):
                if binner not in active:
                    continue
                candidate_dir = candidates_root / binner / f"iter_{iteration:02d}"
                ensure_dir(candidate_dir)
                target_fasta = target_fasta_by_binner[binner]
                bins_dir = runners[binner](args, iteration, target_fasta, candidate_dir, log_handle)
                if not fasta_files(bins_dir):
                    active.remove(binner)
                    continue
                metrics = run_pigeon(
                    args,
                    binner,
                    iteration,
                    bins_dir,
                    candidate_dir,
                    metrics_by_binner[binner],
                    log_handle,
                )
                metrics_by_binner[binner].append(metrics)
                t = telemetry_from_metrics(metrics)
                residual_summary: Dict[str, Any] = {}
                residual_fasta = candidate_dir / "residual" / f"{binner}.iter_{iteration:02d}.residual.fa"
                if args.residual_rebinning:
                    residual_summary = write_residual_target(
                        args,
                        target_fasta,
                        bins_dir,
                        residual_fasta,
                        candidate_dir / "residual" / "residual_summary.json",
                    )
                    t.update({
                        "residual_contig_fraction": residual_summary["residual_contig_fraction"],
                        "residual_base_fraction": residual_summary["residual_base_fraction"],
                        "residual_kmer_fraction": residual_summary["residual_kmer_fraction"],
                    })
                else:
                    t.update({
                        "residual_contig_fraction": 1.0,
                        "residual_base_fraction": 1.0,
                        "residual_kmer_fraction": 1.0,
                    })

                score = score_from_telemetry(t, args.selection_score_mode)
                score_components = score_components_from_telemetry(t)
                previous_best = best_score_by_binner[binner]
                if previous_best is None:
                    delta_score = None
                    best_score_by_binner[binner] = score
                    no_gain_by_binner[binner] = 0
                    gain_state = "improved"
                else:
                    delta_score = score - previous_best
                    if score > previous_best + args.selection_min_delta_score:
                        best_score_by_binner[binner] = score
                        no_gain_by_binner[binner] = 0
                        gain_state = "improved"
                    else:
                        no_gain_by_binner[binner] += 1
                        gain_state = "plateau"

                stop_reason = ""
                if no_gain_by_binner[binner] >= args.selection_patience:
                    stop_reason = "score_plateau"
                    active.remove(binner)
                elif (
                    args.residual_rebinning and
                    residual_summary and
                    residual_summary["retained_contigs"] < args.residual_min_contigs
                ):
                    stop_reason = "residual_exhausted"
                    active.remove(binner)
                elif iteration == args.max_iterations:
                    stop_reason = "max_iterations"
                elif t.get("decision") == "stop_candidate":
                    stop_reason = "pigeon_stop_candidate"
                    active.remove(binner)
                elif args.residual_rebinning and residual_summary:
                    target_fasta_by_binner[binner] = residual_fasta

                records.append({
                    "sample": args.sample_id,
                    "binner": binner,
                    "iteration": iteration,
                    "seed": iteration_seed(
                        getattr(args, f"{binner}_seed"),
                        iteration,
                        args.seed_step,
                    ),
                    "target_fasta": str(target_fasta),
                    "residual_fasta": str(residual_fasta) if args.residual_rebinning else "",
                    "metrics": str(metrics),
                    "bins": str(bins_dir),
                    "telemetry": t,
                    "score": score,
                    "score_components": score_components,
                    "delta_score": delta_score,
                    "best_score_so_far": best_score_by_binner[binner],
                    "no_gain_count": no_gain_by_binner[binner],
                    "gain_state": gain_state,
                    "stop_reason": stop_reason,
                })

        manifest = args.outdir / "candidate_manifest.tsv"
        write_manifest(manifest, records)
        selection_out = write_selection_outputs(args, records)

        if args.consensus_tool == "dastool":
            consensus_bins = run_dastool(args, selection_out, log_handle)
            materialize_final_bins([consensus_bins], final_bins)
        else:
            run_no_consensus(selection_out, final_bins)

    summary = {
        "schema_version": "pigeon.iterative_binning.v1",
        "sample_id": args.sample_id,
        "binners": binners,
        "max_iterations": args.max_iterations,
        "n_candidates": len(records),
        "consensus_tool": args.consensus_tool,
        "score_mode": args.selection_score_mode,
        "residual_rebinning": args.residual_rebinning,
        "residual_ksize": args.residual_ksize,
        "residual_min_unexplained_fraction": args.residual_min_unexplained_fraction,
        "residual_min_contig_len": args.residual_min_contig_len,
        "final_bins": str(final_bins),
        "candidate_manifest": str(args.outdir / "candidate_manifest.tsv"),
        "selection_summary": str(args.outdir / "selection" / "selection_summary.json"),
    }
    with open(args.outdir / "iterative_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)


if __name__ == "__main__":
    main()
