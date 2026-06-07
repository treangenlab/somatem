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
import html
import json
import math
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Set, Tuple

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:  # graceful standalone fallback; Nextflow env includes plotly
    go = None
    make_subplots = None


FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
DEFAULT_BINNER_ORDER = ("semibin2", "metabat2", "vamb")
BAM_REFERENCE_MATCH_BINNER = {"vamb"}
PIGEON_LOSS_MODEL = "adaptive_focus_v1"
PIGEON_LOSS_COMPONENTS = [
    "missing_bin_recovery",
    "unexplained_unitigs",
    "assembly_bin_discordance",
    "fragmented_recovery",
    "support_imbalance",
    "novelty_partition_distance",
]


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


def _clamp01(value: Any) -> float:
    x = _as_float(value)
    if math.isnan(x) or math.isinf(x):
        return 0.0
    return max(0.0, min(1.0, x))


def script_dir() -> Path:
    return Path(__file__).resolve().parent


def run_command(cmd: List[str], cwd: Path, log_handle):
    log_handle.write(f"{cwd}\t{' '.join(shlex.quote(str(c)) for c in cmd)}\n")
    log_handle.flush()
    subprocess.run(cmd, cwd=str(cwd), check=True)


def log_note(log_handle, cwd: Path, message: str):
    log_handle.write(f"{cwd}\t# {message}\n")
    log_handle.flush()


def mark_binner_stopped(
    records: List[Dict[str, Any]],
    active: Set[str],
    binner: str,
    reason: str,
):
    active.discard(binner)
    for record in reversed(records):
        if record["binner"] == binner:
            record["stop_reason"] = reason
            return


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
    if "pigeon_loss" not in t:
        loss, focus, temperature = pigeon_loss_from_telemetry(t)
        t["pigeon_loss"] = metrics.get("pigeon_loss", loss)
        t["pigeon_loss_score"] = metrics.get("pigeon_loss_score", 1.0 - _as_float(t.get("pigeon_loss")))
        t["pigeon_loss_components"] = metrics.get("pigeon_loss_components", pigeon_loss_components_from_telemetry(t))
        t["pigeon_loss_focus"] = metrics.get("pigeon_loss_focus", focus)
        t["pigeon_loss_temperature"] = metrics.get("pigeon_loss_temperature", temperature)
        t["pigeon_loss_model"] = metrics.get("pigeon_loss_model", PIGEON_LOSS_MODEL)
        t["pigeon_loss_dominant_component"] = metrics.get("pigeon_loss_dominant_component", dominant_loss_component(t["pigeon_loss_components"], t["pigeon_loss_focus"]))
    else:
        t.setdefault("pigeon_loss_score", metrics.get("pigeon_loss_score", 1.0 - _as_float(t.get("pigeon_loss"))))
        t.setdefault("pigeon_loss_components", metrics.get("pigeon_loss_components", pigeon_loss_components_from_telemetry(t)))
        loss, focus, temperature = pigeon_loss_from_telemetry(t)
        t.setdefault("pigeon_loss_focus", metrics.get("pigeon_loss_focus", focus))
        t.setdefault("pigeon_loss_temperature", metrics.get("pigeon_loss_temperature", temperature))
        t.setdefault("pigeon_loss_model", metrics.get("pigeon_loss_model", PIGEON_LOSS_MODEL))
        t.setdefault("pigeon_loss_dominant_component", metrics.get("pigeon_loss_dominant_component", dominant_loss_component(t.get("pigeon_loss_components", {}), t.get("pigeon_loss_focus", {}))))
    return t


def pigeon_loss_components_from_telemetry(t: Dict[str, Any]) -> Dict[str, float]:
    both_fraction = _clamp01(t.get("both_fraction"))
    novelty_partition_distance = math.sqrt(_clamp01(1.0 - math.sqrt(both_fraction)))
    components = {
        "missing_bin_recovery": _clamp01(1.0 - _as_float(t.get("bin_explained_fraction"))),
        "unexplained_unitigs": _clamp01(t.get("unexplained_fraction")),
        "assembly_bin_discordance": _clamp01(
            _as_float(t.get("assembly_only_fraction")) +
            _as_float(t.get("bin_only_fraction"))
        ),
        "fragmented_recovery": _clamp01(1.0 - _as_float(t.get("auc_norm"))),
        "support_imbalance": _clamp01(
            abs(
                _as_float(t.get("assembly_only_fraction")) -
                _as_float(t.get("bin_only_fraction"))
            )
        ),
        "novelty_partition_distance": novelty_partition_distance,
    }
    return components


def adaptive_focus_loss(components: Dict[str, float]) -> Tuple[float, Dict[str, float], float]:
    clean = {key: _clamp01(value) for key, value in components.items()}
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
    return _clamp01(loss), focus, temperature


def pigeon_loss_from_telemetry(t: Dict[str, Any]) -> Tuple[float, Dict[str, float], float]:
    return adaptive_focus_loss(pigeon_loss_components_from_telemetry(t))


def dominant_loss_component(components: Dict[str, float], focus: Dict[str, float]) -> str:
    if not components:
        return ""
    return max(
        components,
        key=lambda key: (_clamp01(components.get(key, 0.0)) * _clamp01(focus.get(key, 0.0)), _clamp01(components.get(key, 0.0))),
    )


def score_from_telemetry(t: Dict[str, Any], mode: str) -> float:
    components = score_components_from_telemetry(t)
    if mode == "pam":
        return components["pam"]
    if mode == "bin_explained":
        return components["bin_explained"]
    if mode == "both":
        return components["both"]
    if mode in {"loss", "pigeon_loss"}:
        return components["pigeon_loss_score"]
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
        "pigeon_loss": _clamp01(t.get("pigeon_loss", pigeon_loss_from_telemetry(t)[0])),
        "pigeon_loss_score": _clamp01(t.get("pigeon_loss_score", 1.0 - pigeon_loss_from_telemetry(t)[0])),
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
        "pigeon_loss", "pigeon_loss_score", "pigeon_loss_dominant_component",
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
                "pigeon_loss": t.get("pigeon_loss"),
                "pigeon_loss_score": t.get("pigeon_loss_score"),
                "pigeon_loss_dominant_component": t.get("pigeon_loss_dominant_component"),
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
        "pigeon_loss", "pigeon_loss_score", "pigeon_loss_dominant_component",
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
                "pigeon_loss": t.get("pigeon_loss"),
                "pigeon_loss_score": t.get("pigeon_loss_score"),
                "pigeon_loss_dominant_component": t.get("pigeon_loss_dominant_component"),
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
            "pigeon_loss_model": PIGEON_LOSS_MODEL,
            "pigeon_loss_components": PIGEON_LOSS_COMPONENTS,
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

    selected_rows = selected_dirs(selection_out / "selected_iterations.tsv")
    if not selected_rows:
        raise SystemExit("DAS Tool consensus requested, but no selected binner outputs were available.")

    for row in selected_rows:
        label = row["binner"]
        bins_dir = Path(row["selected_bins_dir"])
        normalized_bins = consensus / "normalized_bins" / label
        if not fasta_files(bins_dir):
            raise SystemExit(f"Selected bins directory for {label} is empty: {bins_dir}")
        normalize_bins_for_dastool(bins_dir, normalized_bins)
        table = consensus / f"{label}.contigs2bin.tsv"
        n_contigs = write_contigs2bin_table(normalized_bins, table)
        log_note(
            log_handle,
            consensus,
            f"wrote {n_contigs} contig-to-bin rows for {label} to {table}",
        )
        if n_contigs == 0:
            raise SystemExit(f"Selected bins directory for {label} had no FASTA records: {bins_dir}")
        table_paths.append(table)
        labels.append(label)

    if not table_paths:
        raise SystemExit("DAS Tool consensus requested, but no contig-to-bin tables were created.")

    assembly_for_dastool = uncompressed_fasta_for_dastool(
        args.assembly,
        consensus / f"{args.sample_id}.assembly.fa",
    )
    if assembly_for_dastool != args.assembly:
        log_note(log_handle, consensus, f"decompressed assembly for DAS Tool to {assembly_for_dastool}")

    cmd = [
        "DAS_Tool",
        "-i", ",".join(str(p) for p in table_paths),
        "-l", ",".join(labels),
        "-c", str(assembly_for_dastool),
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


def bin_label_from_fasta(path: Path) -> str:
    stem = path.name
    if stem.endswith(".gz"):
        stem = stem[:-3]
    return stem.rsplit(".", 1)[0]


def normalize_bins_for_dastool(source_dir: Path, target_dir: Path):
    ensure_dir(target_dir)
    for fasta in fasta_files(source_dir):
        target = target_dir / f"{bin_label_from_fasta(fasta)}.fa"
        if str(fasta).endswith(".gz"):
            with gzip.open(fasta, "rt") as reader, open(target, "w") as writer:
                shutil.copyfileobj(reader, writer)
        else:
            shutil.copy2(fasta, target)


def write_contigs2bin_table(bins_dir: Path, table: Path) -> int:
    n_contigs = 0
    with open(table, "w") as handle:
        for fasta in fasta_files(bins_dir):
            bin_label = bin_label_from_fasta(fasta)
            for contig, _sequence in iter_fasta(fasta):
                handle.write(f"{contig}\t{bin_label}\n")
                n_contigs += 1
    return n_contigs


def uncompressed_fasta_for_dastool(source: Path, target: Path) -> Path:
    if not str(source).endswith(".gz"):
        return source
    ensure_dir(target.parent)
    with gzip.open(source, "rt") as reader, open(target, "w") as writer:
        shutil.copyfileobj(reader, writer)
    return target


def run_no_consensus(selection_out: Path, final_dir: Path):
    dirs = [Path(row["selected_bins_dir"]) for row in selected_dirs(selection_out / "selected_iterations.tsv")]
    materialize_final_bins(dirs, final_dir)


def summarize_binner_convergence(records: List[Dict[str, Any]], binners: List[str]) -> Dict[str, Dict[str, Any]]:
    summary: Dict[str, Dict[str, Any]] = {}
    convergence_reasons = {"score_plateau", "residual_exhausted", "pigeon_stop_candidate"}

    for binner in binners:
        rows = sorted(
            [record for record in records if record["binner"] == binner],
            key=lambda record: _as_int(record.get("iteration")),
        )
        if not rows:
            summary[binner] = {
                "status": "no_bins",
                "converged": False,
                "iterations_run": 0,
                "optimal_iteration": None,
                "optimal_score": None,
                "optimal_bins": "",
                "selected_bins_dir": "",
            }
            continue

        best = max(
            rows,
            key=lambda record: (
                _as_float(record.get("score")),
                _as_float(record.get("telemetry", {}).get("bin_explained_fraction")),
                -_as_int(record.get("iteration")),
            ),
        )
        final = rows[-1]
        stop_reason = final.get("stop_reason") or "max_iterations"
        summary[binner] = {
            "status": stop_reason,
            "converged": stop_reason in convergence_reasons,
            "iterations_run": len(rows),
            "optimal_iteration": _as_int(best.get("iteration")),
            "optimal_score": _as_float(best.get("score")),
            "optimal_bins": best.get("bins", ""),
            "selected_bins_dir": best.get("selected_bins_dir", ""),
            "final_score": _as_float(final.get("score")),
        }

    return summary


def _fmt(value: Any, digits: int = 4) -> str:
    if value is None or value == "":
        return ""
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return html.escape(str(value))


def _html_table(headers: List[str], rows: List[List[Any]]) -> str:
    head = "".join(f"<th>{html.escape(str(header))}</th>" for header in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def write_iterative_report(args: argparse.Namespace,
                           records: List[Dict[str, Any]],
                           summary: Dict[str, Any],
                           binner_convergence: Dict[str, Dict[str, Any]]):
    report_path = args.outdir / "iterative_report.html"
    if not records:
        report_path.write_text("<html><body><h1>Pigeon Iterative Binning</h1><p>No candidate records were produced.</p></body></html>")
        return report_path

    ordered = sorted(records, key=lambda row: (row["binner"], _as_int(row["iteration"])))
    binners = sorted({row["binner"] for row in ordered})
    fig = None
    if make_subplots is not None and go is not None:
        fig = make_subplots(
            rows=2,
            cols=2,
            subplot_titles=(
                "Optimization Score by Iteration",
                "Pigeon Loss by Iteration",
                "Graph Recovery by Iteration",
                "Dominant Loss Component Counts",
            ),
        )

    selected_keys = {
        (row["binner"], _as_int(row["iteration"]))
        for row in ordered
        if row.get("selected_bins_dir")
    }

    for binner in binners:
        rows = [row for row in ordered if row["binner"] == binner]
        x = [_as_int(row["iteration"]) for row in rows]
        score = [_as_float(row.get("score")) for row in rows]
        loss = [_as_float(row.get("telemetry", {}).get("pigeon_loss")) for row in rows]
        pam = [_as_float(row.get("telemetry", {}).get("pam")) for row in rows]
        recovery = [_as_float(row.get("telemetry", {}).get("bin_explained_fraction")) for row in rows]
        selected_x = [_as_int(row["iteration"]) for row in rows if (row["binner"], _as_int(row["iteration"])) in selected_keys]
        selected_score = [_as_float(row.get("score")) for row in rows if (row["binner"], _as_int(row["iteration"])) in selected_keys]

        if fig is not None:
            fig.add_trace(go.Scatter(x=x, y=score, mode="lines+markers", name=f"{binner} score"), row=1, col=1)
            if selected_x:
                fig.add_trace(go.Scatter(x=selected_x, y=selected_score, mode="markers", marker=dict(size=13, symbol="star"), name=f"{binner} selected"), row=1, col=1)
            fig.add_trace(go.Scatter(x=x, y=loss, mode="lines+markers", name=f"{binner} loss"), row=1, col=2)
            fig.add_trace(go.Scatter(x=x, y=recovery, mode="lines+markers", name=f"{binner} bin recovery"), row=2, col=1)
            fig.add_trace(go.Scatter(x=x, y=pam, mode="lines+markers", name=f"{binner} PAM", line=dict(dash="dot")), row=2, col=1)

    component_counts: Dict[str, int] = {}
    for row in ordered:
        component = row.get("telemetry", {}).get("pigeon_loss_dominant_component") or "unknown"
        component_counts[component] = component_counts.get(component, 0) + 1
    if fig is not None:
        fig.add_trace(
            go.Bar(x=list(component_counts.keys()), y=list(component_counts.values()), name="dominant component"),
            row=2,
            col=2,
        )

        fig.update_xaxes(title_text="Iteration", row=1, col=1)
        fig.update_xaxes(title_text="Iteration", row=1, col=2)
        fig.update_xaxes(title_text="Iteration", row=2, col=1)
        fig.update_yaxes(title_text="Score", row=1, col=1)
        fig.update_yaxes(title_text="Loss", row=1, col=2)
        fig.update_yaxes(title_text="Fraction / score", row=2, col=1)
        fig.update_yaxes(title_text="Candidate count", row=2, col=2)
        fig.update_layout(template="plotly_white", height=850, title=f"Pigeon Iterative Binning Report: {args.sample_id}")

    selected_rows = []
    for binner in binners:
        info = binner_convergence.get(binner, {})
        best = next(
            (row for row in ordered if row["binner"] == binner and _as_int(row["iteration"]) == _as_int(info.get("optimal_iteration"))),
            None,
        )
        telemetry = best.get("telemetry", {}) if best else {}
        selected_rows.append([
            binner,
            info.get("optimal_iteration", ""),
            _fmt(info.get("optimal_score")),
            _fmt(telemetry.get("pigeon_loss")),
            _fmt(telemetry.get("pam")),
            _fmt(telemetry.get("bin_explained_fraction")),
            telemetry.get("pigeon_loss_dominant_component", ""),
            info.get("status", ""),
        ])

    candidate_rows = []
    for row in ordered:
        telemetry = row.get("telemetry", {})
        candidate_rows.append([
            row.get("binner", ""),
            row.get("iteration", ""),
            row.get("seed", ""),
            _fmt(row.get("score")),
            _fmt(row.get("delta_score")),
            _fmt(telemetry.get("pigeon_loss")),
            _fmt(telemetry.get("pam")),
            _fmt(telemetry.get("bin_explained_fraction")),
            telemetry.get("pigeon_loss_dominant_component", ""),
            row.get("gain_state", ""),
            row.get("stop_reason", ""),
        ])

    css = """
    body { font-family: Arial, sans-serif; margin: 28px; color: #222; }
    h1, h2 { color: #1f2933; }
    .summary { display: flex; flex-wrap: wrap; gap: 12px; margin: 18px 0; }
    .metric { border: 1px solid #d5dce3; border-radius: 6px; padding: 10px 14px; background: #f8fafc; min-width: 150px; }
    .metric strong { display: block; font-size: 18px; margin-top: 4px; }
    table { border-collapse: collapse; width: 100%; margin: 14px 0 28px; font-size: 13px; }
    th, td { border: 1px solid #d7dde4; padding: 6px 8px; text-align: left; }
    th { background: #eef3f8; }
    code { background: #eef3f8; padding: 1px 4px; border-radius: 4px; }
    """
    html_doc = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Pigeon Iterative Binning Report - {html.escape(args.sample_id)}</title>
<style>{css}</style>
</head>
<body>
<h1>Pigeon Iterative Binning Report: {html.escape(args.sample_id)}</h1>
<div class="summary">
  <div class="metric">Candidates<strong>{len(records)}</strong></div>
  <div class="metric">Binners<strong>{len(binners)}</strong></div>
  <div class="metric">Score mode<strong>{html.escape(str(args.selection_score_mode))}</strong></div>
  <div class="metric">Loss model<strong>{html.escape(PIGEON_LOSS_MODEL)}</strong></div>
  <div class="metric">Consensus<strong>{html.escape(str(args.consensus_tool))}</strong></div>
</div>
<p>This report shows how Pigeon selected candidate bins across binners and seed iterations. Higher optimization score is better; lower Pigeon loss is better. Star markers indicate the candidate selected for each binner before consensus.</p>
{fig.to_html(full_html=False, include_plotlyjs='cdn') if fig is not None else '<p><strong>Plotly is not available;</strong> graphical panels were skipped, but tables below contain the full trajectory.</p>'}
<h2>Selected Candidates</h2>
{_html_table(['Binner', 'Selected iteration', 'Score', 'Pigeon loss', 'PAM', 'Bin recovery', 'Dominant loss component', 'Stop reason'], selected_rows)}
<h2>All Candidate Iterations</h2>
{_html_table(['Binner', 'Iteration', 'Seed', 'Score', 'Delta score', 'Pigeon loss', 'PAM', 'Bin recovery', 'Dominant loss component', 'Gain state', 'Stop reason'], candidate_rows)}
<h2>Output Files</h2>
<ul>
  <li><code>selection/selection_trajectory.tsv</code>: machine-readable trajectory behind this plot.</li>
  <li><code>selection/selected_iterations.tsv</code>: selected candidate per binner.</li>
  <li><code>iterative_summary.json</code>: summary metadata and selected output locations.</li>
</ul>
</body>
</html>
"""
    report_path.write_text(html_doc)
    return report_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run bounded iterative MAG binning and Pigeon-based selection."
    )
    parser.add_argument("--gfa", required=True, type=Path)
    parser.add_argument("--assembly", required=True, type=Path)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--binners", default=",".join(DEFAULT_BINNER_ORDER),
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
                        choices=["composite", "unitig_recovery", "pam", "pigeon_loss", "loss", "bin_explained", "both"])
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

    if args.max_iterations < 1:
        raise SystemExit("--max-iterations must be at least 1")
    if args.threads < 1:
        raise SystemExit("--threads must be at least 1")
    if args.seed_step < 1:
        raise SystemExit("--seed-step must be at least 1")
    for label, path in [("--gfa", args.gfa), ("--assembly", args.assembly), ("--bam", args.bam)]:
        if not path.exists():
            raise SystemExit(f"{label} does not exist: {path}")
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
                if (
                    binner in BAM_REFERENCE_MATCH_BINNER
                    and target_fasta.resolve() != args.assembly
                ):
                    mark_binner_stopped(records, active, binner, "residual_bam_mismatch")
                    log_note(
                        log_handle,
                        candidate_dir,
                        f"stopping {binner} before iteration {iteration}: residual FASTA "
                        "does not match BAM reference set",
                    )
                    continue
                try:
                    bins_dir = runners[binner](args, iteration, target_fasta, candidate_dir, log_handle)
                except subprocess.CalledProcessError as exc:
                    mark_binner_stopped(records, active, binner, "binner_failed")
                    log_note(
                        log_handle,
                        candidate_dir,
                        f"stopping {binner} after failed iteration {iteration}: {exc}",
                    )
                    continue
                if not fasta_files(bins_dir):
                    mark_binner_stopped(records, active, binner, "no_bins")
                    log_note(
                        log_handle,
                        candidate_dir,
                        f"stopping {binner} after iteration {iteration}: no bin FASTA files found",
                    )
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

        consensus_bins = None
        if args.consensus_tool == "dastool":
            consensus_bins = run_dastool(args, selection_out, log_handle)
            materialize_final_bins([consensus_bins], final_bins)
        else:
            run_no_consensus(selection_out, final_bins)

        binner_convergence = summarize_binner_convergence(records, binners)

    summary = {
        "schema_version": "pigeon.iterative_binning.v1",
        "sample_id": args.sample_id,
        "binners": binners,
        "max_iterations": args.max_iterations,
        "n_candidates": len(records),
        "consensus_tool": args.consensus_tool,
        "score_mode": args.selection_score_mode,
        "pigeon_loss_model": PIGEON_LOSS_MODEL,
        "pigeon_loss_components": PIGEON_LOSS_COMPONENTS,
        "residual_rebinning": args.residual_rebinning,
        "residual_ksize": args.residual_ksize,
        "residual_min_unexplained_fraction": args.residual_min_unexplained_fraction,
        "residual_min_contig_len": args.residual_min_contig_len,
        "final_bins": str(final_bins),
        "dastool_bins": str(consensus_bins) if consensus_bins else "",
        "binner_convergence": binner_convergence,
        "candidate_manifest": str(args.outdir / "candidate_manifest.tsv"),
        "selection_summary": str(args.outdir / "selection" / "selection_summary.json"),
    }
    with open(args.outdir / "iterative_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    write_iterative_report(args, records, summary, binner_convergence)


if __name__ == "__main__":
    main()
