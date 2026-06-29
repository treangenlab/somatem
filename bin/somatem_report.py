#!/usr/bin/env python3
"""Create polished, self-contained HTML summary reports for somatem workflows."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import html
import json
import math
import mimetypes
import os
import re
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable
from urllib.parse import quote


WORKFLOW_COPY = {
    "pre_processing": {
        "title": "Pre-processing and quality control",
        "summary": (
            "This report summarizes sample intake, read quality assessment, optional host "
            "depletion, read filtering, and the cleaned read files passed to downstream "
            "somatem analyses."
        ),
        "methods": [
            "Raw long-read quality was assessed with NanoPlot before filtering.",
            "When enabled, host-derived reads were removed with Hostile using the configured host index.",
            "Reads were filtered with Chopper using the configured minimum quality and length thresholds.",
            "Final NanoPlot quality summaries were generated from the cleaned read set.",
        ],
        "accent": "#007c89",
    },
    "16S": {
        "title": "16S taxonomic profiling",
        "summary": (
            "This report summarizes cleaned 16S long-read samples, EMU abundance profiling, "
            "and interactive taxonomy visualization outputs generated for collaborator-facing review."
        ),
        "methods": [
            "Cleaned reads from the pre-processing workflow were profiled with EMU.",
            "Per-sample abundance tables were converted for TaxBurst visualization.",
            "TaxBurst HTML reports provide an interactive view of taxonomic composition.",
        ],
        "accent": "#7c3aed",
    },
    "taxonomic_profiling": {
        "title": "Metagenomic taxonomic profiling",
        "summary": (
            "This report summarizes marker-based taxonomic profiling of cleaned metagenomic reads, "
            "interactive taxonomy visualization, and optional MAGNET presence/absence refinement."
        ),
        "methods": [
            "Cleaned reads were profiled with LEMUR against the configured marker database and taxonomy.",
            "Taxonomic abundance outputs were converted for interactive TaxBurst visualization.",
            "MAGNET was used to refine low-abundance or low-coverage calls when metagenomic profiling was selected.",
        ],
        "accent": "#0f766e",
    },
    "assembly_mags": {
        "title": "Assembly and MAG recovery",
        "summary": (
            "This report summarizes de novo assembly, assembly graph visualization, read mapping, "
            "coverage, genome binning, bin quality control, taxonomic appraisal, Pigeon analyses, "
            "and annotation outputs."
        ),
        "methods": [
            "Cleaned reads were profiled with SingleM and assembled with Flye in metagenome mode.",
            "Assembly graphs were visualized and reads were mapped back to contigs with minimap2 and samtools.",
            "Genome bins were recovered with the configured SemiBin/Pigeon strategy and evaluated with CheckM2.",
            "Recovered bins were appraised taxonomically with SingleM and annotated with Bakta when they passed the completeness threshold.",
        ],
        "accent": "#b45309",
    },
    "genome_dynamics": {
        "title": "Genome dynamics",
        "summary": (
            "This report summarizes the longitudinal genome-dynamics workflow, including joint read "
            "assembly, graph construction, graph visualization, and supporting Rhea output tables."
        ),
        "methods": [
            "Cleaned reads across samples were collected for joint graph-aware analysis.",
            "Rhea generated assembly graph, coverage, and structural-variant supporting outputs.",
            "Bandage rendered assembly graph images for visual inspection and discussion.",
        ],
        "accent": "#be123c",
    },
    "isolate_analysis": {
        "title": "Isolate analysis",
        "summary": (
            "This report summarizes long-read-first bacterial isolate assembly, optional hybrid polishing, "
            "read classification, assembly quality assessment, genome annotation, and optional B. cereus "
            "typing outputs."
        ),
        "methods": [
            "Long reads were classified with Kraken2 against the configured standard-8 database.",
            "Long-read assemblies were generated with Autocycler and Flye, then polished with short reads when hybrid assembly was enabled.",
            "Final isolate assemblies were evaluated with CheckM2 to estimate completeness and contamination.",
            "Assemblies were annotated with Bakta, and optional BTyper3 typing was run from the Bakta genome FASTA output.",
        ],
        "accent": "#2563eb",
    },
}


TYPE_LABELS = [
    ("fastq", "Cleaned FASTQ files", lambda p: p.name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))),
    ("taxonomy", "Taxonomy tables", lambda p: any(x in p.name.lower() for x in ["abundance", "taxonomic_profile", "relative_abundance", "cluster_representative"])),
    ("assembly", "Assembly files", lambda p: p.name.endswith((".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz", ".gfa", ".gfa.gz"))),
    ("mapping", "Mapping and coverage", lambda p: p.name.endswith((".bam", ".bai", ".coverage.txt")) or "coverage" in p.name.lower()),
    ("binning", "Binning outputs", lambda p: any(x in p.name.lower() for x in ["bin", "semibin", "dastool", "iterative", "checkm2", "completeness"])),
    ("annotation", "Annotation outputs", lambda p: p.name.endswith((".gff", ".gff3", ".gbff", ".embl", ".faa", ".ffn")) or "bakta" in p.name.lower()),
    ("visual", "Visual reports and images", lambda p: p.name.endswith((".html", ".png", ".svg"))),
    ("versions", "Software version files", lambda p: p.name == "versions.yml" or p.name.endswith("_versions.yml")),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow", required=True)
    parser.add_argument("--workflow-label", required=True)
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--analysis-type", default="unknown")
    parser.add_argument("--data-type", default="unknown")
    parser.add_argument("--sample-environment", default="unknown")
    parser.add_argument("--sequencing-technology", default="unknown")
    parser.add_argument("--pipeline-version", default="unknown")
    parser.add_argument("--run-name", default="unknown")
    return parser.parse_args()


def esc(value: object) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def slugify(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")


def read_text(path: Path, limit: int = 200_000) -> str:
    try:
        if path.suffix == ".gz":
            with gzip.open(path, "rt", errors="replace") as handle:
                return handle.read(limit)
        with path.open("rt", errors="replace") as handle:
            return handle.read(limit)
    except Exception:
        return ""


def human_bytes(size: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    return f"{size} B"


def collect_files(input_dir: Path) -> list[Path]:
    if not input_dir.exists():
        return []
    ignored_parts = {".git", "__pycache__", ".nextflow", "work"}
    return sorted(
        path
        for path in input_dir.rglob("*")
        if path.is_file() and not ignored_parts.intersection(path.parts)
    )


def file_size(path: Path) -> int:
    try:
        return path.stat().st_size
    except OSError:
        return 0


def read_table(path: Path, max_rows: int = 12) -> tuple[list[str], list[list[str]]]:
    text = read_text(path, 200_000)
    if not text.strip():
        return [], []
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return [], []
    sample = "\n".join(lines[:20])
    if "\t" in lines[0]:
        dialect = csv.excel_tab
    else:
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;")
        except Exception:
            dialect = csv.excel
    rows = list(csv.reader(lines[: max_rows + 1], dialect))
    if not rows:
        return [], []
    header = [cell.strip() for cell in rows[0]]
    body = [[cell.strip() for cell in row] for row in rows[1:]]
    if len(header) == 1 and len(body) <= 1:
        return [], []
    return header, body


def read_table_rows(path: Path, max_rows: int = 10000) -> tuple[list[str], list[dict[str, str]]]:
    header, body = read_table(path, max_rows=max_rows)
    rows = []
    for values in body:
        row = {}
        for idx, name in enumerate(header):
            row[name] = values[idx] if idx < len(values) else ""
        rows.append(row)
    return header, rows


def read_samplesheet(path: Path) -> tuple[list[str], list[list[str]]]:
    header, rows = read_table(path, max_rows=5000)
    return header, rows


def classify_file(path: Path) -> str:
    for key, _label, predicate in TYPE_LABELS:
        if predicate(path):
            return key
    return "other"


def parse_versions(files: Iterable[Path]) -> dict[str, str]:
    versions: dict[str, str] = {}
    for path in files:
        if not (path.name == "versions.yml" or path.name.endswith("_versions.yml")):
            continue
        text = read_text(path)
        current_process = None
        for raw_line in text.splitlines():
            line = raw_line.rstrip()
            if not line.strip():
                continue
            process_match = re.match(r'^"?([^":]+)"?:\s*$', line.strip())
            if process_match and not raw_line.startswith(" "):
                current_process = process_match.group(1).split(":")[-1]
                continue
            tool_match = re.match(r"\s+([^:]+):\s*(.+?)\s*$", line)
            if tool_match:
                tool = tool_match.group(1).strip().strip('"')
                version = tool_match.group(2).strip().strip('"').strip("'")
                label = f"{current_process}: {tool}" if current_process else tool
                versions[label] = version
    return dict(sorted(versions.items()))


def parse_json_summaries(files: Iterable[Path]) -> list[tuple[str, dict[str, object]]]:
    summaries = []
    for path in files:
        if path.suffix.lower() != ".json":
            continue
        text = read_text(path, 300_000)
        try:
            data = json.loads(text)
        except Exception:
            continue
        if isinstance(data, dict):
            flat = {}
            for key, value in data.items():
                if isinstance(value, (str, int, float, bool)) or value is None:
                    flat[key] = value
                elif isinstance(value, list):
                    flat[key] = f"{len(value)} items"
                elif isinstance(value, dict):
                    flat[key] = f"{len(value)} fields"
            summaries.append((path.name, flat))
    return summaries[:8]


def numeric_column_stats(header: list[str], rows: list[list[str]]) -> list[tuple[str, float, float, float]]:
    stats = []
    for col_idx, name in enumerate(header):
        values = []
        for row in rows:
            if col_idx >= len(row):
                continue
            try:
                values.append(float(row[col_idx].replace(",", "")))
            except Exception:
                pass
        if len(values) >= 2:
            stats.append((name, min(values), statistics.median(values), max(values)))
    return stats[:8]


def to_float(value: object) -> float | None:
    try:
        text = str(value).strip().replace(",", "")
        if text in {"", "NA", "NaN", "nan", "None"}:
            return None
        return float(text)
    except Exception:
        return None


def pct(value: float) -> str:
    return f"{value * 100:.1f}%"


def table_html(header: list[str], rows: list[list[str]], max_cols: int = 8) -> str:
    if not header:
        return ""
    clipped_header = header[:max_cols]
    extra_cols = max(0, len(header) - max_cols)
    thead = "".join(f"<th>{esc(col)}</th>" for col in clipped_header)
    if extra_cols:
        thead += "<th>...</th>"
    body = []
    for row in rows:
        clipped = row[:max_cols]
        cells = "".join(f"<td>{esc(cell)}</td>" for cell in clipped)
        if extra_cols:
            cells += f"<td>{extra_cols} more</td>"
        body.append(f"<tr>{cells}</tr>")
    return f"<div class=\"table-wrap\"><table><thead><tr>{thead}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"


def chart_palette() -> list[str]:
    return ["#007c89", "#7c3aed", "#b45309", "#0f766e", "#be123c", "#2563eb", "#ca8a04", "#4b5563"]


def sample_name_from_path(path: Path) -> str:
    name = path.name
    suffixes = [
        ".fastq_rel-abundance.tsv",
        ".fq_rel-abundance.tsv",
        "_rel-abundance.tsv",
        "_relative_abundance.tsv",
        "-relative_abundance.tsv",
        "_taxonomic_profile.tsv",
        ".taxonomic_profile.tsv",
        ".tsv",
        ".csv",
        ".txt",
    ]
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def first_present(header: list[str], candidates: list[str]) -> str | None:
    lookup = {name.lower(): name for name in header}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def taxonomy_table_paths(files: list[Path]) -> list[Path]:
    paths = []
    for path in files:
        name = path.name.lower()
        if path.suffix.lower() not in {".tsv", ".csv", ".txt"}:
            continue
        if path.name == "versions.yml" or path.name.endswith("_versions.yml"):
            continue
        if any(token in name for token in ["abundance", "relative_abundance", "taxonomic_profile", "otu_table"]):
            paths.append(path)
    return paths


def taxonomy_profiles(files: list[Path], rank: str = "phylum") -> dict[str, dict[str, float]]:
    profiles: dict[str, dict[str, float]] = {}
    for path in taxonomy_table_paths(files):
        header, rows = read_table_rows(path)
        if not header or not rows:
            continue
        abundance_col = first_present(header, ["abundance", "F", "relative_abundance", "fraction", "percent", "percentage"])
        rank_col = first_present(header, [rank, rank.capitalize(), "species", "genus", "taxon", "name"])
        if not abundance_col or not rank_col:
            continue
        sample = sample_name_from_path(path)
        values: dict[str, float] = defaultdict(float)
        for row in rows:
            abundance = to_float(row.get(abundance_col))
            taxon = str(row.get(rank_col, "")).strip() or "Unclassified"
            if abundance is None or abundance <= 0:
                continue
            values[taxon] += abundance
        total = sum(values.values())
        if total <= 0:
            continue
        if total > 1.5:
            values = defaultdict(float, {taxon: value / total for taxon, value in values.items()})
        profiles[sample] = dict(values)
    return profiles


def top_taxa(profiles: dict[str, dict[str, float]], limit: int = 7) -> list[str]:
    totals: Counter[str] = Counter()
    for composition in profiles.values():
        totals.update(composition)
    return [taxon for taxon, _value in totals.most_common(limit)]


def stacked_composition_svg(profiles: dict[str, dict[str, float]], title: str) -> str:
    if not profiles:
        return ""
    samples = sorted(profiles)
    taxa = top_taxa(profiles, limit=7)
    palette = chart_palette()
    width = 900
    row_h = 42
    left = 170
    top = 36
    bar_w = 560
    height = top + len(samples) * row_h + 86
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-label="{esc(title)}">',
        f'<text x="{left}" y="20" class="chart-title">{esc(title)}</text>',
    ]
    for row_idx, sample in enumerate(samples):
        y = top + row_idx * row_h
        parts.append(f'<text x="0" y="{y + 20}" class="axis-label">{esc(sample)}</text>')
        x = left
        shown = 0.0
        composition = profiles[sample]
        for idx, taxon in enumerate(taxa):
            value = composition.get(taxon, 0.0)
            if value <= 0:
                continue
            w = max(0.0, min(bar_w * value, bar_w - (x - left)))
            if w > 0:
                parts.append(
                    f'<rect x="{x:.2f}" y="{y}" width="{w:.2f}" height="24" fill="{palette[idx % len(palette)]}">'
                    f'<title>{esc(sample)} - {esc(taxon)}: {pct(value)}</title></rect>'
                )
                x += w
                shown += value
        other = max(0.0, 1.0 - shown)
        if other > 0.005 and x < left + bar_w:
            w = min(bar_w * other, left + bar_w - x)
            parts.append(
                f'<rect x="{x:.2f}" y="{y}" width="{w:.2f}" height="24" fill="#cbd5e1">'
                f'<title>{esc(sample)} - Other: {pct(other)}</title></rect>'
            )
        parts.append(f'<text x="{left + bar_w + 12}" y="{y + 18}" class="axis-label">100%</text>')
    legend_y = top + len(samples) * row_h + 24
    legend_items = taxa + ["Other"]
    for idx, taxon in enumerate(legend_items):
        x = left + (idx % 4) * 170
        y = legend_y + (idx // 4) * 22
        color = palette[idx % len(palette)] if taxon != "Other" else "#cbd5e1"
        parts.append(f'<rect x="{x}" y="{y - 11}" width="12" height="12" fill="{color}"></rect>')
        parts.append(f'<text x="{x + 18}" y="{y}" class="legend-label">{esc(taxon)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def horizontal_bar_svg(values: list[tuple[str, float]], title: str, value_suffix: str = "") -> str:
    if not values:
        return ""
    values = values[:10]
    max_value = max(value for _name, value in values) or 1.0
    width = 900
    left = 230
    bar_w = 520
    row_h = 34
    top = 36
    height = top + len(values) * row_h + 28
    palette = chart_palette()
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-label="{esc(title)}">',
        f'<text x="{left}" y="20" class="chart-title">{esc(title)}</text>',
    ]
    for idx, (name, value) in enumerate(values):
        y = top + idx * row_h
        w = bar_w * (value / max_value if max_value else 0)
        label = f"{value:.3g}{value_suffix}"
        parts.append(f'<text x="0" y="{y + 18}" class="axis-label">{esc(name)}</text>')
        parts.append(f'<rect x="{left}" y="{y}" width="{w:.2f}" height="20" fill="{palette[idx % len(palette)]}"></rect>')
        parts.append(f'<text x="{left + w + 10}" y="{y + 16}" class="axis-label">{esc(label)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def histogram_svg(values: list[float], title: str, value_suffix: str = "", bins: int = 10) -> str:
    values = [value for value in values if value is not None]
    if not values:
        return ""
    low = min(values)
    high = max(values)
    if low == high:
        low = 0.0
        high = max(high, 1.0)
    bin_count = max(3, min(bins, 12))
    step = (high - low) / bin_count
    counts = [0] * bin_count
    for value in values:
        idx = min(int((value - low) / step), bin_count - 1)
        counts[idx] += 1
    width = 900
    left = 70
    top = 46
    chart_w = 700
    chart_h = 220
    max_count = max(counts) or 1
    bar_gap = 8
    bar_w = (chart_w - bar_gap * (bin_count - 1)) / bin_count
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} 340" role="img" aria-label="{esc(title)}">',
        f'<text x="{left}" y="22" class="chart-title">{esc(title)}</text>',
        f'<line x1="{left}" y1="{top + chart_h}" x2="{left + chart_w}" y2="{top + chart_h}" class="chart-axis"></line>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_h}" class="chart-axis"></line>',
    ]
    for idx, count in enumerate(counts):
        x = left + idx * (bar_w + bar_gap)
        h = chart_h * count / max_count
        y = top + chart_h - h
        start = low + idx * step
        end = start + step
        label = f"{start:.3g}-{end:.3g}{value_suffix}"
        parts.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{h:.2f}" fill="{chart_palette()[idx % len(chart_palette())]}">'
            f'<title>{esc(label)}: {count}</title></rect>'
        )
        parts.append(f'<text x="{x + bar_w / 2:.2f}" y="{top + chart_h + 18}" class="axis-label" text-anchor="middle">{esc(label)}</text>')
        parts.append(f'<text x="{x + bar_w / 2:.2f}" y="{max(y - 6, top + 12):.2f}" class="axis-label" text-anchor="middle">{count}</text>')
    parts.append("</svg>")
    return "".join(parts)


def scatter_svg(points: list[tuple[str, float, float]], title: str, x_label: str, y_label: str) -> str:
    if not points:
        return ""
    width = 900
    left = 82
    top = 44
    chart_w = 660
    chart_h = 260
    x_max = max(100.0, max(point[1] for point in points) * 1.08)
    y_max = max(10.0, max(point[2] for point in points) * 1.18)
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} 370" role="img" aria-label="{esc(title)}">',
        f'<text x="{left}" y="22" class="chart-title">{esc(title)}</text>',
        f'<line x1="{left}" y1="{top + chart_h}" x2="{left + chart_w}" y2="{top + chart_h}" class="chart-axis"></line>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_h}" class="chart-axis"></line>',
        f'<text x="{left + chart_w / 2}" y="{top + chart_h + 48}" class="axis-label" text-anchor="middle">{esc(x_label)}</text>',
        f'<text x="18" y="{top + chart_h / 2}" class="axis-label" transform="rotate(-90 18 {top + chart_h / 2})" text-anchor="middle">{esc(y_label)}</text>',
    ]
    palette = chart_palette()
    for idx, (label, x_value, y_value) in enumerate(points):
        x = left + chart_w * x_value / x_max
        y = top + chart_h - (chart_h * y_value / y_max)
        parts.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="6" fill="{palette[idx % len(palette)]}" opacity="0.82">'
            f'<title>{esc(label)}: completeness {x_value:.1f}%, contamination {y_value:.1f}%</title></circle>'
        )
    parts.append("</svg>")
    return "".join(parts)


def diversity_metrics(profiles: dict[str, dict[str, float]]) -> list[list[str]]:
    rows = []
    for sample, composition in sorted(profiles.items()):
        values = [value for value in composition.values() if value > 0]
        total = sum(values)
        if total <= 0:
            continue
        proportions = [value / total for value in values]
        shannon = -sum(value * math.log2(value) for value in proportions if value > 0)
        top_taxon, top_value = max(composition.items(), key=lambda item: item[1])
        rows.append([sample, len(values), f"{shannon:.2f}", top_taxon, pct(top_value / total)])
    return rows


def taxonomy_figure_section(files: list[Path], workflow: str) -> str:
    if workflow not in {"16S", "taxonomic_profiling", "assembly_mags"}:
        return ""
    phylum_profiles = taxonomy_profiles(files, rank="phylum")
    genus_profiles = taxonomy_profiles(files, rank="genus")
    species_profiles = taxonomy_profiles(files, rank="species")
    if not phylum_profiles and not genus_profiles and not species_profiles:
        return ""
    body = (
        "<p>These quick-look summaries are derived from the staged taxonomic abundance tables. "
        "They are intended to highlight broad composition patterns and dominant taxa before opening "
        "the full interactive TaxBurst reports.</p>"
    )
    if phylum_profiles:
        body += stacked_composition_svg(phylum_profiles, "Relative abundance by phylum")
    if workflow == "16S" and genus_profiles:
        body += stacked_composition_svg(genus_profiles, "Relative abundance by genus")
    if workflow == "taxonomic_profiling" and species_profiles:
        body += stacked_composition_svg(species_profiles, "Relative abundance by species")
    if species_profiles:
        totals: Counter[str] = Counter()
        for composition in species_profiles.values():
            total = sum(composition.values()) or 1.0
            for taxon, value in composition.items():
                totals[taxon] += value / total
        averaged = [(taxon, 100 * value / max(len(species_profiles), 1)) for taxon, value in totals.most_common(10)]
        body += horizontal_bar_svg(averaged, "Most abundant taxa across staged samples", value_suffix="%")
        rows = diversity_metrics(species_profiles)
        if rows:
            body += "<h3>Sample-level Diversity Summary</h3>"
            body += (
                "<p>Observed taxa counts and Shannon diversity are calculated from non-zero relative "
                "abundance entries in each staged taxonomic profile. These values are descriptive "
                "screening metrics, not formal ecological inference.</p>"
            )
            body += table_html(["Sample", "Observed taxa", "Shannon diversity", "Dominant taxon", "Dominant taxon abundance"], rows, max_cols=5)
    return body


def checkm2_rows(files: list[Path]) -> list[dict[str, str]]:
    rows_out = []
    for path in files:
        name = path.name.lower()
        if "checkm2" not in name and "completeness" not in name:
            continue
        header, rows = read_table_rows(path)
        if not header:
            continue
        completeness_col = first_present(header, ["Completeness", "completeness"])
        contamination_col = first_present(header, ["Contamination", "contamination"])
        if not completeness_col:
            continue
        id_col = first_present(header, ["Name", "Bin Id", "Bin", "Genome", "sample", "id"]) or header[0]
        for row in rows:
            completeness = to_float(row.get(completeness_col))
            contamination = to_float(row.get(contamination_col)) if contamination_col else None
            if completeness is None:
                continue
            rows_out.append(
                {
                    "id": row.get(id_col, "unknown"),
                    "completeness": str(completeness),
                    "contamination": "" if contamination is None else str(contamination),
                }
            )
    return rows_out


def mag_quality_section(files: list[Path], workflow: str) -> str:
    if workflow not in {"assembly_mags", "isolate_analysis"}:
        return ""
    is_isolate = workflow == "isolate_analysis"
    entity_label = "assembly" if is_isolate else "bin"
    entity_label_title = "Assembly" if is_isolate else "Bin"
    rows = checkm2_rows(files)
    body = ""
    contig_lengths = []
    contig_coverages = []
    for path in files:
        name = path.name.lower()
        if path.suffix.lower() not in {".txt", ".tsv", ".csv"}:
            continue
        if not any(token in name for token in ["assembly", "flye", "contig", "coverage"]):
            continue
        header, table_rows = read_table_rows(path)
        length_col = first_present(header, ["length", "len", "contig_length", "seq_len"])
        coverage_col = first_present(header, ["cov", "coverage", "mean_coverage", "coverage_mean", "meandepth", "depth", "avg_depth"])
        id_col = first_present(header, ["seq_name", "contig", "contig_id", "name", "id"]) or (header[0] if header else None)
        if length_col and id_col:
            for row in table_rows:
                value = to_float(row.get(length_col))
                if value is not None and value > 0:
                    contig_lengths.append((row.get(id_col, "unknown"), value))
        if coverage_col and id_col:
            for row in table_rows:
                value = to_float(row.get(coverage_col))
                if value is not None and value > 0:
                    contig_coverages.append((row.get(id_col, "unknown"), value))
    if contig_lengths:
        body += (
            "<h3>Assembly Contiguity Overview</h3>"
            "<p>Contig length summaries are calculated from staged assembly tables when available. "
            "These plots provide a quick check of assembly fragmentation and identify the longest assembled sequences.</p>"
        )
        body += horizontal_bar_svg(sorted(contig_lengths, key=lambda item: item[1], reverse=True), "Longest assembled contigs", value_suffix=" bp")
        body += histogram_svg([value for _name, value in contig_lengths], "Contig length distribution", value_suffix=" bp")
    if rows:
        high_quality = 0
        medium_quality = 0
        low_quality = 0
        chart_values = []
        contamination_values = []
        scatter_points = []
        table_rows = []
        for row in rows:
            completeness = to_float(row["completeness"]) or 0.0
            contamination = to_float(row["contamination"])
            if completeness >= 90 and (contamination is None or contamination <= 5):
                high_quality += 1
                quality_label = "High-quality"
            elif completeness >= 50 and (contamination is None or contamination <= 10):
                medium_quality += 1
                quality_label = "Medium-quality"
            else:
                low_quality += 1
                quality_label = "Lower-quality"
            chart_values.append((row["id"], completeness))
            if contamination is not None:
                contamination_values.append((row["id"], contamination))
                scatter_points.append((row["id"], completeness, contamination))
            table_rows.append([
                row["id"],
                f"{completeness:.1f}%",
                "not reported" if contamination is None else f"{contamination:.1f}%",
                quality_label,
            ])
        body += (
            "<p>Quality summaries use CheckM2 completeness and contamination estimates. "
            f"This run staged {len(rows)} evaluated {entity_label}{'' if len(rows) == 1 else 's'}, including "
            f"{high_quality} high-quality {entity_label}{'' if high_quality == 1 else 's'} "
            f"(at least 90% complete and at most 5% contamination), {medium_quality} "
            f"medium-quality {entity_label}{'' if medium_quality == 1 else 's'} "
            "(at least 50% complete and at most 10% contamination), and "
            f"{low_quality} lower-quality {entity_label}{'' if low_quality == 1 else 's'}.</p>"
        )
        body += horizontal_bar_svg(
            [
                ("High-quality", high_quality),
                ("Medium-quality", medium_quality),
                ("Lower-quality", low_quality),
            ],
            f"{entity_label_title} quality class counts",
        )
        body += horizontal_bar_svg(sorted(chart_values, key=lambda item: item[1], reverse=True), f"{entity_label_title} completeness estimates", value_suffix="%")
        body += histogram_svg([value for _name, value in chart_values], "Completeness distribution", value_suffix="%")
        if contamination_values:
            body += horizontal_bar_svg(
                sorted(contamination_values, key=lambda item: item[1], reverse=True),
                f"{entity_label_title} contamination estimates",
                value_suffix="%",
            )
            body += scatter_svg(scatter_points, "Completeness versus contamination", "Completeness (%)", "Contamination (%)")
        body += table_html([entity_label_title, "Completeness", "Contamination", "Quality class"], table_rows[:20], max_cols=4)
    coverage_values = list(contig_coverages)
    for path in files:
        if "coverage" not in path.name.lower() or path.suffix.lower() not in {".txt", ".tsv", ".csv"}:
            continue
        header, table_rows = read_table_rows(path)
        coverage_col = first_present(header, ["coverage", "mean_coverage", "coverage_mean", "meandepth", "depth", "avg_depth"])
        id_col = first_present(header, ["rname", "contig", "contig_id", "name", "id"]) or (header[0] if header else None)
        if not coverage_col or not id_col:
            continue
        for row in table_rows:
            value = to_float(row.get(coverage_col))
            if value is not None:
                coverage_values.append((row.get(id_col, "unknown"), value))
    if coverage_values:
        body += (
            "<h3>Coverage Overview</h3>"
            "<p>The coverage plot highlights the highest-coverage contigs or bins in the staged mapping summary. "
            "Coverage is useful for spotting dominant genomes and uneven assembly support.</p>"
        )
        body += horizontal_bar_svg(sorted(coverage_values, key=lambda item: item[1], reverse=True), "Highest coverage contigs or bins", value_suffix="x")
        body += histogram_svg([value for _name, value in coverage_values], "Coverage distribution", value_suffix="x")
    binning_counts = Counter()
    for path in files:
        name = path.name.lower()
        if not any(token in name for token in ["bin", "pigeon", "manifest", "trajectory", "selected"]):
            continue
        if path.suffix.lower() not in {".txt", ".tsv", ".csv"}:
            continue
        header, table_rows = read_table_rows(path)
        if not header or not table_rows:
            continue
        status_col = first_present(header, ["status", "selected", "decision", "source", "binner", "method", "iteration"])
        if status_col:
            for row in table_rows:
                value = str(row.get(status_col, "")).strip() or "reported"
                binning_counts[value] += 1
        else:
            binning_counts[path.name] += len(table_rows)
    if binning_counts:
        body += (
            "<h3>Binning Output Overview</h3>"
            "<p>This summary counts entries from staged binning, refinement, or Pigeon metadata tables. "
            "It provides a compact view of which binning decisions or output classes are represented in the run.</p>"
        )
        body += horizontal_bar_svg(binning_counts.most_common(10), "Binning and refinement record counts")
    annotation_counts = Counter()
    for path in files:
        name = path.name.lower()
        if not any(token in name for token in ["bakta", "annotation"]):
            continue
        annotation_counts[path.suffix.lower() or "reported"] += 1
    if annotation_counts:
        body += (
            "<h3>Annotation Output Overview</h3>"
            "<p>Annotation outputs are counted by file type to show which Bakta result artifacts were staged "
            "for downstream review.</p>"
        )
        body += horizontal_bar_svg(annotation_counts.most_common(), "Annotation output files by type")
    return body


def table_description(path: Path, header: list[str], workflow: str) -> str:
    name = path.name.lower()
    columns = {col.lower() for col in header}
    if any(token in name for token in ["abundance", "relative_abundance", "taxonomic_profile"]):
        rank_terms = [rank for rank in ["species", "genus", "family", "phylum", "superkingdom"] if rank in columns]
        ranks = ", ".join(rank_terms) if rank_terms else "taxonomic ranks"
        return (
            f"<p>This table previews taxonomic abundance calls from <strong>{esc(path.name)}</strong>. "
            f"It includes relative abundance estimates and taxonomy fields ({esc(ranks)}) that support "
            "sample composition review and downstream visualization.</p>"
        )
    if "checkm2" in name or "completeness" in name:
        return (
            f"<p>This table previews MAG quality estimates from <strong>{esc(path.name)}</strong>. "
            "Completeness and contamination values help prioritize bins for interpretation, annotation, "
            "and follow-up quality control.</p>"
        )
    if "coverage" in name:
        return (
            f"<p>This table previews read-mapping coverage from <strong>{esc(path.name)}</strong>. "
            "Coverage summaries help assess assembly support and identify contigs or bins with notably "
            "high or low read depth.</p>"
        )
    if any(token in name for token in ["pigeon", "bin", "manifest", "trajectory"]):
        return (
            f"<p>This table previews binning and genome-recovery metadata from <strong>{esc(path.name)}</strong>. "
            "These fields document which bins were selected, how they moved through refinement steps, "
            "and which outputs are available for review.</p>"
        )
    return (
        f"<p>This table previews <strong>{esc(path.name)}</strong>, one of the staged workflow outputs. "
        "Only the first rows and columns are shown here; the complete file is available in the published results.</p>"
    )




def inventory_table(files: list[Path]) -> str:
    rows = []
    for path in files[:120]:
        rows.append(
            [
                path.name,
                classify_file(path).replace("_", " "),
                human_bytes(file_size(path)),
                esc("/".join(path.parts[-3:])),
            ]
        )
    return table_html(["File", "Category", "Size", "Staged path"], rows, max_cols=4)


def table_previews(files: list[Path]) -> list[tuple[Path, str, list[tuple[str, float, float, float]]]]:
    candidates = []
    skip_suffixes = {".html", ".png", ".svg", ".bam", ".bai", ".gz", ".fa", ".fasta", ".fna", ".gfa", ".gbff", ".embl", ".faa", ".ffn"}
    for path in files:
        if path.suffix.lower() in skip_suffixes:
            continue
        if file_size(path) > 3_000_000:
            continue
        header, rows = read_table(path, max_rows=10)
        if header and rows:
            stats = numeric_column_stats(header, rows)
            candidates.append((path, table_html(header, rows), stats))
    return candidates[:10]


def image_previews(files: list[Path]) -> list[tuple[str, str]]:
    previews = []
    for path in files:
        if path.suffix.lower() not in {".png", ".svg"}:
            continue
        size_limit = 5_000_000 if path.suffix.lower() == ".png" else 800_000
        if file_size(path) > size_limit:
            continue
        mime = mimetypes.guess_type(path.name)[0] or "image/png"
        try:
            payload = base64.b64encode(path.read_bytes()).decode("ascii")
        except Exception:
            continue
        previews.append((path.name, f"data:{mime};base64,{payload}"))
        if len(previews) >= 4:
            break
    return previews


def taxburst_report_links(files: list[Path]) -> list[tuple[str, str, str]]:
    links = []
    for path in files:
        if path.suffix.lower() != ".html" or "taxburst" not in path.name.lower():
            continue
        sample_id = re.sub(r"_taxburst\.html$", "", path.name, flags=re.IGNORECASE)
        if sample_id == path.name:
            sample_id = path.stem
        href = f"../taxonomy/{quote(sample_id)}/{quote(path.name)}"
        links.append((sample_id, path.name, href))
    return sorted(links)


def interactive_report_links(files: list[Path]) -> str:
    taxburst_links = taxburst_report_links(files)
    if not taxburst_links:
        return ""
    rows = "".join(
        "<tr>"
        f"<td>{esc(sample_id)}</td>"
        f"<td>{esc(filename)}</td>"
        f"<td><a href=\"{esc(href)}\">Open TaxBurst report</a></td>"
        f"<td>{esc(href)}</td>"
        "</tr>"
        for sample_id, filename, href in taxburst_links
    )
    table = (
        "<div class=\"table-wrap\"><table><thead><tr>"
        "<th>Sample</th><th>Report file</th><th>Link</th><th>Published path</th>"
        f"</tr></thead><tbody>{rows}</tbody></table></div>"
    )
    return (
        "<p>Interactive TaxBurst reports are published alongside this summary under the taxonomy results directory.</p>"
        + table
    )


def metric_cards(sample_count: int, files: list[Path], versions: dict[str, str]) -> str:
    counts = Counter(classify_file(path) for path in files)
    cards = [
        ("Samples", sample_count),
        ("Staged outputs", len(files)),
        ("Result tables", sum(1 for p in files if p.suffix.lower() in {".tsv", ".csv", ".txt"})),
        ("Visual outputs", counts["visual"]),
        ("Version entries", len(versions)),
    ]
    return "<div class=\"metric-grid\">" + "".join(
        f"<div class=\"metric\"><span>{esc(label)}</span><strong>{esc(value)}</strong></div>"
        for label, value in cards
    ) + "</div>"


def category_summary(files: list[Path]) -> str:
    counts = Counter(classify_file(path) for path in files)
    rows = []
    label_by_key = {key: label for key, label, _predicate in TYPE_LABELS}
    for key, count in sorted(counts.items()):
        label = label_by_key.get(key, "Other outputs")
        total_size = sum(file_size(path) for path in files if classify_file(path) == key)
        rows.append([label, count, human_bytes(total_size)])
    return table_html(["Output class", "Files", "Total staged size"], rows, max_cols=3)


def params_table(args: argparse.Namespace) -> str:
    rows = [
        ["Workflow", args.workflow_label],
        ["Analysis type", args.analysis_type],
        ["Data type", args.data_type],
        ["Sequencing technology", args.sequencing_technology],
        ["Sample environment", args.sample_environment],
        ["Pipeline version", args.pipeline_version],
        ["Nextflow run name", args.run_name],
        ["Report generated", datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
    ]
    return table_html(["Field", "Value"], rows, max_cols=2)


def build_toc(sections: list[tuple[str, str]]) -> str:
    links = "".join(f"<a href=\"#{esc(anchor)}\">{esc(title)}</a>" for anchor, title in sections)
    return f"<nav class=\"toc\"><div class=\"toc-title\">Contents</div>{links}</nav>"


def section(title: str, body: str) -> tuple[str, str, str]:
    anchor = slugify(title)
    return anchor, title, f"<section id=\"{esc(anchor)}\"><h2>{esc(title)}</h2>{body}</section>"


def main() -> int:
    args = parse_args()
    files = collect_files(args.input_dir)
    sample_header, sample_rows = read_samplesheet(args.samplesheet)
    sample_count = len(sample_rows)
    versions = parse_versions(files)
    workflow_info = WORKFLOW_COPY.get(
        args.workflow,
        {
            "title": args.workflow_label,
            "summary": "This report summarizes outputs generated by a somatem workflow.",
            "methods": [],
            "accent": "#007c89",
        },
    )

    sections = []
    sections.append(section(
        "Abstract",
        f"<p class=\"lead\">{esc(workflow_info['summary'])}</p>"
        + metric_cards(sample_count, files, versions),
    ))
    sections.append(section(
        "Data Input And Metadata",
        "<p>The table below records the samplesheet rows used to seed this workflow run.</p>"
        + table_html(sample_header, sample_rows[:100], max_cols=10),
    ))
    sections.append(section(
        "Workflow Outputs",
        "<p>This section inventories the files staged into the summary report. The categories distinguish "
        "taxonomic profiles, genome assembly outputs, binning or annotation files, visual reports, and "
        "software version records so collaborators can quickly see which evidence files are available.</p>"
        + category_summary(files)
        + "<h3>Output Inventory</h3>"
        + inventory_table(files),
    ))

    taxonomy_figures = taxonomy_figure_section(files, args.workflow)
    if taxonomy_figures:
        sections.append(section("Microbiome Composition Overview", taxonomy_figures))

    mag_figures = mag_quality_section(files, args.workflow)
    if mag_figures:
        quality_title = "Assembly Quality Overview" if args.workflow == "isolate_analysis" else "MAG Quality And Coverage Overview"
        sections.append(section(quality_title, mag_figures))

    preview_blocks = []
    for path, preview, stats in table_previews(files):
        stat_rows = [[metric, f"{low:.3g}", f"{median:.3g}", f"{high:.3g}"] for metric, low, median, high in stats]
        stat_html = table_html(["Metric", "Min", "Median", "Max"], stat_rows, max_cols=4) if stat_rows else ""
        header, _rows = read_table(path, max_rows=2)
        preview_blocks.append(
            f"<article class=\"result-block\"><h3>{esc(path.name)}</h3>"
            + table_description(path, header, args.workflow)
            + preview
            + stat_html
            + "</article>"
        )
    sections.append(section(
        "Result Tables",
        "<p>Selected tables are previewed below with a short description of how each file contributes to interpretation. "
        "Large binary files and long sequence files are intentionally excluded from preview.</p>"
        + ("".join(preview_blocks) if preview_blocks else "<p>No small tabular outputs were available for preview in this workflow.</p>"),
    ))

    interactive_links = interactive_report_links(files)
    if interactive_links:
        sections.append(section("Interactive Taxonomy Reports", interactive_links))

    image_blocks = [
        f"<figure><img src=\"{src}\" alt=\"{esc(name)}\"><figcaption>{esc(name)}</figcaption></figure>"
        for name, src in image_previews(files)
    ]
    sections.append(section(
        "Visual Outputs",
        "<div class=\"figure-grid\">" + "".join(image_blocks) + "</div>"
        if image_blocks else "<p>No small PNG or SVG outputs were available for embedding in this report.</p>",
    ))

    json_blocks = []
    for name, values in parse_json_summaries(files):
        rows = [[key, value] for key, value in values.items()]
        json_blocks.append(f"<article class=\"result-block\"><h3>{esc(name)}</h3>{table_html(['Field', 'Value'], rows, max_cols=2)}</article>")
    if json_blocks:
        sections.append(section("Structured Summaries", "".join(json_blocks)))

    methods = "".join(f"<li>{esc(item)}</li>" for item in workflow_info.get("methods", []))
    sections.append(section(
        "Methods",
        "<p>The methods summary below describes the analysis steps represented in this report. "
        "It is written for collaborator review and should be interpreted together with the full workflow outputs.</p>"
        + f"<ul class=\"method-list\">{methods}</ul>"
        + "<h3>Run Configuration</h3>"
        + params_table(args),
    ))

    version_rows = [[tool, version] for tool, version in versions.items()]
    sections.append(section(
        "Software Versions",
        table_html(["Process and tool", "Version"], version_rows, max_cols=2)
        if version_rows else "<p>No software version files were staged for this report.</p>",
    ))

    sections.append(section(
        "Final Notes",
        "<p>Use this report as a high-level index into the workflow outputs. Detailed interactive HTML outputs "
        "such as TaxBurst, Pigeon, NanoPlot, or graph visualizations remain available in their corresponding "
        "published result directories.</p>",
    ))

    section_links = [(anchor, title) for anchor, title, _html in sections]
    sections_html = "".join(block for _anchor, _title, block in sections)
    title = workflow_info["title"]
    accent = workflow_info["accent"]

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)} | somatem summary report</title>
<style>
:root {{
  --accent: {accent};
  --ink: #1f2933;
  --muted: #5f6c7b;
  --line: #d8dee7;
  --paper: #ffffff;
  --soft: #f6f8fb;
  --warm: #f59e0b;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  color: var(--ink);
  background: #edf2f7;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  line-height: 1.55;
}}
.shell {{
  display: grid;
  grid-template-columns: 260px minmax(0, 1fr);
  min-height: 100vh;
}}
.toc {{
  position: sticky;
  top: 0;
  height: 100vh;
  overflow: auto;
  padding: 28px 22px;
  background: #17212b;
  color: #f8fafc;
}}
.toc-title {{
  font-size: 0.78rem;
  text-transform: uppercase;
  letter-spacing: 0.09em;
  color: #b8c2cc;
  margin-bottom: 14px;
}}
.toc a {{
  display: block;
  color: #f8fafc;
  text-decoration: none;
  padding: 9px 0;
  border-bottom: 1px solid rgba(255,255,255,0.09);
  font-size: 0.94rem;
}}
.report {{
  max-width: 1180px;
  width: 100%;
  margin: 0 auto;
  background: var(--paper);
  min-height: 100vh;
}}
header {{
  padding: 46px 54px 38px;
  border-top: 8px solid var(--accent);
  background: linear-gradient(180deg, #ffffff 0%, #f8fafc 100%);
  border-bottom: 1px solid var(--line);
}}
.eyebrow {{
  color: var(--accent);
  font-weight: 700;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  font-size: 0.78rem;
}}
h1 {{
  margin: 8px 0 10px;
  font-size: clamp(2rem, 5vw, 3.2rem);
  line-height: 1.05;
  letter-spacing: 0;
}}
.subtitle {{
  max-width: 820px;
  color: var(--muted);
  font-size: 1.08rem;
}}
main {{ padding: 32px 54px 60px; }}
section {{
  padding: 28px 0;
  border-bottom: 1px solid var(--line);
}}
h2 {{
  margin: 0 0 14px;
  font-size: 1.45rem;
  letter-spacing: 0;
}}
h3 {{
  margin: 22px 0 10px;
  font-size: 1rem;
  color: #293845;
}}
.lead {{
  font-size: 1.08rem;
  max-width: 900px;
}}
.metric-grid {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: 12px;
  margin: 24px 0 4px;
}}
.metric {{
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 14px 16px;
  background: var(--soft);
}}
.metric span {{
  display: block;
  color: var(--muted);
  font-size: 0.82rem;
}}
.metric strong {{
  display: block;
  margin-top: 4px;
  font-size: 1.5rem;
  color: var(--ink);
}}
.table-wrap {{
  width: 100%;
  overflow: auto;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fff;
  margin: 12px 0;
}}
table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 0.9rem;
}}
th, td {{
  padding: 9px 11px;
  border-bottom: 1px solid #edf1f5;
  text-align: left;
  vertical-align: top;
  white-space: nowrap;
}}
th {{
  position: sticky;
  top: 0;
  background: #f1f5f9;
  color: #26323f;
  font-weight: 700;
}}
.result-block {{
  border-left: 4px solid var(--accent);
  padding: 2px 0 8px 16px;
  margin: 18px 0 26px;
}}
.chart {{
  width: 100%;
  height: auto;
  display: block;
  margin: 18px 0 26px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fff;
  padding: 12px;
}}
.chart-title {{
  font-size: 16px;
  font-weight: 700;
  fill: #26323f;
}}
.axis-label, .legend-label {{
  font-size: 12px;
  fill: #4b5563;
}}
.chart-axis {{
  stroke: #94a3b8;
  stroke-width: 1;
}}
.figure-grid {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 18px;
}}
figure {{
  margin: 0;
  padding: 12px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fff;
}}
figure img {{
  width: 100%;
  max-height: 360px;
  object-fit: contain;
  display: block;
}}
figcaption {{
  margin-top: 8px;
  color: var(--muted);
  font-size: 0.86rem;
}}
.method-list li {{ margin: 7px 0; }}
@media (max-width: 900px) {{
  .shell {{ display: block; }}
  .toc {{ position: static; height: auto; }}
  header, main {{ padding-left: 24px; padding-right: 24px; }}
  th, td {{ white-space: normal; }}
}}
</style>
</head>
<body>
<div class="shell">
{build_toc(section_links)}
<div class="report">
<header>
  <div class="eyebrow">somatem summary report</div>
  <h1>{esc(title)}</h1>
  <div class="subtitle">{esc(args.workflow_label)} report for run {esc(args.run_name)}. Generated {esc(datetime.now().strftime("%B %d, %Y"))}.</div>
</header>
<main>
{sections_html}
</main>
</div>
</div>
</body>
</html>
"""

    args.output.write_text(document, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
