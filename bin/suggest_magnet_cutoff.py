#!/usr/bin/env python3
"""
Suggest MAGnet --min-abundance cutoff from Lemur abundance table.

This script analyzes a Lemur classifier abundance table and suggests an optimal
--min-abundance cutoff for MAGnet reference selection. The goal is to capture
a target fraction (default 99%) of the total community abundance while minimizing
the number of references to download and process.

Algorithm:
  1. Read Lemur abundance table (TSV format)
  2. Extract relative abundance column (tries common names)
  3. Keep only positive abundances
  4. Sort descending and compute cumulative sum
  5. Find the abundance threshold where cumulative sum reaches target_frac
  6. Output: cutoff value, number of taxa kept, fraction explained

Author: Austin Marshall
Date: 2026-01-07
"""

import argparse
import sys
from pathlib import Path
from typing import Tuple

import pandas as pd


def infer_abundance_col_name(classifier_table: Path) -> str:
    """
    Infer the name of the relative abundance column from common patterns.
    Falls back to 'F' (seen in some Lemur outputs) and then a numeric-only column heuristic.
    """
    df0 = pd.read_csv(classifier_table, sep="\t", nrows=0)
    cols = df0.columns.tolist()
    cols_lower = [c.lower() for c in cols]

    # Priority order of column names to check
    candidates = [
        "relative_abundance",
        "rel_abundance",
        "abundance",
        "fraction",
        "relative abundance",
        "rel abundance",
        # Lemur sometimes uses a single-letter column for fractions
        "f",
    ]

    for candidate in candidates:
        cand = candidate.lower()
        if cand in cols_lower:
            return cols[cols_lower.index(cand)]

    # Heuristic fallback: pick a single numeric-only column (excluding obvious IDs / taxonomy)
    df = pd.read_csv(classifier_table, sep="\t", nrows=50)
    numeric_cols = []
    for c in df.columns:
        s = pd.to_numeric(df[c], errors="coerce")
        # treat as numeric if most values parse as numbers and at least one is non-null
        nonnull = s.notna().sum()
        if nonnull > 0 and nonnull / len(df) >= 0.8:
            numeric_cols.append(c)

    if len(numeric_cols) == 1:
        return numeric_cols[0]

    raise ValueError(
        f"Could not infer abundance column from: {cols}. "
        f"Expected one of: {candidates} (or a single numeric column). "
        f"Numeric-like columns detected: {numeric_cols}"
    )


def suggest_magnet_min_abundance(
    classifier_table: Path,
    target_frac: float = 0.99,
) -> Tuple[float, int, float]:
    """
    Suggest a MAGnet --min-abundance cutoff from a Lemur abundance table.

    Logic:
      - Take the abundance column (relative_abundance or similar).
      - Keep only positive abundances.
      - Sort descending, compute cumulative sum.
      - Choose the abundance where cumulative abundance first reaches
        target_frac * total_abundance.

    Args:
        classifier_table: Path to Lemur TSV output file
        target_frac: Target fraction of total abundance to capture (default: 0.99)

    Returns:
        Tuple of:
            cutoff (float): abundance threshold
            n_kept (int): number of taxa with abundance >= cutoff
            frac_explained (float): cumulative abundance fraction at that cutoff

    Raises:
        FileNotFoundError: If classifier table doesn't exist
        ValueError: If no abundance column found or no positive abundances
    """
    classifier_table = Path(classifier_table)
    if not classifier_table.exists():
        raise FileNotFoundError(f"Classifier table not found: {classifier_table}")

    # Infer which column contains abundance values
    col_name = infer_abundance_col_name(classifier_table)
    
    # Read the full table
    df = pd.read_csv(classifier_table, sep="\t")

    if col_name not in df.columns:
        raise ValueError(
            f"Abundance column '{col_name}' not found in {classifier_table}. "
            f"Columns: {list(df.columns)}"
        )

    # Extract abundances and filter to positive values only
    abundances = df[col_name].astype(float)
    abundances = abundances[abundances > 0]

    if abundances.empty:
        raise ValueError("No positive abundances found in classifier table.")

    # Sort descending and compute cumulative sum
    abund_sorted = abundances.sort_values(ascending=False).to_numpy()
    cum = abund_sorted.cumsum()
    total = cum[-1]

    # Find index where we first reach target_frac * total
    threshold = target_frac * total
    cutoff_idx = (cum >= threshold).argmax()
    cutoff = float(abund_sorted[cutoff_idx])

    n_kept = int(cutoff_idx + 1)
    frac_explained = float(cum[cutoff_idx] / total)

    return cutoff, n_kept, frac_explained


def main():
    parser = argparse.ArgumentParser(
        description="Suggest MAGnet --min-abundance cutoff from Lemur abundance table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Suggest cutoff to capture 99%% of abundance
  %(prog)s --input lemur_profile.tsv --output cutoff.txt

  # Capture 95%% of abundance (more permissive)
  %(prog)s --input lemur_profile.tsv --target-frac 0.95 --output cutoff.txt

  # Capture 99.9%% (stricter, more references)
  %(prog)s --input lemur_profile.tsv --target-frac 0.999 --output cutoff.txt
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        type=Path,
        required=True,
        help='Input Lemur abundance table (TSV format)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=Path,
        required=True,
        help='Output file containing just the cutoff value'
    )
    
    parser.add_argument(
        '--stats',
        type=Path,
        help='Optional output file for detailed statistics (TSV format)'
    )
    
    parser.add_argument(
        '--target-frac',
        type=float,
        default=0.99,
        help='Target fraction of total abundance to capture (default: 0.99)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print detailed information to stderr'
    )
    
    args = parser.parse_args()
    
    # Validate target fraction
    if not 0 < args.target_frac <= 1:
        parser.error(f"--target-frac must be between 0 and 1, got {args.target_frac}")
    
    try:
        # Calculate suggested cutoff
        cutoff, n_kept, frac_explained = suggest_magnet_min_abundance(
            args.input,
            args.target_frac
        )
        
        # Infer column name for reporting
        col_name = infer_abundance_col_name(args.input)
        
        # Read table to get total taxa count
        df = pd.read_csv(args.input, sep="\t")
        total_taxa = len(df[df[col_name] > 0])
        
        # Write cutoff value to output file (just the number for easy parsing)
        args.output.write_text(f"{cutoff}\n")
        
        # Write detailed stats if requested
        if args.stats:
            stats_lines = [
                f"abundance_column\t{col_name}",
                f"total_taxa_with_abundance\t{total_taxa}",
                f"target_fraction\t{args.target_frac:.4f}",
                f"suggested_cutoff\t{cutoff:.6g}",
                f"taxa_kept\t{n_kept}",
                f"taxa_dropped\t{total_taxa - n_kept}",
                f"fraction_explained\t{frac_explained:.6f}",
                f"fraction_explained_pct\t{frac_explained*100:.2f}",
            ]
            args.stats.write_text("\n".join(stats_lines) + "\n")
        
        # Print info to stderr if verbose
        if args.verbose:
            print("[INFO] Auto MAGnet cutoff from Lemur profile:", file=sys.stderr)
            print(f"       abundance column        : {col_name}", file=sys.stderr)
            print(f"       taxa with >0 abundance  : {total_taxa}", file=sys.stderr)
            print(f"       target fraction         : {args.target_frac:.3f}", file=sys.stderr)
            print(f"       suggested cutoff        : {cutoff:.6g}", file=sys.stderr)
            print(f"       taxa kept               : {n_kept}", file=sys.stderr)
            print(f"       taxa dropped            : {total_taxa - n_kept}", file=sys.stderr)
            print(f"       fraction explained      : {frac_explained*100:.2f}%", file=sys.stderr)
        
        return 0
        
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
