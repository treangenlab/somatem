#!/usr/bin/env python3
"""
Select Variant Caller Based on Coverage

Analyzes depth coverage from samtools depth output and decides whether
to use Clair3 (high coverage ≥8x) or bcftools call (low coverage <8x).

Based on marlin.py compute_depth() and variant calling logic.
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Select variant caller based on coverage depth"
    )
    parser.add_argument(
        "--depth",
        required=True,
        type=Path,
        help="Input depth file from samtools depth",
    )
    parser.add_argument(
        "--min-cov",
        type=float,
        default=8.0,
        help="Minimum coverage threshold for Clair3 (default: 8.0)",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output file to write caller choice (clair3 or bcftools)",
    )
    parser.add_argument(
        "--stats",
        required=True,
        type=Path,
        help="Output file to write coverage statistics",
    )
    return parser.parse_args()


def calculate_coverage_stats(depth_file):
    """
    Calculate coverage statistics from samtools depth output.
    
    Args:
        depth_file: Path to samtools depth output file
        
    Returns:
        dict: Coverage statistics including mean, median, and covered bases
    """
    depths = []
    total_bases = 0
    covered_bases = 0
    
    try:
        with open(depth_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                    
                try:
                    depth = int(parts[2])
                    depths.append(depth)
                    total_bases += 1
                    if depth > 0:
                        covered_bases += 1
                except (ValueError, IndexError):
                    continue
                    
    except FileNotFoundError:
        print(f"Error: Depth file not found: {depth_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading depth file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if not depths:
        print("Error: No valid depth data found in input file", file=sys.stderr)
        sys.exit(1)
    
    # Calculate statistics
    mean_depth = sum(depths) / len(depths) if depths else 0
    median_depth = sorted(depths)[len(depths) // 2] if depths else 0
    max_depth = max(depths) if depths else 0
    coverage_fraction = covered_bases / total_bases if total_bases > 0 else 0
    
    return {
        "mean_depth": mean_depth,
        "median_depth": median_depth,
        "max_depth": max_depth,
        "total_bases": total_bases,
        "covered_bases": covered_bases,
        "coverage_fraction": coverage_fraction,
    }


def select_caller(stats, min_coverage):
    """
    Select variant caller based on coverage statistics.
    
    Args:
        stats: Coverage statistics dictionary
        min_coverage: Minimum coverage threshold for Clair3
        
    Returns:
        str: Caller choice ('clair3' or 'bcftools')
    """
    mean_depth = stats["mean_depth"]
    
    if mean_depth >= min_coverage:
        return "clair3"
    else:
        return "bcftools"


def write_output(caller, output_file):
    """Write caller choice to output file."""
    try:
        with open(output_file, 'w') as f:
            f.write(f"{caller}\n")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)


def write_stats(stats, caller, min_cov, stats_file):
    """Write coverage statistics to output file."""
    try:
        with open(stats_file, 'w') as f:
            f.write("metric\tvalue\n")
            f.write(f"mean_depth\t{stats['mean_depth']:.2f}\n")
            f.write(f"median_depth\t{stats['median_depth']:.2f}\n")
            f.write(f"max_depth\t{stats['max_depth']}\n")
            f.write(f"total_bases\t{stats['total_bases']}\n")
            f.write(f"covered_bases\t{stats['covered_bases']}\n")
            f.write(f"coverage_fraction\t{stats['coverage_fraction']:.4f}\n")
            f.write(f"min_cov_threshold\t{min_cov}\n")
            f.write(f"selected_caller\t{caller}\n")
    except Exception as e:
        print(f"Error writing stats file: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """Main function."""
    args = parse_args()
    
    # Calculate coverage statistics
    stats = calculate_coverage_stats(args.depth)
    
    # Select caller based on coverage
    caller = select_caller(stats, args.min_cov)
    
    # Write outputs
    write_output(caller, args.output)
    write_stats(stats, caller, args.min_cov, args.stats)
    
    # Print summary to stderr for logging
    print(f"Coverage analysis complete:", file=sys.stderr)
    print(f"  Mean depth: {stats['mean_depth']:.2f}x", file=sys.stderr)
    print(f"  Median depth: {stats['median_depth']:.2f}x", file=sys.stderr)
    print(f"  Coverage: {stats['coverage_fraction']*100:.2f}%", file=sys.stderr)
    print(f"  Selected caller: {caller}", file=sys.stderr)


if __name__ == "__main__":
    main()
