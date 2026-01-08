#!/usr/bin/env python3
"""
Parse samtools depth output and compute statistics.
"""
import argparse
import json
from collections import defaultdict
from statistics import median
from typing import Dict, List

def parse_depth(depth_file: str) -> tuple:
    """
    Parse samtools depth -a output.
    
    Returns:
        tuple: (depths_dict, chrom_meds_dict, global_median)
    """
    depths: Dict[str, List[int]] = defaultdict(list)
    
    with open(depth_file) as f:
        for line in f:
            chrom, pos_str, d_str = line.rstrip().split("\t")
            d = int(d_str)
            depths[chrom].append(d)
    
    chrom_meds: Dict[str, float] = {}
    chrom_med_list: List[float] = []
    
    for chrom, arr in depths.items():
        if not arr:
            chrom_meds[chrom] = 0.0
        else:
            m = float(median(arr))
            chrom_meds[chrom] = m
            chrom_med_list.append(m)
    
    global_median = float(median(chrom_med_list)) if chrom_med_list else 0.0
    
    return depths, chrom_meds, global_median

def find_coverage_gaps(
    depths: Dict[str, List[int]],
    min_cov: int,
    min_gap_len: int,
) -> Dict[str, List[tuple]]:
    """
    Identify coverage gaps.
    
    Returns:
        dict: chrom -> list of (start, end) 1-based inclusive coordinates
    """
    gaps: Dict[str, List[tuple]] = {}
    
    for chrom, arr in depths.items():
        chrom_gaps: List[tuple] = []
        in_gap = False
        gap_start = 0
        
        for i, d in enumerate(arr):
            pos = i + 1
            if d < min_cov:
                if not in_gap:
                    in_gap = True
                    gap_start = pos
            else:
                if in_gap:
                    gap_end = pos - 1
                    if gap_end - gap_start + 1 >= min_gap_len:
                        chrom_gaps.append((gap_start, gap_end))
                    in_gap = False
        
        if in_gap:
            gap_end = len(arr)
            if gap_end - gap_start + 1 >= min_gap_len:
                chrom_gaps.append((gap_start, gap_end))
        
        gaps[chrom] = chrom_gaps
    
    return gaps

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--depth', required=True, help='Depth file from samtools depth')
    parser.add_argument('--min-cov', type=int, required=True, help='Minimum coverage for gap detection')
    parser.add_argument('--min-gap-len', type=int, required=True, help='Minimum gap length')
    parser.add_argument('--output-json', required=True, help='Output JSON file')
    args = parser.parse_args()
    
    depths, chrom_meds, global_med = parse_depth(args.depth)
    gaps = find_coverage_gaps(depths, args.min_cov, args.min_gap_len)
    
    # Convert depths to serializable format
    depths_serializable = {k: v for k, v in depths.items()}
    
    result = {
        'global_median': global_med,
        'chrom_medians': chrom_meds,
        'gaps': {k: [(s, e) for s, e in v] for k, v in gaps.items()},
        'depths': depths_serializable
    }
    
    with open(args.output_json, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"Global median depth: {global_med:.2f}")
    for chrom, m in chrom_meds.items():
        print(f"  {chrom}: {m:.2f}")

if __name__ == '__main__':
    main()
