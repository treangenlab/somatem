#!/usr/bin/env python3
"""
Build contigs from consensus by breaking at coverage gaps.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


def parse_depth(depth_file: Path) -> Dict[str, List[int]]:
    """Parse samtools depth output."""
    depths: Dict[str, List[int]] = defaultdict(list)
    
    with depth_file.open() as f:
        for line in f:
            parts = line.rstrip().split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            depth_val = int(parts[2])
            depths[chrom].append(depth_val)
    
    return depths


def find_coverage_gaps(
    depths: Dict[str, List[int]],
    min_cov: int,
    min_gap_len: int,
) -> Dict[str, List[Tuple[int, int]]]:
    """Find coverage gaps (runs of depth < min_cov with length >= min_gap_len)."""
    gaps: Dict[str, List[Tuple[int, int]]] = {}
    
    for chrom, arr in depths.items():
        chrom_gaps: List[Tuple[int, int]] = []
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


def load_fasta(path: Path) -> Dict[str, str]:
    """Load FASTA sequences."""
    seqs: Dict[str, List[str]] = {}
    name = None
    buf: List[str] = []
    
    with path.open() as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = ''.join(buf)
    
    return seqs


def write_fasta(seqs: Dict[str, str], path: Path):
    """Write FASTA sequences."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as f:
        for name, seq in seqs.items():
            f.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def build_contigs(
    consensus_fa: Path,
    depth_file: Path,
    min_cov: int,
    min_gap_len: int,
    min_contig_len: int,
    contigs_out: Path,
    layout_out: Path,
):
    """Build contigs from consensus by breaking at coverage gaps."""
    
    # Parse depth and find gaps
    depths = parse_depth(depth_file)
    gaps = find_coverage_gaps(depths, min_cov, min_gap_len)
    
    print(f"[INFO] Found gaps in sequences:", file=sys.stderr)
    for chrom, gap_list in gaps.items():
        print(f"  {chrom}: {len(gap_list)} gaps", file=sys.stderr)
    
    # Load consensus
    seqs = load_fasta(consensus_fa)
    
    # Build contigs
    contigs: Dict[str, str] = {}
    rows: List[str] = []
    
    for chrom, seq in seqs.items():
        L = len(seq)
        chrom_gaps = gaps.get(chrom, [])
        
        # Build intervals between gaps
        intervals: List[Tuple[int, int]] = []
        last_end = 0
        
        for (gstart, gend) in chrom_gaps:
            start = last_end + 1
            end = gstart - 1
            if end >= start:
                intervals.append((start, end))
            last_end = gend
        
        if last_end < L:
            intervals.append((last_end + 1, L))
        
        # Extract contigs
        for idx, (start, end) in enumerate(intervals, start=1):
            length = end - start + 1
            if length < min_contig_len:
                continue
            
            contig_id = f"{chrom}_seg{idx}:{start}-{end}"
            subseq = seq[start - 1:end]
            
            # Compute mean depth
            mean_depth_val = 0.0
            if chrom in depths:
                arr = depths[chrom]
                s_idx = max(0, start - 1)
                e_idx = min(len(arr), end)
                subd = arr[s_idx:e_idx]
                if subd:
                    mean_depth_val = float(sum(subd)) / len(subd)
            
            contigs[contig_id] = subseq
            rows.append(
                f"{contig_id}\t{chrom}\t{start}\t{end}\t{length}\t{mean_depth_val:.3f}"
            )
    
    # Write outputs
    write_fasta(contigs, contigs_out)
    
    layout_out.parent.mkdir(parents=True, exist_ok=True)
    with layout_out.open('w') as f:
        f.write("contig_id\tchrom\tstart\tend\tlength\tmean_depth\n")
        for row in rows:
            f.write(row + '\n')
    
    print(f"[INFO] Wrote {len(contigs)} contigs to {contigs_out}", file=sys.stderr)
    print(f"[INFO] Wrote layout to {layout_out}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Build contigs from consensus by breaking at coverage gaps'
    )
    parser.add_argument('--consensus', required=True, help='Input consensus FASTA')
    parser.add_argument('--depth', required=True, help='Depth file from samtools depth')
    parser.add_argument('--min-cov', type=int, default=2, help='Minimum coverage for non-gap')
    parser.add_argument('--min-gap-len', type=int, default=200, help='Minimum gap length')
    parser.add_argument('--min-contig-len', type=int, default=1000, help='Minimum contig length')
    parser.add_argument('--contigs-out', required=True, help='Output contigs FASTA')
    parser.add_argument('--layout-out', required=True, help='Output layout TSV')
    
    args = parser.parse_args()
    
    build_contigs(
        consensus_fa=Path(args.consensus),
        depth_file=Path(args.depth),
        min_cov=args.min_cov,
        min_gap_len=args.min_gap_len,
        min_contig_len=args.min_contig_len,
        contigs_out=Path(args.contigs_out),
        layout_out=Path(args.layout_out),
    )


if __name__ == '__main__':
    main()
