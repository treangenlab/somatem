#!/usr/bin/env python3
"""
Classify variants by confidence level based on coverage and allele frequency.
Outputs classified VCF and high-confidence VCF for consensus generation.
"""

import sys
import gzip
import argparse
import subprocess
from statistics import median
from collections import defaultdict


def parse_depth_file(depth_file):
    """Parse samtools depth output and calculate coverage statistics."""
    depths = defaultdict(list)
    
    with open(depth_file) as f:
        for line in f:
            chrom, pos_str, d_str = line.rstrip().split("\t")
            depths[chrom].append(int(d_str))
    
    return depths


def calculate_thresholds(depths, high_cov_frac):
    """Calculate low and high coverage thresholds."""
    chrom_meds = []
    for chrom, arr in depths.items():
        if arr:
            chrom_meds.append(float(median(arr)))
    
    global_median_cov = float(median(chrom_meds)) if chrom_meds else 0.0
    
    # Calculate thresholds
    low_cov = max(2, int(global_median_cov * 0.1))
    high_cov = max(10, int(global_median_cov * high_cov_frac))
    
    print(f"[INFO] Global median depth: {global_median_cov:.2f}", file=sys.stderr)
    print(f"[INFO] Low coverage threshold: {low_cov}", file=sys.stderr)
    print(f"[INFO] High coverage threshold: {high_cov}", file=sys.stderr)
    
    return low_cov, high_cov, global_median_cov


def open_vcf(vcf_file):
    """Open VCF file (handles both gzipped and plain text)."""
    if vcf_file.endswith('.gz'):
        return gzip.open(vcf_file, 'rt')
    else:
        return open(vcf_file, 'r')


def classify_variant(site_cov, af, low_cov, high_cov, af_high, af_med):
    """Classify a variant based on coverage and allele frequency."""
    if site_cov >= high_cov and af >= af_high:
        return "HIGH"
    elif site_cov >= low_cov and af >= af_med:
        return "MED"
    else:
        return "LOW"


def parse_info_field(info):
    """Parse VCF INFO field into a dictionary."""
    info_dict = {}
    if info and info != ".":
        for item in info.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                info_dict[k] = v
    return info_dict


def get_allele_frequency(info_dict):
    """Extract allele frequency from INFO field."""
    af = 1.0
    if "AF" in info_dict:
        try:
            af_vals = [float(x) for x in info_dict["AF"].split(",")]
            af = max(af_vals)
        except ValueError:
            pass
    return af


def process_vcf(vcf_file, depth_file, prefix, af_high, af_med, high_cov_frac):
    """Main processing function."""
    
    # Parse depth file
    depths = parse_depth_file(depth_file)
    
    # Calculate thresholds
    low_cov, high_cov, global_median_cov = calculate_thresholds(depths, high_cov_frac)
    
    # Read VCF
    header_lines = []
    data_lines = []
    
    with open_vcf(vcf_file) as fin:
        for line in fin:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                data_lines.append(line)
    
    # Check if CONF header already exists
    has_conf = any("ID=CONF" in h for h in header_lines if h.startswith("##INFO"))
    
    # Write classified VCFs
    classified_vcf = f"{prefix}_classified.vcf"
    highconf_vcf = f"{prefix}_highconf.vcf"
    
    with open(classified_vcf, "w") as fout_all, \
         open(highconf_vcf, "w") as fout_high:
        
        # Write headers
        for line in header_lines:
            if line.startswith("#CHROM") and not has_conf:
                # Add CONF INFO header before #CHROM line
                fout_all.write('##INFO=<ID=CONF,Number=1,Type=String,Description="Variant confidence: HIGH|MED|LOW">\n')
                fout_high.write('##INFO=<ID=CONF,Number=1,Type=String,Description="Variant confidence: HIGH|MED|LOW">\n')
            fout_all.write(line)
            fout_high.write(line)
        
        # Process variants
        for line in data_lines:
            cols = line.rstrip().split("\t")
            if len(cols) < 8:
                continue
            
            chrom, pos_str, vid, ref, alt, qual, flt, info = cols[:8]
            pos = int(pos_str)
            
            # Parse INFO
            info_dict = parse_info_field(info)
            
            # Get allele frequency
            af = get_allele_frequency(info_dict)
            
            # Get site coverage
            site_cov = 0
            if chrom in depths:
                arr = depths[chrom]
                idx = pos - 1
                if 0 <= idx < len(arr):
                    site_cov = arr[idx]
            
            # Classify variant
            conf = classify_variant(site_cov, af, low_cov, high_cov, af_high, af_med)
            
            # Update INFO field
            if info == "." or info == "":
                info_str = f"CONF={conf}"
            else:
                info_str = f"{info};CONF={conf}"
            
            cols[7] = info_str
            new_line = "\t".join(cols) + "\n"
            
            # Write to classified VCF
            fout_all.write(new_line)
            
            # Write to high-confidence VCF if HIGH
            if conf == "HIGH":
                fout_high.write(new_line)
    
    # Compress and index both VCFs
    for vcf_name in [classified_vcf, highconf_vcf]:
        vcf_gz = f"{vcf_name}.gz"
        
        # Sort and compress
        subprocess.run([
            'bcftools', 'sort', '-Oz',
            '-o', vcf_gz,
            vcf_name
        ], check=True)
        
        # Index
        subprocess.run([
            'tabix', '-p', 'vcf',
            vcf_gz
        ], check=True)
        
        print(f"[INFO] Created {vcf_gz}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Classify variants by confidence')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--depth', required=True, help='Samtools depth file')
    parser.add_argument('--prefix', required=True, help='Output prefix')
    parser.add_argument('--af-high', type=float, default=0.9, help='AF threshold for HIGH confidence')
    parser.add_argument('--af-med', type=float, default=0.7, help='AF threshold for MED confidence')
    parser.add_argument('--high-cov-frac', type=float, default=0.5, help='Fraction of median for high coverage')
    
    args = parser.parse_args()
    
    process_vcf(
        args.vcf,
        args.depth,
        args.prefix,
        args.af_high,
        args.af_med,
        args.high_cov_frac
    )


if __name__ == '__main__':
    main()
