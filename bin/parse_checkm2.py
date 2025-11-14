#!/usr/bin/env python3

# Script to parse CheckM2 output and create a mapping of bin names to completeness
# Usage: python parse_checkm2.py <checkm2_tsv>

import pandas as pd
import csv
import os
import sys

if len(sys.argv) < 2:
    print("Usage: python parse_checkm2.py <checkm2_tsv>")
    sys.exit(1)

checkm2_tsv = sys.argv[1]

# Read CheckM2 results
checkm2_df = pd.read_csv(checkm2_tsv, sep='\t')

# Create a mapping of bin name to completeness
completeness_map = {}
for _, row in checkm2_df.iterrows():
    bin_name = row['Name']
    completeness = float(row['Completeness'])
    completeness_map[bin_name] = completeness
    print(f"Bin {bin_name}: {completeness}% complete")

# Write the completeness mapping to a CSV file
with open('bins_with_completeness.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['bin_name', 'completeness'])
    for bin_name, completeness in completeness_map.items():
        writer.writerow([bin_name, completeness])

print(f"Created completeness mapping for {len(completeness_map)} bins")