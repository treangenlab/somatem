#!/usr/bin/env python
import csv
import sys
import argparse
from collections import namedtuple, defaultdict
csv.field_size_limit(sys.maxsize)

def main():
    parser = argparse.ArgumentParser(description="Get threats from bowtie2 and rapsearch")
    parser.add_argument("--rapsearch-results", required=True)
    parser.add_argument("--bowtie2-results", required=True)
    parser.add_argument("-o", "--out", help="Path to output file")
    args = parser.parse_args()

    rapsearch_file = args.rapsearch_results
    bowtie_file = args.bowtie2_results
    output_file = args.out

    bsat_hits = set()
    # for filename in [rapsearch_file, bowtie_file]:
    #     with open(filename) as fd:
    #         for line in fd:
    #             if line[0] == "#":
    #                 continue
    #             else:
    #                 bsat_hits.add(line.split('\t', 1)[0])
                    
    with open(bowtie_file) as fd:
        for line in fd:
            if line[0] == "#":
                continue
            else:
                bsat_hits.add(line.split('\t', 1)[0])
    
    with open(rapsearch_file) as fd:
        for line in fd:
            if line[0] == "#":
                continue
            else:
                vals = line.split('\t')
                identity = float(vals[2])
                length = float(vals[3])
                
                if identity >= 80 and length >= 30:
                    bsat_hits.add(vals[0])

    with open(output_file, 'w') as out_fd:
        out_fd.writelines(line + '\n' for line in ["query"] + list(bsat_hits))

main()

