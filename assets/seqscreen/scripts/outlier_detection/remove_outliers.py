#!/usr/bin/env python3

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply SeqScreen outlier detection calls to a BLAST BTAB file."
    )
    parser.add_argument("-ol", "--outliers", required=True)
    parser.add_argument("-b", "--btab", required=True)
    parser.add_argument("-o", "--out", required=True)
    return parser.parse_args()


def load_outliers(path):
    outliers = {}
    with open(path, "r") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields or fields[0] == "#query_sequence":
                continue
            if len(fields) > 1 and fields[1] != "NA":
                outliers.setdefault(fields[0], set()).update(fields[1].split(";"))
    return outliers


def main():
    args = parse_args()
    outliers = load_outliers(args.outliers)

    with open(args.out, "w") as out_handle, open(args.btab, "r") as btab_handle:
        for line in btab_handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in outliers:
                if len(fields) > 1 and fields[1] in outliers[fields[0]]:
                    out_handle.write(line)
            else:
                out_handle.write(line)


if __name__ == "__main__":
    main()
