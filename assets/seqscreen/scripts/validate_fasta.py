#!/usr/bin/env python3

import argparse
import gzip
import sys


NT_BASES = set("ATGCWSMKRYBDHVN")
AA_BASES = set("ARNDCQEGHILKMFPSTWYVX")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate a nucleotide or peptide FASTA file."
    )
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-co", "--cleaned_output")
    parser.add_argument(
        "-mx",
        "--max_seq_size",
        type=int,
        default=1_000_000_000,
    )
    parser.add_argument("-mn", "--max_N", type=int, default=1_000_000_000)
    parser.add_argument("-p", "--peptide", action="store_true")
    args = parser.parse_args()

    if args.max_seq_size <= 0 or args.max_seq_size > 1_000_000_000:
        parser.error("--max_seq_size must be an integer > 0 and <= 1 billion")
    if args.max_N <= 0 or args.max_N > 1_000_000_000:
        parser.error("--max_N must be an integer > 0 and <= 1 billion")
    return args


def open_fasta(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def fail(message):
    raise SystemExit(f"\n FAIL: {message}\n")


def evaluate_seq(header, seq, args, seen, counters, out_handle):
    seq = seq.upper()
    keep = True

    if len(seq) == 0:
        if args.cleaned_output:
            keep = False
        else:
            fail(f"Your input FASTA file contains empty sequences (no base pairs): {header}")
        counters["no_seq"] += 1

    if len(seq) > args.max_seq_size:
        if args.cleaned_output:
            keep = False
        else:
            fail(
                "Your input FASTA file contains sequences larger than your "
                f"cut-off of {args.max_seq_size} : {header}"
            )
        counters["seqs_too_big"] += 1

    if args.peptide:
        invalid = "".join(base for base in seq if base not in AA_BASES)
    else:
        number_of_ns = seq.count("N")
        if number_of_ns >= args.max_N:
            if args.cleaned_output:
                keep = False
            else:
                fail(
                    "Your input FASTA file contains a sequence with too many "
                    f"'N' bases.\n Limit = {args.max_N}\n Actual = {number_of_ns}\n Seq: {header}"
                )
            counters["too_many_N"] += 1
        invalid = "".join(base for base in seq if base not in NT_BASES)

    if invalid:
        if args.cleaned_output:
            keep = False
        else:
            fail(
                "Your input FASTA file contains invalid bases (valid NT bases: "
                "A,T,G,C,W,S,M,K,R,Y,B,D,H,V,N; valid AA residues: "
                "A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,X) the following "
                f"sequence contains the characters '{invalid}':\n\n{header}"
            )
        counters["bad_bases"] += 1

    if header in seen:
        if args.cleaned_output:
            keep = False
        else:
            fail(
                "Your input FASTA file does not contain unique identifiers. "
                f"Offending sequence ID:\n\n{header}"
            )
        counters["dup_header"] += 1
    seen.add(header)

    if args.cleaned_output and keep:
        out_handle.write(f"{header}\n{seq}\n")


def main():
    args = parse_args()
    counters = {
        "no_seq": 0,
        "dup_header": 0,
        "bad_bases": 0,
        "seqs_too_big": 0,
        "too_many_N": 0,
    }
    seen = set()
    out_handle = open(args.cleaned_output, "w") if args.cleaned_output else None

    try:
        with open_fasta(args.fasta) as handle:
            header = ""
            seq_parts = []
            for raw_line in handle:
                line = raw_line.replace("\r", "\n").strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        evaluate_seq(header, "".join(seq_parts), args, seen, counters, out_handle)
                    header = line
                    seq_parts = []
                else:
                    seq_parts.append(line)

            evaluate_seq(header, "".join(seq_parts), args, seen, counters, out_handle)
    finally:
        if out_handle:
            out_handle.close()

    if args.cleaned_output and any(counters.values()):
        sys.stderr.write(
            "\n WARNING: Cleaned output written to: "
            f"{args.cleaned_output}\n However, some sequences were thrown out "
            "for the following reasons:\n"
            f" Duplicate headers = {counters['dup_header']}\n"
            f" Header with no sequence data = {counters['no_seq']}\n"
            f" Sequence contains invalid bases = {counters['bad_bases']}\n"
            f" Sequences were too big: {counters['seqs_too_big']}\n"
            f" Sequences with too many 'N's: {counters['too_many_N']}\n\n"
        )


if __name__ == "__main__":
    main()
