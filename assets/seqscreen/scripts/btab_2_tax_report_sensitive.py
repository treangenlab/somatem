#!/usr/bin/env python3

import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert sensitive-mode BLASTx/BLASTn BTAB files to SeqScreen taxonomy format."
    )
    parser.add_argument("-x", "--blastx", required=True)
    parser.add_argument("-n", "--blastn", required=True)
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("-c", "--cutoff", type=int, default=0)
    args = parser.parse_args()
    if args.cutoff < 0 or args.cutoff > 50:
        parser.error(f"--cutoff must be between 0 and 50, not {args.cutoff}")
    return args


def base_query(query):
    return re.sub(r"_unambig_\d+$", "", query)


def find_top_hits(path, top_hits):
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            query = base_query(fields[0])
            bitscore = float(fields[11])
            if bitscore < 0:
                raise ValueError(
                    "Error: parsing this BTAB file did not go as expected. "
                    "Might be grabbing the wrong field for Bit Score."
                )
            if query not in top_hits or top_hits[query] < bitscore:
                top_hits[query] = bitscore


def add_result(results, sources, confidence, query, taxid, source, bitscore, percent_identity):
    if query in results and taxid in results[query]:
        if results[query][taxid] < bitscore:
            results[query][taxid] = bitscore
    else:
        results.setdefault(query, {})[taxid] = bitscore

    sources.setdefault(query, {})[taxid] = source
    confidence.setdefault(query, {})[taxid] = percent_identity / 100.0


def collect_blastn(path, cutoff, top_hits, results, sources, confidence):
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            query = base_query(fields[0])
            bitscore = float(fields[11])
            taxid = fields[14]
            if bitscore >= cutoff * top_hits[query]:
                add_result(
                    results,
                    sources,
                    confidence,
                    query,
                    taxid,
                    "blastn",
                    bitscore,
                    float(fields[2]),
                )


def collect_blastx(path, cutoff, top_hits, results, sources, confidence):
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            query = base_query(fields[0])
            bitscore = float(fields[11])
            taxid = re.sub(r"^.*TaxID=", "", fields[15])
            taxid = re.sub(r" .*$", "", taxid)
            if bitscore >= cutoff * top_hits[query]:
                source = "all" if query in results and taxid in results[query] else "blastx"
                add_result(
                    results,
                    sources,
                    confidence,
                    query,
                    taxid,
                    source,
                    bitscore,
                    float(fields[2]),
                )


def main():
    args = parse_args()
    cutoff = (100 - int(args.cutoff)) / 100
    top_hits = {}
    results = {}
    sources = {}
    confidence = {}

    find_top_hits(args.blastx, top_hits)
    find_top_hits(args.blastn, top_hits)
    collect_blastn(args.blastn, cutoff, top_hits, results, sources, confidence)
    collect_blastx(args.blastx, cutoff, top_hits, results, sources, confidence)

    with open(args.out, "w") as out:
        out.write("#query\ttaxid\tsource\tconfidence\n")
        for query in sorted(results):
            ranked_taxids = sorted(
                results[query],
                key=lambda taxid: results[query][taxid],
                reverse=True,
            )
            out.write(
                f"{query}\t"
                f"{','.join(ranked_taxids)}\t"
                f"{','.join(sources[query][taxid] for taxid in ranked_taxids)}\t"
                f"{','.join(str(confidence[query][taxid]) for taxid in ranked_taxids)}\n"
            )


if __name__ == "__main__":
    main()
