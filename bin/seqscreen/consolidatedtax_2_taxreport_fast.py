#!/usr/bin/env python3
import argparse, sys, csv, os, glob

top_hits = {}
top_hits_centrifuge = {}
centrifuge_results_dict = {}
diamond_results_dict = {}

def base_query(id_):
    return id_.split("_unambig_")[0]

def resolve_input(pathlike, label, debug=False):
    """Return a concrete file path. If a directory is given, try to find a single plausible file inside."""
    if os.path.isfile(pathlike):
        if debug: print(f"[debug] Using {label} file: {pathlike}", file=sys.stderr)
        return pathlike
    if os.path.isdir(pathlike):
        patterns = ("*.btab","*.m8","*.tsv","*.txt")
        matches = []
        for pat in patterns:
            matches.extend(glob.glob(os.path.join(pathlike, pat)))
        matches = sorted(set(matches))
        if len(matches) == 0:
            sys.exit(f"Error: {label} path is a directory and no files matched {patterns}: {pathlike}")
        if len(matches) > 1:
            sys.exit(f"Error: {label} path is a directory with multiple candidate files: {matches}\n"
                     f"Please pass the exact file via --{label.lower()}.")
        if debug: print(f"[debug] Resolved {label} file in directory: {matches[0]}", file=sys.stderr)
        return matches[0]
    sys.exit(f"Error: {label} path does not exist: {pathlike}")

def findTopHitsDiamond(diamond_f):
    with open(diamond_f,'r') as inpf:
        for i,line in enumerate(inpf,1):
            if not line.strip() or line.startswith('#'): 
                continue
            tokens = line.rstrip('\n').split('\t')
            if len(tokens) < 12:
                sys.exit(f"Error: DIAMOND line {i} has <12 columns. Check format or output fields.")
            query = base_query(tokens[0])
            try:
                bitscore = float(tokens[11])  # m8 col 12
            except ValueError:
                sys.exit(f"Error: DIAMOND bitscore not a float on line {i}: {tokens[11]}")
            prev = top_hits.get(query)
            if prev is None or bitscore > prev:
                top_hits[query] = bitscore

def findTopHitsCentrifuge(centrifuge_f):
    with open(centrifuge_f,"r") as inpf:
        inpf_reader = csv.DictReader(inpf, delimiter="\t")
        required = {"readID","score"}
        if not required.issubset(inpf_reader.fieldnames or []):
            sys.exit(f"Error: Centrifuge header missing required fields {required}. Got: {inpf_reader.fieldnames}")
        for row in inpf_reader:
            query = base_query(row["readID"])
            try:
                score = float(row["score"])
            except ValueError:
                # allow empty or unparseable scores to be skipped
                continue
            prev = top_hits_centrifuge.get(query)
            if prev is None or score > prev:
                top_hits_centrifuge[query] = score

def processCentrifuge(centrifuge_f, taxlimit):
    with open(centrifuge_f,'r') as inpf:
        inpf_reader = csv.DictReader(inpf, delimiter="\t")
        required = {"readID","taxID","score","seqID","hitLength","queryLength"}
        if not required.issubset(inpf_reader.fieldnames or []):
            sys.exit(f"Error: Centrifuge header missing required fields {required}. Got: {inpf_reader.fieldnames}")
        for row in inpf_reader:
            query = base_query(row["readID"])
            if not query:
                continue
            if query in centrifuge_results_dict and len(centrifuge_results_dict[query][0]) >= taxlimit:
                continue
            taxid = row["taxID"]
            if taxid == "0":
                continue
            try:
                score = float(row["score"])
                hit_length = int(row["hitLength"])
                query_length = int(row["queryLength"])
            except ValueError:
                continue
            centrifuge_conf = hit_length / query_length if query_length > 0 else 0.0
            top = top_hits_centrifuge.get(query)
            if top is not None and score == top:
                if query in centrifuge_results_dict:
                    if taxid not in centrifuge_results_dict[query][0]:
                        centrifuge_results_dict[query][0].add(taxid)
                        centrifuge_results_dict[query][1].append(score)
                        centrifuge_results_dict[query][2].append(f"{centrifuge_conf:.6g}")
                        centrifuge_results_dict[query][3].append("centrifuge")
                else:
                    centrifuge_results_dict[query] = (set([taxid]), [score], [f"{centrifuge_conf:.6g}"], ["centrifuge"])

def processDiamond(diamond_f, cutoff, taxlimit):
    with open(diamond_f,'r') as inpf:
        for i,line in enumerate(inpf,1):
            if not line.strip() or line.startswith('#'):
                continue
            tokens = line.rstrip('\n').split('\t')
            if len(tokens) < 12:
                # Not m8: bail with message that shows the first bad line
                sys.exit(f"Error: DIAMOND line {i} has <12 columns. Check format or output fields.")
            query = base_query(tokens[0])
            try:
                bitscore = float(tokens[11])
            except ValueError:
                sys.exit(f"Error: DIAMOND bitscore not a float on line {i}: {tokens[11]}")
            # taxid extraction: robust fallback if extra columns missing
            taxid = None
            if len(tokens) > 15 and "ID=" in tokens[15]:
                taxid = tokens[15].split("ID=",1)[1].split(" ",1)[0]
            else:
                # If your .btab stores taxid elsewhere, set it here or fail clearly
                sys.exit(f"Error: Cannot parse taxid on DIAMOND line {i}. Expected 'ID=' in column 16.")
            # cutoff guard
            top = top_hits.get(query)
            if top is None:
                # If we never saw this query when scanning top hits, skip safely
                continue
            bitscore_cutoff = cutoff * top
            if query in diamond_results_dict and len(diamond_results_dict[query][0]) >= taxlimit:
                continue
            if bitscore > bitscore_cutoff:
                if query in diamond_results_dict:
                    if taxid not in diamond_results_dict[query][0]:
                        diamond_results_dict[query][0].add(taxid)
                        diamond_results_dict[query][1].append(bitscore)
                        # tokens[2] is percent identity in m8
                        pid = float(tokens[2]) if tokens[2] else 0.0
                        diamond_results_dict[query][2].append(f"{pid/100:.6g}")
                        diamond_results_dict[query][3].append("diamond")
                else:
                    pid = float(tokens[2]) if tokens[2] else 0.0
                    diamond_results_dict[query] = (set([taxid]), [bitscore], [f"{pid/100:.6g}"], ["diamond"])

def writeResults(out_f):
    with open(out_f,'w') as outf:
        outf.write("#query\ttaxid\tsource\tconfidence\tcentrifuge_top_priority\tcombined_taxid\tcentrifuge_multi_tax\tdiamond_multi_tax\n")
        keys = set(centrifuge_results_dict.keys()) | set(diamond_results_dict.keys())
        for k in sorted(keys):
            if k in centrifuge_results_dict:
                taxa = sorted(centrifuge_results_dict[k][0])
                source = centrifuge_results_dict[k][3]
                confidence = centrifuge_results_dict[k][2]
                top1 = [taxa[0]] if taxa else []
            else:
                taxa = sorted(diamond_results_dict[k][0])
                source = diamond_results_dict[k][3]
                confidence = diamond_results_dict[k][2]
                top1 = taxa
            c_multi = ",".join(sorted(centrifuge_results_dict.get(k, (set(),))[0])) if k in centrifuge_results_dict else ""
            d_multi = ",".join(sorted(diamond_results_dict.get(k, (set(),))[0])) if k in diamond_results_dict else ""
            combined = ",".join(sorted(set(centrifuge_results_dict.get(k, (set(),))[0]) |
                                       set(diamond_results_dict.get(k, (set(),))[0])))
            outf.write(
                f"{k}\t{','.join(taxa)}\t{','.join(source)}\t{','.join(confidence)}\t"
                f"{','.join(top1)}\t{combined}\t{c_multi}\t{d_multi}\n"
            )

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-d","--diamond", required=True, help="Path to DIAMOND .btab/.m8 file OR a directory containing it")
    p.add_argument("-c","--centrifuge", required=True, help="Path to Centrifuge TSV file OR a directory containing it")
    p.add_argument("-o","--out", required=True, help="Output TSV")
    p.add_argument("-t","--cutoff", type=int, default=1, help="Top-N percent cutoff from max bitscore (1..100)")
    p.add_argument("--taxlimit", type=int, default=25, help="Max number of taxids to report per query")
    p.add_argument("--debug", action="store_true")
    args = p.parse_args()

    cutoff = (100 - args.cutoff) / 100.0
    diamond_f = resolve_input(args.diamond, "DIAMOND", args.debug)
    centrifuge_f = resolve_input(args.centrifuge, "CENTRIFUGE", args.debug)

    findTopHitsCentrifuge(centrifuge_f)
    findTopHitsDiamond(diamond_f)
    processCentrifuge(centrifuge_f, args.taxlimit)
    processDiamond(diamond_f, cutoff, args.taxlimit)
    writeResults(args.out)

if __name__ == "__main__":
    main()
