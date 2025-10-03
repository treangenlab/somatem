#!/usr/bin/env python
import csv
import re
import argparse
import pickle
from Bio import SeqIO
from collections import defaultdict, OrderedDict, namedtuple

def parse_header(header):
    SeqInfo = namedtuple("SeqInfo", ["uniref", "gene_name", "go_terms", "patric_terms"])
    uniref, header_minus_uniref = header.split(' ', 1)
    uniref = uniref.split("_")[1]
    gene_name = re.search(r"(.*?)\sn=\d+", header_minus_uniref).group(1)
    go_groups = [m.group(1).split(';') for m in re.finditer("GO=(.*?)(\s|$)", header_minus_uniref)]
    go_terms = {x for go_group in go_groups for x in go_group }
    patric_terms = [m.group(1) for m in re.finditer(r"PATRIC=(.*?)(\s|$)", header_minus_uniref)]
    return SeqInfo(uniref, gene_name, go_terms, patric_terms)

def get_ancestral_go(go_terms, go_dict):
    ret_terms = set(go_terms)
    queue = list(go_terms)
    while queue:
        gt = queue.pop(0)
        for parent_gt in go_dict[gt]:
            if parent_gt not in ret_terms:
                ret_terms.add(parent_gt)
                queue.append(parent_gt)
    return ret_terms

def main():
    parser = argparse.ArgumentParser(description="Functional report generation")
    parser.add_argument("--fasta", required=True, help="Input fasta")
    parser.add_argument("-u", "--urbtab", required=True, help="Path to btab output")
    parser.add_argument("-g", "--go", required=True, help="Path to go database")
    parser.add_argument("--include-cc", action="store_true", help="Store cellular component go terms as well")
    parser.add_argument("-o", "--out", required=True, help="Path to output tsv file")
    parser.add_argument("-c", "--cutoff", required=True, type=float, help="Bitscore cutoff percentage")
    parser.add_argument("-a", "--ancestral", action="store_true", help="Include all ancestral GO terms")
    parser.add_argument("-r", "--annotation", help="Path to file containing annotation scores")
    args = parser.parse_args()

    fasta_file = args.fasta
    urbtab_file = args.urbtab
    go_file = args.go
    annotation_score_file = args.annotation
    output_file = args.out
    ancestral = args.ancestral
    bitscore_cutoff = (100 - args.cutoff) / 100
    include_cc = args.include_cc

    seq_ids = [seq.id for seq in SeqIO.parse(fasta_file, "fasta")]

    # Load annotation scores
    with open(annotation_score_file, 'rb') as scores:
        score_dict = pickle.load(scores)

    # Load go network
    go_dict = defaultdict(set)
    with open(go_file) as go_file_fd:
        go_reader = csv.DictReader(go_file_fd, delimiter='\t', fieldnames=["Parent", "Child"])
        for line in go_reader:
            go_dict[line["Child"]].add(line["Parent"])

    # Parse blast/diamond output
    Hit = namedtuple("Hit", ["uniref", "evalue", "bitscore", "header"])
    tophits = OrderedDict()
    urbtab_fieldnames = ["query", 'uniref', "2", "3", "4", "5", "6", "7", "8", "9", "evalue", "bitscore", "12", "13", "14", "sequence_header"] # Not sure what numbered columns are, but they are unused
    with open(urbtab_file) as urbtab_file_fd:
        urbtab_reader = csv.DictReader(urbtab_file_fd, delimiter='\t', fieldnames=urbtab_fieldnames)
        for line in urbtab_reader:
            query = line["query"]
            uniref = line["uniref"].split("_")[1]
            
            header = line["sequence_header"]
            if not header.startswith("UniRef"):
                header = line["uniref"] + " " + header
                
            linehit = Hit(
                uniref=uniref,
                evalue=float(line["evalue"]),
                bitscore=float(line["bitscore"]),
                header=header)
            if query in tophits:
                tophit = tophits[query]
                if (linehit.evalue < tophit.evalue) \
                        or (linehit.evalue == tophit.evalue and score_dict[linehit.uniref] > score_dict[tophit.uniref]):
                    tophits[query] = linehit
            else:
                tophits[query] = linehit


    # Go through file again and get hits within cutoff of top score
    tophits_sets = defaultdict(set) # Sets of all hits within cutofff
    tophits_sets.update({query: {hit} for query, hit in tophits.items()})
    with open(urbtab_file) as urbtab_file_fd:
        urbtab_reader = csv.DictReader(urbtab_file_fd, delimiter='\t', fieldnames=urbtab_fieldnames)
        for line in urbtab_reader:
            query = line["query"]
            tophit = tophits[query]
            uniref = line["uniref"].split("_")[1]
            
            header = line["sequence_header"]
            if not header.startswith("UniRef"):
                header = line["uniref"] + " " + header
            
            linehit = Hit(
                uniref=uniref,
                evalue=float(line["evalue"]),
                bitscore=float(line["bitscore"]),
                header=header)
            if (query in tophits) and (bitscore_cutoff * tophits[query].bitscore <= linehit.bitscore):
                tophits_sets[query].add(linehit)
                # Replace tophit because we have a better annotation score
                if score_dict[linehit.uniref] > score_dict[tophit.uniref]:
                    tophits[query] = linehit
                # Replace tophit because we have equal annotation scores but a better bitscore
                elif score_dict[linehit.uniref] == score_dict[tophit.uniref] and linehit.bitscore > tophit.bitscore:
                    tophits[query] = linehit
                # Replace tophit because we have equal annotation scores and bitscores but a better evalue
                elif score_dict[linehit.uniref] == score_dict[tophit.uniref] and linehit.bitscore == tophit.bitscore and linehit.evalue < tophit.evalue:
                    tophits[query] = linehit


    # Write output
    output_fieldnames = ["#query", "evalue", "gene_name", "uniprot", "patric", "go", "pfam", "other_uniprots", "annotation_score"]
    with open(output_file, 'w') as output_fd:
        output_writer = csv.DictWriter(output_fd, fieldnames=output_fieldnames, delimiter='\t', restval='-')
        output_writer.writeheader()
        for query in seq_ids:
            if query not in tophits:
                output_writer.writerow({"#query": query })
                continue
            tophit = tophits[query]
            top_seqinfo = parse_header(tophit.header)
            go_terms = {go for go in top_seqinfo.go_terms if "GO:0005575" not in get_ancestral_go((go), go_dict) or include_cc}
            patric_terms = top_seqinfo.patric_terms
            other_unirefs = set()
            for cutoffhit in tophits_sets[query]:
                hitinfo = parse_header(cutoffhit.header)
                # We are no longer interested in reporting *all* go terms of uniprots within the bitscore
                # go_terms.update(hitinfo.go_terms)
                # patric_terms.update(hitinfo.patric_terms)
                other_unirefs.add(hitinfo.uniref)
            other_unirefs.discard(tophit.uniref)
            output_writer.writerow({
                "#query": query,
                "evalue": tophit.evalue,
                "gene_name": top_seqinfo.gene_name,
                "uniprot": tophit.uniref,
                "patric": ';'.join(patric_terms),
                "go": ';'.join(get_ancestral_go(go_terms, go_dict) if ancestral else go_terms),
                "pfam": "",
                "other_uniprots": ';'.join(other_unirefs),
                "annotation_score": score_dict[tophit.uniref]
            })

if __name__ == "__main__":
    main()

