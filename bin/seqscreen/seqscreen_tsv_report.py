#!/usr/bin/env python3
# coding: utf-8

import argparse
import sys
import csv
from bitarray import bitarray
from Bio import SeqIO
from collections import defaultdict
import pickle
csv.field_size_limit(sys.maxsize)
# from tqdm.notebook import tqdm

def find_superkingdom(taxid, organism, tax_to_parent:dict, merged_tax:dict):
    """Returns the super kingdom (virus, bacteria, eukaryota, other, or not found) given a taxid

    Args:
        taxid (int): taxid 
        tax_to_parent (dict): dictionary of taxid to parent mappings
    """
    if taxid == '-':
        return '-'
    else:
        taxid = int(taxid)
        
    if "phage" in organism.lower():
        return "p"
        
    virus_id = 10239
    phage_ids = [38018, 12333, 77920, 2731619, 10860, 10841, 585893, 2946170, 2842320]
    bacteria_id = 2
    eukaryota_id = 2759
    root_id = 1
    synthetic_id = 32630
    stop_tax = [virus_id, bacteria_id, eukaryota_id, synthetic_id, root_id] + phage_ids
    
    curr_tax = taxid
    while curr_tax not in stop_tax:

        if curr_tax in tax_to_parent:
            curr_tax = tax_to_parent[curr_tax]
        elif curr_tax in merged_tax:
            curr_tax = merged_tax[curr_tax]
        else:
            return "-"
    
    
    if curr_tax == virus_id:
        return 'v'
    elif curr_tax == bacteria_id:
        return 'b'
    elif curr_tax == eukaryota_id:
        return 'e'
    elif curr_tax == synthetic_id:
        return 's'
    elif curr_tax in phage_ids:
        return 'p'
    else:
        return 'o'
    
def flag_sequence(bsat, vfdb, funsoc, superkingdom):
    """Logic for flagging sequences based on output. The logic is that sequences
    are flagged if they have a funsoc hit/vfdb hit, or if they are a virus/eukaryote/synthetic and have a bsat
    hit. If the sequence is deemed a phage, it is ignored and not flagged, even if it has hits. 

    Args:
        bsat (str): "yes" or "no"
        vfdb (str): "yes" or "no"
        funsoc (bool): Has a funsoc been hit for this sequence
        superkingdom (str): one of v, b, e, p, or o for (virus, bacteria, eukaryote, phage, or other)

    Returns:
        flag: 1 if should be flagged, 0 otherwise
    """
    if superkingdom == "p":
        return "0"
    if vfdb == 'yes' or funsoc:
        return "1"
    if superkingdom in ["v", "e", "s"] and bsat=="yes":
        return "1"
    
    return "0"
    

def main():
    parser = argparse.ArgumentParser(description ="tsv report generation")
    parser.add_argument("--functional", required=True, help="Path to functional annotation output")
    parser.add_argument("--taxonomy", required=True, help="Path to taxonomy output")
    parser.add_argument("--taxlookup", required=True, help="Path to NCBI taxonomy lookup file. Contains a mapping of taxid to organism")
    parser.add_argument("--funsocs", required=True, help="Prefix to the funsoc files")
    parser.add_argument("--fasta", required=True, help="Path to input fasta file")
    parser.add_argument("--out", required=True, help="Path to output file to be created")
    parser.add_argument("--mode", required=True, help='Default/Sensitive Mode')
    parser.add_argument("--bsat", nargs="?", help="Path to bsat file")
    parser.add_argument("--vfdb", nargs="?", help="Path to vfdb file")
    parser.add_argument("--parenttax", required=True, help="Path to the pickle file containing dictionary of tax to parent id mappings")
    parser.add_argument("--mergedtax", required=True, help="Path to merged_taxa.pck file")
    args = parser.parse_args()

    tax_lookup_file = args.taxlookup
    tax_results_file = args.taxonomy
    func_results_file = args.functional
    funsocs_file = args.funsocs
    fasta_file = args.fasta
    output_file = args.out
    mode = True if args.mode == "fast" else False
    bsat_file = args.bsat
    vfdb_file = args.vfdb
    parenttax_file = args.parenttax
    mergedtax_file = args.mergedtax

    funsoc_headers = [
        "disable_organ", "cytotoxicity", "degrade_ecm", "induce_inflammation", 
        "bacterial_counter_signaling", "viral_counter_signaling", "resist_complement", 
        "counter_immunoglobulin", "plant_rna_silencing", "resist_oxidative", "suppress_detection", 
        "avirulence_plant", "host_gtpase", "host_transcription", "host_translation", 
        "host_ubiquitin", "host_xenophagy", "nonviral_invasion", 
        "viral_invasion", "viral_movement", "virulence_activity", "host_cell_cycle", 
        "host_cell_death", "host_cytoskeleton", "secreted_effector", "antibiotic_resistance", 
        "develop_in_host", "nonviral_adhesion", "secretion", "toxin_synthase", "viral_adhesion", 
        "virulence_regulator"
    ]
    if mode: ## fast mode output
        if bsat_file and vfdb_file: 
            output_header = ["query", "taxid", "centrifuge_multi_tax","diamond_multi_tax","go", "multi_taxids_confidence", "go_id_confidence", "bsat_hit", "vfdb_hit"] \
            + funsoc_headers \
            + ["vbeo", "size", "organism", "gene_name", "uniprot", "uniprot evalue", "flag"]
        else: ## ont mode
            output_header = ["query", "taxid", "centrifuge_multi_tax","diamond_multi_tax","go", "multi_taxids_confidence", "go_id_confidence"] + funsoc_headers \
            + ["vbeo", "size", "organism", "gene_name", "uniprot", "uniprot evalue"]
            
    else: ## sensitive mode
                
        if bsat_file and vfdb_file:
            output_header = ["query", "taxid","go", "multi_taxids_confidence", "go_id_confidence","bsat_hit", "vfdb_hit"] + funsoc_headers \
            + ["vbeo", "size", "organism", "gene_name", "uniprot", "uniprot evalue", "flag"]
        else: ## ont mode
            output_header = ["query", "taxid","go", "multi_taxids_confidence", "go_id_confidence"] + funsoc_headers \
            + ["vbeo", "size", "organism", "gene_name", "uniprot", "uniprot evalue"]
            

    taxid_to_taxname = defaultdict(lambda: '-')
    seqid_to_seqrecord = {}
    seqid_to_tax = defaultdict(list) # seq_id : [{source : x, tax_id : y, conf : z} , ...]
    seqid_to_func = defaultdict(lambda: {"go": [], "uniprot" : "-", "evalue": '-', "gene_name": '-'}) # seq_id : {go : [{goid : y, conf : z} ... ], uniprot : x, "evalue"}
    uniprot_to_funsocs = dict() 


    # Read in taxonomy lookup file
    with open(tax_lookup_file) as stream:
        reader = csv.DictReader(stream, delimiter='\t', fieldnames=["tax_id", "tax_name"])
        for line in reader:
            taxid_to_taxname[line["tax_id"]] = line["tax_name"]
        
    # Read in taxonomy to parent map and merged taxa
    with open(parenttax_file, 'rb') as handle:
        tax_to_parent = pickle.load(handle)
    with open(mergedtax_file, 'rb') as handle:
        merged_tax = pickle.load(handle)        

    # Read in funsocs text file to get header info
    with open(funsocs_file + ".tsv") as stream:
        reader = csv.reader(stream, delimiter='\t')
        funsocs_fieldnames = next(reader)
        funsocs_fieldnames = [f for f in funsocs_fieldnames if f in funsoc_headers]

    with open(funsocs_file + ".pck", 'rb') as funsocs_pck_fd:
        uniprot_to_funsocs = pickle.load(funsocs_pck_fd)

    with open(funsocs_file + "_commensal_list.pck", 'rb') as funsocs_pck_fd:
        nofunsoc_ids = pickle.load(funsocs_pck_fd)

    # Read in fasta file
    for record in SeqIO.parse(fasta_file, "fasta"):
        seqid_to_seqrecord[record.id] = record

    # Parse taxonomy results
    with open(tax_results_file) as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        for line in reader:
            sources = line["source"].split(',')
            taxids = line["taxid"].split(',')
            confidences = line["confidence"].split(',')
            if mode:
                centrifuge_multi_tax = line["centrifuge_multi_tax"]
                diamond_multi_tax = line["diamond_multi_tax"]
                seqid_to_tax[line["#query"]] = [[{"source": source, "taxid": taxid, "conf": float(confidence)} 
                for source, taxid, confidence in zip(sources, taxids, confidences)], centrifuge_multi_tax,diamond_multi_tax]
                assert(len(seqid_to_tax[line["#query"]]) == 3)
            else:
                seqid_to_tax[line["#query"]] = [[{"source": source, "taxid": taxid, "conf": float(confidence)} 
                for source, taxid, confidence in zip(sources, taxids, confidences)]]
                assert(len(seqid_to_tax[line["#query"]]) == 1)
               
       
    # Parse functional results
    with open(func_results_file) as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        for line in reader:
            goids = line["go"].split(';')
            seqid_to_func[line["#query"]] = {
                "go" : [{"goid": goid, "conf": 1.0} for goid in line["go"].split(';')if goid.strip() != ''], 
                "uniprot": line["uniprot"], 
                "evalue": line["evalue"], 
                "gene_name": line["gene_name"]
            }

   
   
    # Parse bsat
    seqid_bsat_hits = set()
    if bsat_file:
        with open(bsat_file) as bsat_f:
            reader = csv.DictReader(bsat_f,delimiter="\t")
            for line in reader:
                seqid_bsat_hits.add(line["query"])
                
    # Parse vfdb
    seqid_vfdb_hits = set()
    if vfdb_file:
        with open(vfdb_file) as vfdb_f:
            reader = csv.DictReader(vfdb_f,delimiter="\t")
            for line in reader:
                seqid_vfdb_hits.add(line["query"])


    with open(output_file, 'w') as stream:
        writer = csv.DictWriter(stream, delimiter='\t', fieldnames=output_header)
        writer.writeheader()
        for seqid, record in seqid_to_seqrecord.items():
            seq_dict = {"query": seqid, "size": len(record.seq)}
            if not seqid_to_tax[seqid]:
                tax_info = []
            else:    
                tax_info = seqid_to_tax[seqid][0]
            taxid = tax_info[0]["taxid"] if len(tax_info) > 0 else '-'
            taxconf = tax_info[0]["conf"] if len(tax_info) > 0 else '-'
            if mode:
                if len(tax_info) > 0:
                    centrifuge_multi_tax_res = seqid_to_tax[seqid][1] if not seqid_to_tax[seqid][1] == "" else "-"
                    diamond_multi_tax_res = seqid_to_tax[seqid][2] if not seqid_to_tax[seqid][2] == "" else "-"
                else:
                    centrifuge_multi_tax_res = "-"
                    diamond_multi_tax_res = "-"
                seq_dict["centrifuge_multi_tax"] = centrifuge_multi_tax_res
                seq_dict["diamond_multi_tax"] = diamond_multi_tax_res
            if bsat_file:
                seq_dict["bsat_hit"] = "yes" if seqid in seqid_bsat_hits else "no"
            if vfdb_file:
                seq_dict["vfdb_hit"] = "yes" if seqid in seqid_vfdb_hits else "no"
            seq_dict["taxid"] = taxid
            
            seq_dict["multi_taxids_confidence"] = ",".join(
                entry["taxid"] + ':' + str(round(entry["conf"], 3)) 
                    for entry in tax_info) if len(tax_info) > 0 else '-'
            seq_dict["organism"] = taxid_to_taxname[taxid]
            seq_dict["gene_name"] = seqid_to_func[seqid]["gene_name"]
            uniprot = seqid_to_func[seqid]["uniprot"]
            seq_dict["uniprot"] = uniprot
            seq_dict["uniprot evalue"] = seqid_to_func[seqid]["evalue"]
            seq_dict["go"] = ';'.join(go_hit["goid"] for go_hit in seqid_to_func[seqid]["go"])
            seq_dict["go_id_confidence"] = ';'.join(
                go_hit["goid"] + "[" + str(go_hit["conf"]) + "]"  
                    for go_hit in seqid_to_func[seqid]["go"])
            
            seq_dict["vbeo"] = find_superkingdom(taxid, seq_dict["organism"], tax_to_parent, merged_tax)

            
            funsoc_flagged = False ##for sequence flagging logic
            for idx in range(len(funsoc_headers)):
                if uniprot in uniprot_to_funsocs:
                    if uniprot_to_funsocs[uniprot][idx]:
                        funsoc_name = funsocs_fieldnames[idx]
                        seq_dict[funsoc_name] = '1'
                        if funsoc_name != 'antibiotic_resistance' and funsoc_name != 'secretion': ##ignore antibiotic resistance and secretion funsocs
                            funsoc_flagged = True
                    else: 
                        seq_dict[funsocs_fieldnames[idx]] = '0'
                elif uniprot in nofunsoc_ids:
                    seq_dict[funsocs_fieldnames[idx]] = '0'
                else:
                    seq_dict[funsocs_fieldnames[idx]] = '-'
            
            if bsat_file and vfdb_file:
                seq_dict["flag"] = flag_sequence(seq_dict["bsat_hit"], seq_dict["vfdb_hit"], funsoc_flagged , seq_dict["vbeo"])
            writer.writerow(seq_dict)


if __name__ == "__main__":
    main()
