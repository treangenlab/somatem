#!/usr/bin/env python3
import os
import sys
import argparse
import glob
import numpy as np
import pickle
import pandas as pd


def parseArguments(_args):
    parser = argparse.ArgumentParser(prog="format", usage="%(prog)s [options]",
                                   description=f"SeqScreen-Nano: tool for pathogen identification from ONT data. Formats output file.",
                                   epilog="Author: {__author__}\n")
    parser.add_argument("--report", type=str, required=True, help="SeqScreen report")
    parser.add_argument("--format", type=int, required=True, help="Format type: [1] Original, [2] Hits only, [3] FunSoC only, [4] Gene-Centric,\
    [5] Gene-Centric FunSoC Only [Default: 1]", choices = [1,2,3,4,5])
    parser.add_argument("--mode", type=str, required=True, help="Mode")
    parser.add_argument("--databases", type=str, required=True, help="Databases path")
    parser.add_argument("--taxonomy_confidence_threshold", type=float, default=0.0, help="taxonomy confidence threshold (average for multi-tax)")
    parser.add_argument("--filter-taxon", type=str, required=True, help="Filter out these taxons [Default: None]")
    parser.add_argument("--keep-taxon", type=str, required=True, help="Only keep these taxons [Default: All]")
    
    parseArgs = parser.parse_args(args=_args)
    return parseArgs


def calculate_conf(conf):
    res = []
    conf_list = conf.split(",")
    for c in conf_list:
        if not "N/A" in c:
            res.append(float(c.split(":")[1]))
    
    if not res:
        return np.inf
    return round(sum(res)/len(res), 2)


def create_pathgo(tsv_file, _out, _databases):
    pgodb = os.path.join(_databases, "id2pgo.pickle")
    with open(pgodb, "rb") as f:
        id2pgo = pickle.load(f)
        
    df = pd.read_csv(tsv_file, header=0, sep="\t")
    df["PathGO"] = [id2pgo.get(k.strip(),"-") for k in df["uniprot"].tolist()]
    df.to_csv(_out, sep="\t", header=True, index=False)


def write_F1(_in, _out, FTL, KTL, ftl, ktl, conf, conf_index):
    with open(_in,"r") as f, open(_out, "w+") as wf:                   
        header = f.readline()
        wf.write(header)
        for line in f:
            WRITE = True
            tokens = line.strip().split("\t")
            if tokens[1] == "-":
                continue
                
            if FTL:
                if tokens[1] in ftl:
                    WRITE = False
            if KTL:
                if tokens[1] not in ktl:
                    WRITE = False
            
            if tokens[conf_index] != "-" and calculate_conf(tokens[conf_index]) < conf:
                WRITE = False
            
            if WRITE:
                wf.write(line) 


def write_F2(_in, _out, FTL, KTL, ftl, ktl, conf, conf_index):
    with open(_in,"r") as f, open(_out, "w+") as wf:
        header = f.readline()
        wf.write(header)
        for line in f:
            WRITE = True
            tokens = line.strip().split("\t")
            if tokens[1] == "-":
                continue
            if tokens[-1] != "-" and tokens[-2] != "-":
                if FTL:
                    if tokens[1] in ftl:
                        WRITE = False
                if KTL:
                    if tokens[1] not in ktl:
                        WRITE = False
                
                if tokens[conf_index] != "-" and calculate_conf(tokens[conf_index]) < conf:
                    WRITE = False
                
                if WRITE:
                    wf.write(line)
                    
def write_F3(_in, _out, FTL, KTL, ftl, ktl, idx, conf, conf_index):
    with open(_in, "r") as f, open(_out, "w+") as wf:
        header = f.readline()
        wf.write(header)
        for line in f:
            WRITE = True
            tokens = line.strip().split("\t")
            if tokens[1] == "-":
                continue
            funsocs = tokens[idx:idx+32]
            if funsocs.count("1"):
                if FTL:
                    if tokens[1] in ftl:
                        WRITE = False
                if KTL:
                    if tokens[1] not in ktl:
                        WRITE = False
                        
                if tokens[conf_index] != "-" and calculate_conf(tokens[conf_index]) < conf:
                    WRITE = False
                
                if WRITE:
                    wf.write(line)
                

def write_F4_5(_in, _out, FTL, KTL, ftl, ktl, idx, FUNSOCS, _mode, conf, conf_index):
    funsoc_headers = ["disable_organ", "cytotoxicity", "degrade_ecm", "induce_inflammation", \
                      "bacterial_counter_signaling", "viral_counter_signaling", "resist_complement", \
                      "counter_immunoglobulin", "plant_rna_silencing", "resist_oxidative", "suppress_detection", \
                      "avirulence_plant", "host_gtpase", "host_transcription", "host_translation", \
                      "host_ubiquitin", "host_xenophagy", "nonviral_invasion", \
                      "viral_invasion", "viral_movement", "virulence_activity", "host_cell_cycle", \
                      "host_cell_death", "host_cytoskeleton", "secreted_effector", "antibiotic_resistance", \
                      "develop_in_host", "nonviral_adhesion", "secretion", "toxin_synthase", "viral_adhesion", \
                      "virulence_regulator"]
    
    
    header = ["UniProt", "Taxid", "Gene Name", "UniProt evalues", "GO", "Number of Reads"] + funsoc_headers
    go_idx = 4
    if _mode == "ont":
        header[5] = "Number of ORFs"
        go_idx = 2
    if _mode == "sensitive":
        go_idx = 2 
    
    gene2info = {}
    
    with open(_in, "r") as f:
        h = f.readline()
        for line in f:
            tokens = line.strip().split("\t")
            if tokens[1] == "-":
                continue
            if tokens[-2] not in gene2info:
                if tokens[-1] != "-" and tokens[-2] != "-":
                    WRITE = True

                    if FTL:
                        if tokens[1] in ftl:
                            WRITE = False
                    if KTL:
                        if tokens[1] not in ktl:
                            WRITE = False

                    if FUNSOCS:
                        funsocs = tokens[idx:idx+32]
                        if not funsocs.count("1"):
                            WRITE = False
                            
                    if tokens[conf_index] != "-" and calculate_conf(tokens[conf_index]) < conf:
                        WRITE = False


                    if WRITE:
                        gene2info[tokens[-2]] = [tokens[1], tokens[-3], [tokens[-1]], tokens[go_idx], 1] + tokens[idx:idx+32]
            else:
                gene2info[tokens[-2]][2].append(tokens[-1])
                gene2info[tokens[-2]][4] += 1
                        
                        
    with open(_out, "w+") as wf:
        wf.write("\t".join(header)+"\n")
        for upid in gene2info:
            gene2info[upid][2] = ",".join(gene2info[upid][2])
            gene2info[upid][4] = str(gene2info[upid][4])
            wf.write("\t".join([upid]+gene2info[upid])+"\n")
            
            

            
            
def main(): 
    not_set_list = [[""],["''"],"",[],["true"]]
    parseArgs = parseArguments(sys.argv[1:])
    _report = parseArgs.report
    _format = parseArgs.format
    _mode = parseArgs.mode
    _databases = parseArgs.databases
    _conf = parseArgs.taxonomy_confidence_threshold

    FTL = False
    KTL = False
    ftl = parseArgs.filter_taxon.split(",")
    ktl = parseArgs.keep_taxon.split(",") 

    if ftl not in not_set_list:
        ftl = set(ftl)
        FTL = True
    if  ktl not in not_set_list:
        ktl = set(ktl)
        KTL = True
        
    funsoc_index = 0
    conf_index = 3
    if _mode == "fast":
        funsoc_index = 7
        conf_index = 5
    elif _mode == "sensitive":
        funsoc_index = 6
    else:
        funsoc_index = 5
    
    tsv_file = sorted(glob.glob(os.path.join(_report, "*.tsv")))[0]
    

    output_prefix = os.path.basename(tsv_file).split(".tsv")[0]
    _report_mod_path = os.path.join(_report, output_prefix + "_f"+str(_format)+".tsv")
    _report_pathgo_path = os.path.join(_report, output_prefix + "_pathgo.tsv")
    
    create_pathgo(tsv_file, _report_pathgo_path, _databases)
 
    if _format == 1:
        if FTL or KTL or _conf:
            write_F1(tsv_file, _report_mod_path, FTL, KTL, ftl, ktl, _conf, conf_index)
    elif _format == 2:
        write_F2(tsv_file, _report_mod_path, FTL, KTL, ftl, ktl, _conf, conf_index)
    elif _format == 3:
        write_F3(tsv_file, _report_mod_path, FTL, KTL, ftl, ktl, funsoc_index, _conf, conf_index)
    elif _format == 4:
        write_F4_5(tsv_file, _report_mod_path, FTL, KTL, ftl, ktl, funsoc_index, False, _mode, _conf, conf_index)
    elif _format == 5:
        write_F4_5(tsv_file, _report_mod_path, FTL, KTL, ftl, ktl, funsoc_index, True, _mode, _conf, conf_index)

if __name__ == '__main__':    
    main()
