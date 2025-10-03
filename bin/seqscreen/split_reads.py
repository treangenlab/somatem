#!/usr/bin/env python3

'''
Copyright 2021, Advait Balaji, All rights reserved.
Free for research use. For commercial use, contact
Dr. Todd Treangen (treangen@rice.edu)
'''

import argparse
import os
import sys
from collections import defaultdict, OrderedDict

from Bio import SeqIO
from intervaltree import IntervalTree
from copy import deepcopy


def parseArguments(_args):
    parser = argparse.ArgumentParser(prog="split_reads", usage="%(prog)s [options]",
                                     description=f"SeqScreen-Nano: tool for pathogen identification from ONT data. Splits reads based on blastx.",
                                     epilog="Author: {__author__}\n")
    parser.add_argument("--fasta", type=str, required=True, help="Query file")
    parser.add_argument("--blastx", type=str, required=True, help="Blastx btab file")
    parser.add_argument("--unmapped", action="store_true", help="Running unmapped reads")
    parseArgs = parser.parse_args(args=_args)
    return parseArgs


def main():
    parseArgs = parseArguments(sys.argv[1:])

    if parseArgs.unmapped:
        first_fasta = parseArgs.fasta.split(".unmapped.fasta")[0]+".split.first.fasta"
        splitReads(parseArgs.fasta, parseArgs.blastx, parseArgs.unmapped)
        mergeInfo(os.path.join(os.path.dirname(parseArgs.blastx), os.path.basename(parseArgs.fasta)))
    else:
        splitReads(parseArgs.fasta, parseArgs.blastx, parseArgs.unmapped)


def splitReads(fasta, btab, unmapped):

    first_mapped_dict = defaultdict(int)
    first_mapped_dict_copy = defaultdict(int)
    first_rc_dict = defaultdict(list)
    first_rc_num_dict = defaultdict(int)
    unmapped_count_dict = defaultdict(int)
    
    sbtab = os.path.join(os.path.dirname(btab), os.path.basename(fasta) + ".ur100.split.first.btab")
    split_fasta = os.path.join(os.path.dirname(btab), os.path.basename(fasta) + ".split.first.fasta")
    unmapped_fasta = os.path.join(os.path.dirname(btab), os.path.basename(fasta) + ".unmapped.fasta")
    read_coordinate_file = os.path.join(os.path.dirname(btab), os.path.basename(fasta) + ".rc.first.tsv")

    if unmapped:
        _fasta = os.path.basename(fasta).split(".unmapped.fasta")[0]
        first_fasta = os.path.join(os.path.dirname(btab),_fasta+".split.first.fasta")
        first_rc = os.path.join(os.path.dirname(btab),_fasta+".rc.first.tsv")
        sbtab = os.path.join(os.path.dirname(btab), _fasta + ".ur100.split.second.btab")
        split_fasta = os.path.join(os.path.dirname(btab), _fasta + ".split.second.fasta")
        read_coordinate_file = os.path.join(os.path.dirname(btab), _fasta + ".rc.second.tsv")
        
        for record in SeqIO.parse(first_fasta,"fasta"):
            if not "unmapped" in record.id:
                _id,cnt = "_".join(record.id.split("_")[:-1]),record.id.split("_")[-1]
                first_mapped_dict[_id] = int(cnt)
                
        first_mapped_dict_copy = deepcopy(first_mapped_dict) #create deepcopy
                
        with open(first_rc,"r") as f:
            f.readline() #read header
            _count = 0
            for line in f:
                tokens = line.strip().split("\t")
                if tokens[1] == "Region unmapped":
                    first_rc_dict[tokens[0]].append(int(tokens[-2]))
        

        for k in first_rc_dict:
            for i,v in enumerate(first_rc_dict[k]):
                first_rc_num_dict[k+"_unmapped_"+str(i)] = v

        del first_rc_dict


    with open(btab, "r") as f, open(sbtab, "w+") as wf:
        seen_ids = OrderedDict()  # Set of ids
        seen_trip = defaultdict(OrderedDict)  # Dict of id and qstart, qend OrderedDict()
        count = 0
        for line in f:
            tokens = line.split("\t")
            r_id, qstart, qend = tokens[0], int(tokens[6]), int(tokens[7])
            write_id = r_id
            if unmapped:
                write_id = r_id.split("_unmapped")[0]
            if not r_id in seen_ids:
                count = first_mapped_dict.get(write_id,-1)+1 #adjust count
                seen_ids[r_id] = None
                seen_trip[r_id][(qstart, qend)] = None
                wf.write(write_id + "_" + str(count) + "\t" + "\t".join(tokens[1:]))
                if unmapped:
                    first_mapped_dict[write_id] = count
            else:
                if (qstart, qend) in seen_trip[r_id]:
                    wf.write(write_id + "_" + str(count) + "\t" + "\t".join(tokens[1:]))
                else:
                    count += 1
                    wf.write(write_id + "_" + str(count) + "\t" + "\t".join(tokens[1:]))
                    seen_trip[r_id][(qstart, qend)] = None
                    if unmapped:
                        first_mapped_dict[write_id] = count
                

    with open(split_fasta, "w+") as wf, open(unmapped_fasta, "w+") as uf, open(read_coordinate_file, "w+") as rcf:
        rcf.write("RID\tMapStatus\tQstart\tQend\n")
        for record in SeqIO.parse(fasta, "fasta"):
            # Retrieve the ID and sequence for each record
            _id = record.id
            _seq = record.seq
            write_id = _id

            if unmapped:
                write_id = _id.split("_unmapped")[0]

            if not _id in seen_ids:  
                if not unmapped:  # Write to unmapped file
                    wf.write(">" + write_id + "_unmapped\n")
                    wf.write(str(_seq) + "\n")
                    uf.write(">" + write_id + "_unmapped\n")
                    uf.write(str(_seq) + "\n")
                    rcf.write(write_id + "\t" + "Read unmapped" + "\t0\t" + str(len(_seq) - 1) + "\n")
                else:
                    if not "_unmapped_" in _id: #The whole read unmapped still unmapped
                        wf.write(">" + write_id + "_unmapped\n")
                        wf.write(str(_seq) + "\n")
                        rcf.write(write_id + "\t" + "Read unmapped" + "\t0\t" + str(len(_seq) - 1) + "\n")
                    else: #Region unmapped still unmapped
                        offset = first_rc_num_dict[_id]
                        wf.write(">" + write_id + "_unmapped_"+str(unmapped_count_dict[write_id])+"\n")
                        wf.write(str(_seq) + "\n")
                        unmapped_count_dict[write_id]+=1
                        rcf.write(write_id + "\t" + "Region unmapped" + "\t" + str(offset) + "\t" + str(offset + len(_seq) - 1) + "\n")


            else:
                # Create an interval tree for each mapped record
                t = IntervalTree()
                t[1:len(str(_seq))] = "m"
                count = first_mapped_dict_copy.get(write_id,-1)+1 #adjust count using deepcopy to restart count
                offset = 0

                if "_unmapped_" in _id and unmapped:
                    offset = first_rc_num_dict[_id]

                for qpair in seen_trip[_id]:
                    wf.write(">" + write_id + "_" + str(count) + "\n")
                    _lbound, _ubound = min(qpair), max(qpair)
                    t.chop(_lbound, _ubound)  # chop the interval seen
                    wf.write(str(_seq[_lbound - 1:_ubound]) + "\n")
                    rcf.write(write_id + "\t" + "Region mapped [ORF]" + "\t" + str(offset+_lbound - 1) + 
                              "\t" + str(offset+_ubound - 1) + "\n")
                    count += 1

                first_mapped_dict_copy[write_id] = count-1 #For debug

                count = unmapped_count_dict[write_id] if unmapped else 0 
                for _int in t:  # iterate over remaining intervals and check if they can be output
                    _lbound, _ubound = min((_int[0], _int[1])), max((_int[0], _int[1]))
                    if _ubound - _lbound >= 90:  # at least 90bp (30aa)
                        if _lbound == 1:
                            _lbound -= 1
                        if _ubound == len(str(_seq)):
                            _ubound += 1
                        wf.write(">" + write_id + "_unmapped_" + str(count) + "\n")
                        wf.write(str(_seq[_lbound:_ubound - 1]) + "\n")
                        uf.write(">" + write_id + "_unmapped_" + str(count) + "\n") #write even in unmapped
                        uf.write(str(_seq[_lbound:_ubound - 1]) + "\n")
                        rcf.write(write_id + "\t" + "Region unmapped" + "\t" + str(offset+_lbound) + 
                                  "\t" + str(offset+_ubound - 2) + "\n")
                        count += 1

                unmapped_count_dict[write_id] = count #For debug
                
                del seen_trip[_id]
                del seen_ids[_id]


def mergeInfo(base):
    base = base.split(".unmapped.fasta")[0]

    sbtab = base + ".ur100.split.btab"
    split_file = base + ".split.fasta"
    read_coordinate_file = base + ".rc.tsv"

    with open(base + ".ur100.split.first.btab", "r") as inpf1, open(base + ".ur100.split.second.btab",
                                                                    "r") as inpf2, open(sbtab, "w+") as wf:
        for line in inpf1:
            wf.write(line)
        for line in inpf2:
            wf.write(line)

    with open(base + ".split.first.fasta", "r") as inpf1, open(base + ".split.second.fasta", "r") as inpf2, open(
            split_file, "w+") as wf:
        flag = False
        for line in inpf1:
            if ">" in line:
                if not "unmapped" in line.strip():
                    wf.write(line)
                    flag = True
            else:
                if flag:
                    wf.write(line)
                    flag = False

        for line in inpf2:
            wf.write(line)

    with open(base + ".rc.first.tsv", "r") as inpf1, open(base + ".rc.second.tsv", "r") as inpf2, open(
            read_coordinate_file, "w+") as wf:
        for line in inpf1:
            tokens = line.split("\t")
            if "unmapped" not in tokens[1]:
                wf.write(line)
        
        h2 = inpf2.readline() #skip header
        for line in inpf2:
            wf.write(line)


if __name__ == '__main__':
    main()
