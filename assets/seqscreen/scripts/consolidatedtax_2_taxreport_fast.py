#!/usr/bin/env python3
'''
AUTHOR: Advait Balaji
DATE: 07/07/2020

DESC: This script takes in the centrifuge and diamond outputs 
and produces a tab separated file of form as required by seqscreen

#query taxid source confidence

consolidatedtax_2_taxreport_fast.py --diamond=/Path/to/diamondoutput.txt --centrifuge=/Path/to/centrifugeoutput
                                     --out=/Path/to/output.txt  --cutoff=1 [--help] 

'''

'''
imports
'''
import argparse
import sys
import csv
import numpy as np



top_hits = {} #top hit in diamond
top_hits_centrifuge = {} #top hit in centrifuge
centrifuge_results_dict = {} #store centrifuge results
diamond_results_dict = {}#try storing dimond results separately

def usage():
    print(f"consolidatedtax_2_taxreport_fast.py \
        --diamond=/Path/to/diamondoutput.txt \
        --centrifuge=/Path/to/centrifugeoutput \
        --out=/Path/to/output.txt \
        --cutoff=1 [--help] ")

'''
processes the centrifuge output file
    @params: centrifuge_f -> centrifuge output file (string)
'''

def processCentrifuge(centrifuge_f, taxlimit):
    with open(centrifuge_f,'r') as inpf:
        inpf_reader = csv.DictReader(inpf,delimiter="\t")
        for row in inpf_reader:
            query = base_query(row["readID"]) #read classified
            taxid = row["taxID"] #tax id
            score = row["score"] #centrifuge score
            seqid = row["seqID"] #seqid for hit in DB
            hit_length = row["hitLength"]
            query_length = row["queryLength"]
            if not query:
                continue
            if query in centrifuge_results_dict and len(centrifuge_results_dict[query][0]) >= taxlimit:
                continue
            if taxid != "0":
                centrifuge_conf = int(hit_length)/int(query_length) #get confidence score
                if score == top_hits_centrifuge[query]: #check if this is the top score
                    if query in centrifuge_results_dict: # already present 
                        if not taxid in centrifuge_results_dict[query][0]: #don't take centrifuge duplicates
                            centrifuge_results_dict[query][0].add(taxid)
                            centrifuge_results_dict[query][1].append(score)
                            centrifuge_results_dict[query][2].append(str(centrifuge_conf))
                            centrifuge_results_dict[query][3].append("centrifuge")
                    else:
                        centrifuge_results_dict[query] = (set([taxid]),[score],[str(centrifuge_conf)],["centrifuge"])


'''
processes the diamond output file
    @params: diamond_f -> diamond output file (string)
    @params: cutoff -> cutoff value to consider (int)
'''
def processDiamond(diamond_f, cutoff, taxlimit):
    with open(diamond_f,'r') as inpf:
        for line in inpf:
            tokens = line.split('\t')
            query = base_query(tokens[0])
            max_bit_score = float(tokens[7])
            bitscore = float(tokens[11])
            taxid = tokens[15].split("ID=")[1].split(" ")[0]
            bitscore_cutoff = cutoff * top_hits[query]
            if query in diamond_results_dict and len(diamond_results_dict[query][0]) >= taxlimit:
                continue
            if not query:
                continue
            if bitscore > bitscore_cutoff: #we are good to take this if not already predicted by centrifuge
                if query in diamond_results_dict:    
                    if not taxid in diamond_results_dict[query][0]:
                        diamond_results_dict[query][0].add(taxid)
                        diamond_results_dict[query][1].append(bitscore)
                        diamond_results_dict[query][2].append(str(float(tokens[2])/100))
                        diamond_results_dict[query][3].append("diamond")
                else:
                    diamond_results_dict[query] = (set([taxid]),[bitscore],[str(float(tokens[2])/100)],["diamond"])

'''
finds the top hits in diamond output file and stores  in top_hit
    @params: diamond_f -> diamond output file (string)
'''

def findTopHitsDiamond(diamond_f):
    with open(diamond_f,'r') as inpf:
        for line in inpf:
            tokens = line.split('\t')
            query = base_query(tokens[0])
            bitscore = float(tokens[11])
            if bitscore < 0:
                print(f"Bitscore was found to be less that zero, might be the wrong column")
                exit(1)
            if query in top_hits:
                if top_hits[query] < bitscore:
                    top_hits[query] = bitscore
            else:
                top_hits[query] = bitscore
'''
finds top hits in centrifuge and stores in top_hit_centrifuge
    @params: centrifuge_f -> centrifuge output file (string)
'''

def findTopHitsCentrifuge(centrifuge_f):
    with open(centrifuge_f,"r") as inpf:
        inpf_reader = csv.DictReader(inpf,delimiter="\t")
        for row in inpf_reader:
            query = base_query(row["readID"]) #read classified
            score = row["score"] #centrifuge score
            if query in top_hits_centrifuge:
                if top_hits_centrifuge[query] < score:
                    top_hits_centrifuge[query] = score
            else:
                top_hits_centrifuge[query] = score
        


'''
write the output file that can be read by the report generation workflow
    @params: out_f -> output file (string)
'''
def writeResults(out_f):
    with open(out_f,'w') as outf:
        outf.write("#query\ttaxid\tsource\tconfidence\tcentrifuge_top_priority\tcombined_taxid\tcentrifuge_multi_tax\tdiamond_multi_tax\n")
        keys = set(list(centrifuge_results_dict.keys()) + list(diamond_results_dict.keys()))
        for k in keys:
            centrifuge_multi_priority,source,confidence,centrifuge_top_priority = "","","",""
            if k in centrifuge_results_dict:
                centrifuge_multi_priority = centrifuge_results_dict[k][0]
                source = centrifuge_results_dict[k][3]
                confidence = centrifuge_results_dict[k][2]
                centrifuge_top_priority = [sorted(list(centrifuge_results_dict[k][0]))[0]]
            else:
                centrifuge_multi_priority = diamond_results_dict[k][0]
                source = diamond_results_dict[k][3]
                confidence = diamond_results_dict[k][2]
                centrifuge_top_priority = sorted(diamond_results_dict[k][0])
            centrifuge_multitax_res = ",".join(centrifuge_results_dict.get(k,[[]])[0]) #either string of taxids or ""
            diamond_multitax_res = ",".join(diamond_results_dict.get(k,[[]])[0]) #either string of taxids or ""
            outf.write(k+"\t"+",".join(centrifuge_multi_priority)+"\t"+",".join(source)+"\t"+",".join(confidence)+"\t"+",".join(centrifuge_top_priority)+"\t"+",".join(set(list(centrifuge_results_dict.get(k,[[]])[0])+list(diamond_results_dict.get(k,[[]])[0])))+"\t"+centrifuge_multitax_res+"\t"+diamond_multitax_res+"\n")
'''
fixes readID in diamond and remove __unambig__ token
    @params: readID from diamond file
'''
def base_query(id):
    return id.split("_unambig_")[0]


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-d","--diamond",type=str,required=True,help="Input: /Path/to/Diamondfile.txt")
    parser.add_argument("-c","--centrifuge",type=str,required=True,help="Input: /Path/to/Centrifugefile.txt")
    parser.add_argument("-o","--out",type=str,required=True,help="Output: /Path/to/Outputfile.txt")
    parser.add_argument("-t","--cutoff",type=int,default=1,help="Cutoff for top n bit scores from max bit score (for ex: 1,5,10...)")
    parser.add_argument("--taxlimit", type=int, help="Maximum number of taxids to report for a single input", default=25)

    args = parser.parse_args() #save args

    cutoff = (100-args.cutoff)/100

    findTopHitsCentrifuge(args.centrifuge)
    findTopHitsDiamond(args.diamond)
    processCentrifuge(args.centrifuge, args.taxlimit)
    processDiamond(args.diamond, cutoff, args.taxlimit)
    #processCentrifuge(args.centrifuge)
    writeResults(args.out)

if __name__  == "__main__":
    main()

