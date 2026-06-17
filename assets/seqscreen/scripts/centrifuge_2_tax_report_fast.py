#!/usr/bin/env python3
'''
AUTHOR: Advait Balaji
DATE: 07/07/2020

DESC: This script takes in the centrifuge results tsv and outputs 
a tab separated file of form as required by seqscreen

#query taxid source

centrifuge_2_tax_report_fast.py --blastx=/Path/to/ --out=/Path/to/output.txt 
                     [--help] 

'''


# add the imports
import argparse
import sys
import csv

'''
print usage 
'''
def usage():
    print(f"centrifuge_2_tax_report_fast.py \
        --blastx=/Path/to/blastx.centrifuge \
        --out=/Path/to/output.txt \
        [--help] ")

'''
processes the centrifuge output file
    @params: centrifuge -> centrifuge output file (string)
    @params: out -> taxonomy classification output file (string)
'''

def processCentrifuge(centrifuge,out):
    results_dict = {}
    with open(centrifuge) as assignments_file:
        assignments_reader = csv.DictReader(assignments_file, delimiter="\t")
        for row in assignments_reader:
            sample_id = row["readID"]
            tax_id = row["taxID"]
            score = row["score"]
            if sample_id in results_dict:
                if results_dict[sample_id][1] < score:
                    results_dict[sample_id][1] = score
            else:
                results_dict[sample_id] = (tax_id,score)

    with open(out,"w+") as outf:
        outf.write("#readID\ttaxID\tscore\n")
        for k,v in results_dict.items():
            outf.write(k+"\t"+v[0]+"\t"+str(v[1])+"\n")

'''
main function
    parses args and calls processCentrifuge
    needs exactly 2 cmd line params to run
'''

def main():

    if (len(sys.argv)-1 != 2):
        print(f"Number of arguments {len(sys.argv)-1} is incorrect")
        usage()
        exit(1)

    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--blastx",type=str,required=True,help="Input: Centrifuge output")
    parser.add_argument("-o","--out",type=str,required=True,help="Output: Taxonomic results output")
    args = parser.parse_args()
    processCentrifuge(args.blastx,args.out)

if __name__ == "__main__":
    main()

