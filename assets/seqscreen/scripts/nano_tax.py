#!/usr/bin/env python3

import os
from collections import defaultdict, Counter
import numpy as np
import pickle
import argparse
from ete3 import NCBITaxa


def reverse_map(uid_to_taxid):
    '''
    Reverse the dictionary output by build_dict(), so that taxIDs map to uniprotIDs.
    :param uid_to_taxid: A dictionary.
    :return: A dictionary
    '''
    reverse = defaultdict(set)
    for uid, taxids in uid_to_taxid.items():
        for taxid in taxids:
            reverse[taxid].add(uid)
    return reverse

def calculate_cost(reversed_dict, evalue_per_uid):
    '''
    For each taxid, take its cost to be the average evalue of the uids it maps to.
    :param
    :return: a dictionary that maps a taxid to its cost
    '''
    evalue_dict = {}

    for taxid in reversed_dict:
        list_of_evalues = list(map(lambda uid:evalue_per_uid[uid], reversed_dict[taxid]))
        average = sum(list_of_evalues) / len(list_of_evalues)
        evalue_dict[taxid] = average

    return evalue_dict

def solve_set_cover_greedy(universe_set, subsets, costs):
    """
        universe_set (list): Universe of elements UniProt
        subsets (dict): Subsets of Universe {'taxid1':{uid1,uid5},'taxid2':{uid3,ui7}}
        costs (dict): Costs of each subset in taxid - {'taxid1':cost, 'taxid2':cost...}
    """
    all_elements = set(e for s in subsets.keys() for e in subsets[s])
    # elements don't cover universe -> invalid input for set cover
    if all_elements != universe_set:
        print('Does not match.')
        return None

    # track elements of universe covered - UniProtID
    covered_set = set()
    cover_sets_solution = []

    while covered_set != universe_set:#seen the entire universe set
        min_cost_elem_ratio = float("inf") #initialize to infinity
        min_set = None
        # find set with minimum cost:elements_added ratio
        for s, elements in subsets.items():
            new_elements = len(elements - covered_set) #new elements in this set
            # set may have same elements as already covered -> new_elements = 0
            # check to avoid division by 0 error
            if new_elements != 0:
                cost_elem_ratio = costs[s] / new_elements #is this cost effective?
                if cost_elem_ratio < min_cost_elem_ratio: #if yes, then add
                    min_cost_elem_ratio = cost_elem_ratio #reset min cost
                    min_set = s # add the set
        cover_sets_solution.append(min_set)
        # union with set already covered_set
        covered_set |= subsets[min_set]
    return cover_sets_solution


def parse_btab(line):
    '''
    Given a line in a btab file corresponding to an orf, get the predicted tax and taxID.
    :param line:
    :return:
    '''
    split_line = line.split('\t')
    qseid_orf = split_line[0]
    qseid, orf = "_".join(qseid_orf.split("_")[:-1]), qseid_orf.split('_')[-1]
    UniRef_uid = split_line[1]
    evalue = float(split_line[10])
    salltitles = split_line[-1]
    temp_0 = salltitles.split('Tax=')[1]
    tax = temp_0.split('TaxID=')[0].strip()
    if len(tax.split()) > 1:
        tax_0 = tax.split()[0]
        tax_1 = tax.split()[1]
        tax = tax_0 + ' ' + tax_1
    temp_1 = temp_0.split('TaxID=')[1]
    taxid = temp_1.split('RepID=')[0].strip()

    return qseid, orf, tax, taxid, evalue, UniRef_uid


def analyze_btab(btab, output, full_dict, ncbi_taxa_db):
    info_of_qseids = defaultdict(list)
    evalue_of_uid = defaultdict(list) # Used to calculate the average evalue for each uniprotID
    relevant_dict = {} # Extract from the full dictionary items that are relevant to this sample. i.e. the uids that are hit.
    final_map = {} # This is the output
    with open(btab,"r") as f:
        for line in f:
            qseid, orf, tax, taxid, evalue, UniRef_uid= parse_btab(line)
            uid = UniRef_uid.split('_')[-1]
            # Take down the evalue corresponding to this uid
            evalue_of_uid[uid].append(evalue)
            # Build the relevant dict. Add the uniprotID to dict, if it is not there already
            if uid not in relevant_dict:
                if uid in full_dict:
                    #print(full_dict[uid])
                    relevant_dict[uid] = full_dict[uid]
                else:
                    relevant_dict[uid] = [taxid]
            # Record other info for this qseid
            info_of_qseids[qseid].append((tax, taxid, uid, evalue))

        # For each uid, calculate its average evalue
        evalue_per_uid = {}
        for uid in evalue_of_uid:
            average = sum(evalue_of_uid[uid]) / len(evalue_of_uid[uid])
            evalue_per_uid[uid] = average

        # Find the minset cover for all hit uniprotIDs
        # Extract the universal set. (1st arg)
        universe_set = {uid for uid in relevant_dict}
        count_tax = dict(Counter([t for uid in relevant_dict for t in relevant_dict[uid]]))
        count_prop = {k:v/sum(count_tax.values()) for k,v in count_tax.items()}
        cutoff = np.percentile(list(count_prop.values()),25)
        ignore_taxa = set([k for k,v in count_prop.items() if v <= cutoff])
        # Reverse the dictionary, get taxid-uid. (2nd arg)
        reversed_dict = reverse_map(relevant_dict)
        # Get the costs (3rd arg)
        costs = calculate_cost(reversed_dict, evalue_per_uid)
        minset = solve_set_cover_greedy(universe_set, reversed_dict, costs)
        minset_pruned = [t for t in minset if t not in ignore_taxa]
        
    for qseid in info_of_qseids:
        mixing_taxa = [int(tup[1]) for tup in info_of_qseids[qseid] if tup[1] != "N/A"]
        ranks = ncbi_taxa_db.get_rank(mixing_taxa)
        genus_taxa = [k for k in ranks if  ranks[k] == "genus"]
        species_taxa = [k for k in ranks if ranks[k] == "species"]

        genus_count = Counter(genus_taxa).most_common()
        pre_count = {species:species_taxa.count(species) for species in set(species_taxa)}
        # If "escherichia" appears as a genus hit, we also count it as e.coli
        for genus, num in genus_count:
            for species in pre_count:
                if ncbi_taxa_db.get_lineage(species)[-1] == genus:
                    pre_count[species] += num
        species_count = sorted(pre_count.items(), key=lambda x: x[1], reverse=True)

        if len(species_count):
            # If there is no simple majority species hit
            if (len(species_count) > 1 and species_count[0][1] == species_count[1][1]):
                # Filter out the hits that are not in the minset
                remaining = list(filter(lambda tup:tup[1] in minset_pruned, info_of_qseids[qseid]))
                # If there are hits that are in the minset, then we choose from them the hit with the lowest evalue
                if len(remaining):
                    prediction = min(remaining, key=lambda tup:tup[-1])[1]
                # If no hit is in the minset, then we choose the hit with the lowest evalue
                else:
                    prediction = min(info_of_qseids[qseid], key=lambda tup:tup[-1])[1]

            else: # has simple majority species hit
                prediction = species_count[0][0]

        else: # No species hits
            if genus_taxa:
                prediction = min(info_of_qseids[qseid], key=lambda tup:tup[-1])[1]
            else:
                prediction = 2 #bateria

        final_map[qseid] = prediction

    with open(os.path.join(output,"read_level_taxid.tsv"), 'w+') as f:
        f.write("Read\ttaxid\n")
        for _id, prediction in final_map.items():
            f.write(str(_id) + '\t' + str(prediction) + '\n')




def main():

    parser = argparse.ArgumentParser(prog="nano_tax")
    parser.add_argument("--blastx", type=str, required=True, help="Blastx btab file")
    parser.add_argument("--database", type=str, required=True, help="Database directory")
    parser.add_argument("--output", type=str, required=True, help="Taxonomic output directory")
    parseArgs = parser.parse_args()

    
    ete3db = os.path.join(parseArgs.database, "reference_inference", "taxa.sqlite")
    ncbi_taxa_db = NCBITaxa(dbfile=ete3db)
    
    full_dict = {}
    with open(os.path.join(parseArgs.database, "uniref_multitax.pickle"),"rb") as inpf:
        full_dict = pickle.load(inpf)
        
    analyze_btab(parseArgs.blastx, parseArgs.output, full_dict, ncbi_taxa_db)


if __name__ == "__main__":
    main()
