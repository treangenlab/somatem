#!/usr/bin/env python3

'''
This module provides functionality for splitting input sequences into fragments, with a default length of 200bp and
overlap of 100bp.

Functions:
    parse_command_args(): Parses command-line arguments and returns them as a namespace object.
    fastadict_from_fasta_lines_iterator(fasta_file_line_iter): Reads a FASTA file from an iterator and returns a
        dictionary of sequences.
    fastadict_read_from_file(file_path, with_key_list=False): Reads a FASTA file and returns a dictionary of sequences.
    fastadict_write_to_file(out_file_path, fasta_dict, subset_keys=None, raw=False, quiet=False): Writes a dictionary
        of sequences to a FASTA file.
    split_fasta_file(input_file_path, output_file_path, fragment_length, overlap_window_length): Splits sequences in a
        FASTA file into fragments and writes the output to a new FASTA file.
    main(args): Main function to execute the script based on command-line arguments.

Created by: Mike Nute (10/13/2024)

Copyright 2024, Michael Nute, All rights reserved.
Free for research use. For commercial use, contact
Dr. Todd Treangen (treangen@rice.edu)

'''

import os, sys, argparse, logging

parser = argparse.ArgumentParser()
rich_format = "[%(filename)s (%(lineno)d) %(asctime)s] %(levelname)s: %(message)s"
logging.basicConfig(format=rich_format, level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)

def parse_command_args():
    global parser

    # subparsers = parser.add_subparsers(help='specifies the action to take.')
    # megahit_stats = subparsers.add_parser('megahit_stats',
    #                                      help='Goes through a megahit results folder and gets the summary statistics '
    #                                           'from the `final.contigs.fa` output file.')

    parser.add_argument('-i', '--fasta', dest='input_file_path', type=str, required=True,
                              help='Path to the fasta file containing the intput sequences.')
    parser.add_argument('-o', '--output', dest='output_file_path', type=str, default='',
                               help='Delimiter to use in separating output fields. Defualt is the same as '
                                    'the input file with the suffix "_split_{fragment_length}bp" before the '
                                    'file extension.')
    parser.add_argument('-fo', '--working_dir', dest='working_dir', type=str, default='',
                               help='Path to the working directory where the output fasta will be written. If \'--output\' '
                                    'is given, this argument is ignored.')
    parser.add_argument('-l', '--length', dest='fragment_length', type=int, default=200,
                               help='Length of output sequences.')
    parser.add_argument('-w', '--overlap', dest='overlap_window_length', type=int, default=100,
                               help='Length of the overlap window to use when splitting the sequences.')
    parser.set_defaults(func=main)

    args = parser.parse_args()

    if args.output_file_path == '':
        input, input_ext = os.path.splitext(args.input_file_path)
        if args.working_dir != '':
            args.output_file_path = os.path.join(args.working_dir, os.path.basename(input) + f'_split_{args.fragment_length}bp' + input_ext)
        else:
            args.output_file_path = input + f'_split_{args.fragment_length}bp' + input_ext

    return args


def fastadict_from_fasta_lines_iterator(fasta_file_line_iter):
    '''
    Essentially reads a fasta file and parses it into a dictionary object (a fastadict, specifically),
    although technically the input here is just an iterator over the lines of a fasta file. That means we
    can stick this in the fastadict_read_from_file function as well as the fastadict_from_lines_list method.

    Args:
        fasta_file_line_iter (Iterator of str): an iterator of strings representing the lines of a fasta file.

    Returns:
        fasta_dict (dict): The Fasta-Dictionary (dict keyed by sequence names, valued by sequence).
    '''
    output={}
    first=True
    seq=''
    name = None
    for l in fasta_file_line_iter:
        if l[0]=='>':
            if first!=True:
                output[name]=seq
            else:
                first=False
            name=l[1:].strip()
            seq=''
        else:
            seq=seq + l.strip()
    if name is not None:
        output[name]=seq
    return output

def fastadict_read_from_file(file_path, with_key_list = False):
    '''
    Reads from a fasta file and returns a dictionary where keys are taxon names and values are strings
    representing the sequence. Does NOT accept .gz files.

    Args:
        file_path (str):    The full system-readable path to the fasta file
        with_key_list (bool, optional):   If given, Return value is a tuple of (Fasta-Dict, Key-List) instead of just
                                the Fasta Dict (default: False).

    Returns:
        fasta_dict (dict):  dict keyed by sequence names, valued by sequence.
    '''
    with open(file_path,'r') as fasta:
        output = fastadict_from_fasta_lines_iterator(fasta)

    if not with_key_list:
        return output
    else:
        return output, list(output.keys())

def fastadict_write_to_file(out_file_path, fasta_dict, subset_keys=None, raw=False, quiet=False):
    '''
    Writes a Fasta-Dictionary object to a fasta file (output file given first).

    Takes a fasta dictionary (keys=taxon names, values=alignment strings) and writes it to a file. Contains
    optional variables to specify only a subset of keys to write, or to write strings without blanks (assumed
    to be '-').

    Args:
        out_file_path (str):    system-readable path to write to
        fasta_dict (:obj:`dict`):      Fasta-Dictionary to write
        subset_keys (:obj:`list`, optional): If given, should be a list of keys such that ONLY those keys will be
                                    written to the output file.
        raw (bool, optional):   If True, sequences are written to the output without any gaps.
        quiet (bool):           if True, does not output messages to the console.

    '''

    if subset_keys==None:
        mykeys=fasta_dict.keys()
    else:
        mykeys=list(set(subset_keys).intersection(fasta_dict.keys()))
        leftover_keys=list(set(subset_keys).difference(fasta_dict.keys()))
        if not quiet:
            logging.INFO('There were ' + str(len(leftover_keys)) + ' keys in the subset that were not in the original data.\n')

    fasta=open(out_file_path,'w')
    for i in mykeys:
        fasta.write('>'+i+'\n')
        if raw==False:
            fasta.write(fasta_dict[i]+'\n')
        else:
            fasta.write(fasta_dict[i].replace('-','')+'\n')
    fasta.close()
    if not quiet:
        logging.INFO('wrote file: ' + out_file_path + ' with the specified keys.')

def split_fasta_file(input_file_path, output_file_path, fragment_length, overlap_window_length):
    '''
    Splits a fasta file into fragments of a given length, with a given overlap window.

    Args:
        input_file_path (str):  Path to the fasta file containing the intput sequences.
        output_file_path (str): Path to the output file, fed straight to np.savez(...)
        fragment_length (int):  Length of output sequences.
        overlap_window_length (str): Length of the overlap window to use when splitting the sequences.
    '''
    fasta_dict = fastadict_read_from_file(input_file_path)
    fasta_new = {}

    for seq_name, seq in fasta_dict.items():
        seq_len = len(seq)
        num_fragments = seq_len // overlap_window_length
        # if seq_len % overlap_window_length > 0:
        #     num_fragments += 1
        start=0; end = min(fragment_length, seq_len);
        i=0;
        fragment = seq[start:end]
        fragment_name = f'{seq_name}_fragment_{i + 1}_{start}-{end}'
        fasta_new[fragment_name] = fragment
        while end < seq_len:
            start = end-overlap_window_length
            end = min(start + fragment_length, seq_len)
            if end == seq_len:
                start = end - fragment_length
            i += 1
            fragment = seq[start:end]
            fragment_name = f'{seq_name}_fragment_{i+1}_{start}-{end}'
            fasta_new[fragment_name] = fragment

    fastadict_write_to_file(output_file_path, fasta_new, quiet=True)
    print(output_file_path)

def main(args):
    split_fasta_file(args.input_file_path, args.output_file_path, args.fragment_length, args.overlap_window_length)

if __name__ == '__main__':
    cmd_args = parse_command_args()
    cmd_args.func(cmd_args)
