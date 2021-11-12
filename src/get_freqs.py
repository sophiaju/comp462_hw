import os, sys
import argparse, json
from Bio.SeqIO.SffIO import _get_read_region
import pandas as pd
import numpy as np
from Bio import SeqIO

# initialize codon dictionary
def init_codons(codon_dict):
    nucleotides = ['A','T','G','C']
    for nuc1 in nucleotides:
        for nuc2 in nucleotides:
            for nuc3 in nucleotides:
                codon_dict[nuc1+nuc2+nuc3] = 0

# populate the nucleotide count dictionary
def pop_nuc_dict(string, nuc_dict):
    nuc_dict['A'] += string.count('A')
    nuc_dict['C'] += string.count('C')
    nuc_dict['T'] += string.count('T')
    nuc_dict['G'] += string.count('G')

# deal with all intergenic regions in a contig
def get_intergenic(contig_bounds, contig_seq, nuc_dict):
    # if whole contig is intergenic, populate nuc_dict with entire sequence
    if contig_bounds.empty:
        pop_nuc_dict(contig_seq, nuc_dict)
    # if contig has genic regions
    else:
        all_inter = ''

        # get beginning intergenic
        first_bound = contig_bounds.start[0]
        # add first inter chunk
        all_inter += contig_seq[0:first_bound-1]
        # get last intergenic
        last_bound = contig_bounds['end'].iloc[-1]
        # add last inter chunk
        all_inter += contig_seq[last_bound+1:]

        # get all middle intergenic nucleotides and add to total
        for i in range(1, len(contig_bounds)):
            int_start, int_end = contig_bounds['start'].loc[i], contig_bounds['end'].loc[i-1]
            all_inter += contig_seq[int_start+1:int_end]

        # then populate dictionary with all counts
        pop_nuc_dict(all_inter, nuc_dict)
        # print(nuc_dict)

def get_genic(contig_bounds, contig_seq, start_dict, mid_dict, stop_dict):
    if not contig_bounds.empty:
        # for each genic region
        for i in range(0, len(contig_bounds)):
            gen_start, gen_end = contig_bounds['start'].loc[i], contig_bounds['end'].loc[i]
            gen_chunk = contig_seq[gen_start-1: gen_end]

            # update start codon count
            start_codon = gen_chunk[0:3]
            start_dict[start_codon] += 1
            # update stop codon count
            stop_codon = gen_chunk[-3:]
            stop_dict[stop_codon] += 1

            middle_chunk = gen_chunk[3:-3]
            # go codon by codon
            for i in range(0,len(middle_chunk),3):
                codon = middle_chunk[i:i+3]

                # update middle codon count
                mid_dict[codon] += 1


def main():
    # read in command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help='add the fasta file to be analyzed')
    parser.add_argument('-b','--bounds', help='add the file with gene region boundaries')
    parser.add_argument('-o', '--output', help='add output file')
    args = parser.parse_args()
    fasta_file = args.fasta
    bounds = args.bounds
    output_file = args.output

    # parse in bounds as df and fasta sequences as dictionary
    bounds_df = pd.read_csv(bounds, names=['contig', 'start', 'end'], sep=' ')
    fasta_dict = SeqIO.to_dict(SeqIO.parse(open(fasta_file), 'fasta'))

    # initialize empty codon count and nucleotide count dictionaries
    start_dict = {}
    mid_dict = {}
    stop_dict = {}
    init_codons(start_dict)
    init_codons(mid_dict)
    init_codons(stop_dict)
    nuc_dict = {'A':0, 'C':0, 'G':0, 'T':0}

    
    # for each contig, populate codon_dict and nuc_dict
    for contig in fasta_dict:
        contig_name = fasta_dict[contig].id
        contig_seq = fasta_dict[contig].seq

        contig_bounds = bounds_df.loc[bounds_df['contig']==contig_name, ['start','end']].reset_index()

        # get intergenic regions and populate nuc_dict
        get_intergenic(contig_bounds, contig_seq, nuc_dict)
        

        # get genic regions and populate codon_dict
        get_genic(contig_bounds, contig_seq, start_dict, mid_dict, stop_dict)



    # print(fasta_dict["DN38.contig00001"].seq[653:656])

    print(json.dumps(nuc_dict, indent=4))
    print(json.dumps(start_dict, indent=4))
    print(json.dumps(mid_dict, indent=4))
    print(json.dumps(stop_dict, indent=4))



if __name__ == '__main__':
    main()