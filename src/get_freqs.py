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
def pop_nuc_dict(string, nuc_dict, contig_name):
    nuc_dict['A'] += string.count('A')
    nuc_dict['C'] += string.count('C')
    nuc_dict['T'] += string.count('T')
    nuc_dict['G'] += string.count('G')
    # if contig_name == 'DN38.contig00082':
    #     print(string.count('A'))
    #     print(nuc_dict['A'])

def pop_nuc_dict_list(list_of_strings, nuc_dict, contig_name):
    for string in list_of_strings:
        pop_nuc_dict(string, nuc_dict, contig_name)

# deal with all intergenic regions in a contig
def get_intergenic(contig_bounds, contig_seq, nuc_dict, contig_name):
    # if whole contig is intergenic, populate nuc_dict with entire sequence
    if contig_bounds.empty:
        pop_nuc_dict(contig_seq, nuc_dict, contig_name)
        # print('no genic', nuc_dict)
        
    # if contig has genic regions
    else:
        all_inter = []

        # if contig_name == 'DN38.contig00082':
        #     print(contig_seq)
        #     print('empty',all_inter)

        # get beginning intergenic
        first_bound = contig_bounds.start[0]
        # add first inter chunk
        all_inter.append(contig_seq[0:first_bound-1])
        # get last intergenic
        last_bound = contig_bounds['end'].iloc[-1]
        # add last inter chunk
        all_inter.append(contig_seq[last_bound+1:])
        # if contig_name == 'DN38.contig00082':
        #     print(all_inter)

        # get all middle intergenic nucleotides and add to total
        for i in range(1, len(contig_bounds)):
            int_end, int_start = contig_bounds['start'].loc[i], contig_bounds['end'].loc[i-1]
            all_inter.append(contig_seq[int_start+1:int_end])
            # if contig_name == 'DN38.contig00082':
            #     print(int_start, int_end)
            #     print(all_inter)
        # then populate dictionary with all counts
        pop_nuc_dict_list(all_inter, nuc_dict, contig_name)
        # print(nuc_dict)

def get_genic(contig_bounds, contig_seq, start_dict, mid_dict, stop_dict, contig_name):
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

            # if contig_name == 'DN38.contig00082':
            #     print(start_codon, stop_codon)

            middle_chunk = gen_chunk[3:-3]

            # if contig_name == 'DN38.contig00082':
            #     print(middle_chunk)

            # go codon by codon
            for i in range(0,len(middle_chunk),3):
                codon = middle_chunk[i:i+3]
                # if contig_name == 'DN38.contig00082':
                #     print(codon)
                # update middle codon count
                mid_dict[codon] += 1

def count_to_freq(dict):
    # make counts into frequencies
    num_values = sum(dict.values())
    for key in dict:
        dict[key] = dict[key] / num_values


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

        # if contig_name == 'DN38.contig00082':
        #     print(contig_seq)

        contig_bounds = bounds_df.loc[bounds_df['contig']==contig_name, ['start','end']].reset_index()

        # print(contig_name)

        # get intergenic regions and populate nuc_dict
        get_intergenic(contig_bounds, contig_seq, nuc_dict, contig_name)
        
        # get genic regions and populate codon_dict
        get_genic(contig_bounds, contig_seq, start_dict, mid_dict, stop_dict, contig_name)

    # make counts into frequencies
    count_to_freq(nuc_dict)
    count_to_freq(start_dict)
    count_to_freq(mid_dict)
    count_to_freq(stop_dict)


    # output to file :D
    with open(output_file, 'w') as file:
        sys.stdout = file
        json.dump(nuc_dict, file, indent=4)
        print('\nstart codon frequencies')
        json.dump(start_dict, file, indent=4)
        print('\nmiddle codon frequencies')
        json.dump(mid_dict, file, indent=4)
        print('\nstop codon frequencies')
        json.dump(stop_dict, file, indent=4)



if __name__ == '__main__':
    main()