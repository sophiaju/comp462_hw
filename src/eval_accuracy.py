import os, sys
import argparse
import json
from Bio import SeqIO
import math
import pandas as pd

def get_match_counts(first_group, sec_group, dict):

    for first_contig in first_group:

        first_contig_name = first_contig[0]
        first_contig_df = first_contig[1]

        # if the first contig is also in the second file
        if first_contig_name in sec_group.groups:
            sec_contig = sec_group.get_group(first_contig_name)

            # get all start and end positions as lists
            sec_start_list = sec_contig['start'].to_list()
            sec_end_list = sec_contig['end'].to_list()

            # for each gene in first group
            for index, row in first_contig_df.iterrows():
                a_start = row['start']
                a_end =  row['end']

                # matches neither
                if a_start not in sec_start_list and a_end not in sec_end_list:
                    dict['none'] += 1
                # matches start but not end
                elif a_start in sec_start_list and a_end not in sec_end_list:
                    dict['start'] += 1
                # matches end but not start
                elif a_start not in sec_start_list and a_end in sec_end_list:
                    dict['end'] += 1
                # matches both
                else:
                    dict['both'] += 1

        # if first contig not in second, no matches
        else:
            # print('no group', first_contig_name)
            # add all first genes to no matches count
            dict['none'] += len(first_contig)

def count_to_freq(dict):
    # make counts into frequencies
    num_values = sum(dict.values())
    for key in dict:
        dict[key] = dict[key] / num_values

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bounds', help='add the annotated bounds file for prediction')
    parser.add_argument('-g','--gff3', help='add the gff3 file with predictions')
    args = parser.parse_args()
    bounds_file = args.bounds
    gff3_file = args.gff3

    annotated_df = pd.read_csv(bounds_file, names=['contig', 'start', 'end'], sep=' ')
    predicted_df = pd.read_csv(gff3_file, names=['contig',2, 'start', 'end', 5, 6, 7, 8], sep='\t')
    predicted_df.drop(columns=[2,5,6,7,8], inplace=True)

    annot_gr = annotated_df.groupby('contig')
    pred_gr = predicted_df.groupby('contig')

    # for each contig in the annotated (given) file
    a_matches = {'both':0, 'start':0, 'end':0, 'none':0}

    get_match_counts(annot_gr, pred_gr, a_matches)
    count_to_freq(a_matches)
    print('fraction of annotated genes that: \nperfectly match both ends of one of predicted genes: ', a_matches['both'])
    print('match the start but not the end of a predicted gene:', a_matches['start'])
    print('match the end but not the start of a predicted gene:', a_matches['end'])
    print('do not match either the start or end of a predicted gene:', a_matches['none'])

    print('\n')

    # same with predicted genes
    p_matches = {'both':0, 'start':0, 'end':0, 'none':0}
    get_match_counts(pred_gr, annot_gr, p_matches)
    count_to_freq(p_matches)
    print('fraction of predicted genes that: \nperfectly match both ends of one of annotated genes: ', p_matches['both'])
    print('match the start but not the end of a annotated gene:', p_matches['start'])
    print('match the end but not the start of a annotated gene:', p_matches['end'])
    print('do not match either the start or end of a annotated gene:', p_matches['none'])




if __name__ == '__main__':
    main()