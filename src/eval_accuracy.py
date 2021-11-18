import os, sys
import argparse
import json
from Bio import SeqIO
import math
import pandas as pd


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bounds', help='add the annotated bounds file for prediction')
    parser.add_argument('-g','--gff3', help='add the gff3 file with predictions')
    parser.add_argument('-o', '--output', help='add output file')
    args = parser.parse_args()
    bounds_file = args.bounds
    gff3_file = args.gff3
    output_file = args.output

    bounds_df = pd.read_csv(bounds_file, names=['contig', 'start', 'end'], sep=' ')
    gff3_df = pd.read_csv(gff3_file, names=['contig',2, 'start', 'end', 5, 6, 7, 8], sep='\t')
    gff3_df.drop(columns=[2,5,6,7,8], inplace=True)

    print(bounds_df)
    print(gff3_df)

    # annotated (given) genes that

    # perfectly match both ends of one of predicted

    # match the start but not the end
     
    # match the end but not the start

    # Do not match either the start or end

    # same with predicted genes




if __name__ == '__main__':
    main()