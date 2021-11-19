import pandas as pd
import numpy as np
import json, argparse

def main():

    # read in command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff3', help='add the gff3 file with contig lengths')
    parser.add_argument('-b','--bounds', help='add the file with gene region boundaries')
    parser.add_argument('-o', '--output', help='add output file')
    args = parser.parse_args()
    gff3_file = args.gff3
    bounds = args.bounds
    output_file = args.output

    # read in boundaries as dataframe
    df = pd.read_csv(bounds, sep = ' ', names=['contig','start','end'], dtype={'contig': str, 'start':int, 'end':int})

    # read in lengths of each contig region
    df_gff = pd.read_csv(gff3_file, sep = ' ', skiprows=[0], nrows=151, names=[0,1,2,'contig','start','end'])
    # for analysis question
    # df_gff = pd.read_csv(gff3_file, sep = ' ', skiprows=[0], nrows=388, names=[0,1,2,'contig','start','end'])
    df_gff.drop(columns=[0,1,2], inplace=True)
    # print(df_gff)
    contig_len_dict=df_gff.set_index('contig').to_dict('index')

    # avg lengths dict
    avg_len = {}

    # calculate length of genic region
    df['len'] = df['end'] - df['start']

    # find average length of genic region
    avg_genic = sum(df['len']) / len(df['len'])
    avg_len['genic'] = avg_genic
    print("average genic region length: ", avg_genic)

    # calculate length of intergenic region
    num_intergenic = 0
    sum_intergenic = 0

    # print(contig_len_dict.keys())

    # group dataframe by contig
    grouped = df.groupby(by='contig')
    for contig in grouped:
        contig_name = contig[0]

        # get the length of each contig (how many nucleotides in the full contig)
        contig_len = contig_len_dict[contig_name]['end']
        contig_df = contig[1]
   
        # append contig start at 1 and length to dataframe
        contig_df = contig_df.append({'contig':contig_name, 'start':contig_len, 'end':1}, ignore_index=True)
        # shift the end column down one (wrapping around)
        contig_df['end'] = np.roll(contig_df['end'],1)
        # calculate the intergenic lengths for each section
        contig_df['inter_len'] = contig_df['start'] - contig_df['end']

        sum_intergenic += sum(contig_df['inter_len'])
        num_intergenic += len(contig_df['inter_len'])

    # print(num_intergenic, sum_intergenic)

    # get lengths of contigs with no genic regions and add to intergenic length
    grouped_keys = [key for key, _ in grouped]
    for contig_key in contig_len_dict:

        # if the contig has no genic regions, ie not in grouped dataframe
        if contig_key not in grouped_keys:
            # no_genic.append(contig_key)
            num_intergenic += 1
            sum_intergenic += contig_len_dict[contig_key]['end']
            # print(contig_len_dict[contig_key]['end'])
    # print(num_intergenic, sum_intergenic)


    # print(num_intergenic, sum_intergenic)
    avg_inter = sum_intergenic / num_intergenic
    avg_len['intergenic'] = avg_inter
    print("average intergenic region length: ", avg_inter)

    with open(output_file, 'w') as file:
        json.dump(avg_len, file)


if __name__ == '__main__':
    main()
