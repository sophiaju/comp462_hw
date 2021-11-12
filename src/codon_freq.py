import os, sys
import argparse, json
import pandas as pd
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta')
    parser.add_argument('-c','--cleaned')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    cleaned = args.cleaned
    fasta_file = args.fasta
    output_file = args.output

    df = pd.read_csv(cleaned, names=['contig', 'start', 'end'])
    
    # initialize dict with all possible codons
    codons = {}
    nucleotides = ['A','T','G','C']
    for nuc1 in nucleotides:
        for nuc2 in nucleotides:
            for nuc3 in nucleotides:
                codons[nuc1+nuc2+nuc3] = 0
    
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(fasta_file), 'fasta'))
    indexes = {}
    for key in fasta_sequences:
        temp = df.loc[df['contig']==key,['start','end']].reset_index()
        indexes[key] = []
        for i in range(temp.shape[0]):
            indexes[key].append((int(temp.loc[i,'start']),int(temp.loc[i,'end'])))

    sequences = {}
    for key in fasta_sequences:
        sequence = ''
        for i in range(len(fasta_sequences[key])):
            sequence+=fasta_sequences[key][i]
        sequences[key] = sequence

    total_codons = 0
    total_start_codons = 0
    total_stop_codons = 0
    start_codons = {}
    stop_codons = {
        'TAA': 0,
        'TAG': 0,
        'TGA': 0
    }
    
    for key in sequences:
        for a,b in indexes[key]:
            for i in range(a,b,3):
                # correct the reading frame
                codon = sequences[key][i-1:i+2]

                # start of the genic region
                if i == a:
                    if codon not in start_codons:  
                        start_codons[codon] = 1
                    else:
                        start_codons[codon] += 1
                    total_start_codons +=1

                # do stop codons
                elif codon in stop_codons:
                    stop_codons[codon] +=1
                    total_stop_codons +=1
                else:
                    codons[codon]+=1
                    total_codons +=1

    # calculate the frequencies
    for key in codons:
        codons[key] = codons[key] / total_codons
    for key in start_codons:
        start_codons[key] = start_codons[key] / total_start_codons
    for key in stop_codons:
        stop_codons[key] =stop_codons[key] / total_stop_codons    

    with open(output_file, 'w') as fp:
        json.dump(codons, fp, indent=4)
    with open(os.path.join('data', 'start_codon_freq.json'), 'w') as fp:
        json.dump(start_codons, fp, indent=4)
    with open(os.path.join('data', 'stop_codon_freq.json'), 'w') as fp:
        json.dump(stop_codons, fp, indent=4)
if __name__ == '__main__':
    main()