### I ran everything from the src directory

## To get results for Q1.a):
1. run `get_boundaries.sh input.file output.file strand_code`\
exact command I used is:\
`./get_boundaries.sh ../data/Vibrio_cholerae.GFC_11.37.gff3 cholerae/cholerae_bound.txt +`

2. to get average lengths of intergenic and genic regions, run `get_avg_len.py -g gff3.file -b bounds.file -o output.file`\
exact command I used is:\
`python get_avg_len.py -g ../data/Vibrio_cholerae.GFC_11.37.gff3 -b cholerae/cholerae_bound.txt -o cholerae/average_lengths.json`

3. to get nucleotide frequencies for intergenic regions and codon frequencies for genic regions, run `get_freqs.py -f fasta.file -b bounds.file -o output.file`\
exact command I used is:\
`python get_freqs.py -f ../data/Vibrio_cholerae.GFC_11.dna.toplevel.fa -b cholerae/cholerae_bound.txt -o cholerae/codon_frequencies.json`

## To get results for Q1.b & c):
run `python viterbi.py -f fasta.file -c config.json -o output.file`\
exact command I used is:\
`python viterbi.py -f ../data/Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa -c config.json -o vulnificus/predictions.gff3`

## To get results for Q1.d):
run `get_boundaries.sh input.file output.file strand_code` and\
`python eval_accuracy.py -b bounds.file -g gff3.file`\
exact commands I used are:\
`./get_boundaries.sh ../data/Vibrio_vulnificus.ASM74310v1.37.gff3 vulnificus/vulnificus_bound.txt +` and\
`python eval_accuracy.py -b vulnificus/vulnificus_bound.txt -g vulnificus/predictions.gff3`