To get results for Q1.a):
1. run get_boundaries.sh input.file output.file \
exact command I used is: \ 
./get_boundaries.sh Vibrio_cholerae.GFC_11.37.gff3 cholerae_bound.txt

2. to get average lengths of intergenic and genic regions, run get_avg_len.py -g gff3.file -b bounds.file -o output.file \
exact command I used is: \
python get_avg_len.py -g Vibrio_cholerae.GFC_11.37.gff3 -b cholerae_bound.txt -o average_lengths.json

3. to get nucleotide frequencies for intergenic regions and codon frequencies for genic regions, run get_freqs.py -f fasta.file -b bounds.file -o output.file  \ 
exact command I used is: \
python get_freqs.py -f Vibrio_cholerae.GFC_11.dna.toplevel.fa -b cholerae_bound.txt -o codon_frequencies.json