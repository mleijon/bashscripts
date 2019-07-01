#!/usr/bin/env bash

FILES=$PWD/*.fa
nr=0
nr_of_files=$(ls *.fa|wc -l)
for f in ${FILES}; do
    nr=$((nr+1))
    echo -en 'Processing file: '${nr}'/'${nr_of_files}\\r
    blastn -task blastn -query ${f} -db /home/micke/blast_dbs/nt \
    -out ${f/fa/blastn} -evalue 1e-5 -word_size 11 -gapopen 5\
    -gapextend 2 -reward 2 -penalty -3 -num_descriptions 5 -num_alignments 0\
    -outfmt 0 -num_threads 8 -parse_deflines
done