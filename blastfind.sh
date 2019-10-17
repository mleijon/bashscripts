#!/bin/bash

FILES="$PWD/*.blastn"
nr=0
nr_of_files=$(find -- *.blastn|wc -l)
for b in ${FILES}; do
    nr=$((nr+1))
    echo -en 'Processing file: '${nr}'/'"${nr_of_files}"\\r
    f=$(dirname "$b")
    f=$(dirname "$f")/$(basename "$b")
    f=${f/blastn/fa}
    blast_find.py -f "${f}" -b "${b}" -t 'virus';wait
done
mkdir top;mkdir deep;mkdir virus_fasta;wait
mv -- *top.txt top; mv -- *deep.txt deep; mv -- *virus.fa virus_fasta
