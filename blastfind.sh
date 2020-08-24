#!/bin/bash

FILES="$PWD/*.$1"
nr=0
nr_of_files=$(find -- *.$1|wc -l)
for b in ${FILES}; do
    nr=$((nr+1))
    echo -en 'Processing file: '${nr}'/'"${nr_of_files}"\\r
    f=$(dirname "$b")
    f=$(dirname "$f")/$(basename "$b")
    f=${f/$1/fa}
    blast_find.py -f "${f}" -b "${b}" -t 'virus';wait
done
mkdir -p top;mkdir -p deep;mkdir -p virus_fasta;wait
mv -- *top.txt top; mv -- *deep.txt deep; mv -- *virus.fa virus_fasta
wait; gzip *.$1
