#!/usr/bin/env bash
BLASTDB='/storage/micke/gb_release_254/'
export  BLASTDB
FILES=( "${PWD}"/*.fa )
for f in "${FILES[@]}"; do
    echo -en 'Processing file: '"${f}"\\r
    blastn -task megablast -query "${f}" -db /storage/micke/gb_release_254/nt -out "${f/fa/blastn}" -evalue 1e-10 \
    -word_size 28 -reward 1 -penalty -2 -num_threads 72 -outfmt "6 qseqid  sacc evalue sscinames scomname"  \
    -num_alignments 5 -subject_besthit; wait
done
