#!/bin/bash

FILES="/home/micke/unicycler-run2/assemblies_renamed/*.fasta"
for f in $FILES; do
outname=${f##*/}
echo processing: $outname
outname=${outname%.fasta}
wait;rgi main -i $f -o $outname -t contig --exclude_nudge --clean;wait
done