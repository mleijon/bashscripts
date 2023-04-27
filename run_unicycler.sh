#!/bin/bash

dir="/storage/micke/SA3/tempo/*"
odir="/storage/micke/SA3/unicycler/"
for f in $dir; do
   if [ -d "$f" ]; then
      unicycler -1 "$f"/"${f##*/}"_1.fastq.gz -2 "$f"/"${f##*/}"_2.fastq.gz -o "$odir"/"${f##*/}" --verbosity 0 
   fi
done
