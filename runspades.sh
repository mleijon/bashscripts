#!/bin/bash

FILES="/home/fou/Documents/Emily/846394/*/*_1.fastq.gz"
for f in $FILES; do
name=${f##*HJM7LDSXX_}
outdir=$PWD/${name%%_*}/spades
wait;spades.py -o $outdir -1 $f -2 ${f/_1.fastq/_2.fastq} -t 72 -m 450; sleep 1m;sync;
done
