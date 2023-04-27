#!/bin/bash

FILES="/storage/micke/fsk/ipn/*1_001.fastq.gz"
for f in $FILES; do
name=${f##*/ipn/}
outdir=/storage/micke/fsk/ipn/spades/${name%%_L001_*}
wait;spades.py -o $outdir -1 $f -2 ${f/1_001.fastq.gz/2_001.fastq.gz} -t 72 -m 350; sleep 1m;sync;
done
