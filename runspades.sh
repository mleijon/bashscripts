#!/bin/bash

FILES="/data/micke/data_evaluation/kc-ranch/*_R1_001.fastq.gz"
for f in $FILES; do
name=${f##*/}
name=${name%%_*}
path=${f%/*}/
outdir=$path/$name/spades
wait;spades.py -o $outdir -1 $f -2 ${f/_R1_/_R2_} -t 48 -m 160; sleep 1m;sync;
done
