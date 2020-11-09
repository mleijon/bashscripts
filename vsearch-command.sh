#!/bin/bash

FILES="/media/micke/extra/reningsverksdata/00-RawData-k/*R1*"
for f in $FILES; do
  vsearch --fastq_join "$f" --reverse "${f/R1/R2}" --log log-join.txt --no_progress --fastqout /dev/stdout|\
  vsearch --derep_fulllength - --output /dev/stdout --sizeout --no_progress --log log-derep.txt|\
 /home/micke/PycharmProjects/fasta-tools/vsearchjoin_split.py| gzip > test.fa.gz
done
