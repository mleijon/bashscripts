#!/bin/bash

FILES="/storage/MSBvatten/msbvattenJS/SE-2194/190822_A00181_0108_AHLFMYDSXX/Sample_SE-2194-Botvid-Granby-filter/*R1*"
for f in $FILES; do
  outname=${f##*/}
  outname=${outname%%_*}
  outdir=/storage/micke/msbh2o/out_diamond
  vsearch --fastq_join "$f" --reverse "${f/R1/R2}" --log $outdir/$outname-join.log --no_progress --fastqout /dev/stdout|\
  vsearch --derep_fulllength - --output /dev/stdout --sizeout --no_progress --log $outdir/$outname-derep.log|\
  /home/micke/PycharmProjects/fasta-tools/vsearchjoin_split.py|\
  diamond blastx -d /ssd2/diamondDB/nr -q "$f" -o $outdir/$outname.daa -k5 -e1E-5 -f102 -b20 -c1 --compress 1;wait
done