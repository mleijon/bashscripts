#!/bin/bash

FILES="/storage/micke/ku-kc-ranch/01-Cleaned/*.fastq.gz"
for f in $FILES; do
outname=${f##*/}
outname=${outname%%_*}
outdir=/storage/micke/ku-kc-ranch/01-Cleaned/diamond_out
wait;diamond blastx -d /ssd2/diamondDB/nr -q $f -o $outdir/$outname.daa --max-target-seqs 5 --evalue 1E-5 --outfmt 102 -b 12 -c 1 --compress 1;wait
done
