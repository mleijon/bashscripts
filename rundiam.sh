#!/bin/bash

FILES="/storage/micke/msbh2o/HTStream_processed/svartvattnet_kallan/*_R*.fastq.gz"
for f in $FILES; do
outname=${f##*/}
outname=${outname%%_*}
outdir=/storage/micke/msbh2o/HTStream_processed/svartvattnet_kallan/out_diamond
wait;diamond blastx -d /ssd2/diamondDB/nr -q $f -o $outdir/$outname.daa --max-target-seqs 5 --evalue 1E-5 --outfmt 102 -b 20 -c 1 --compress 1;wait
done
