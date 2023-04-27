#!/bin/bash
FILES="/storage/micke/fsk/ipn/spades/*/contigs.fasta"
for f in $FILES; do
outname=${f##*spades/}
outname=${outname%-*}
outdir=/storage/micke/fsk/ipn/diamond_out/
wait;diamond blastx -d /ssd2/classify/nr -q $f -o $outdir/$outname.daa --max-target-seqs 1 --evalue 1E-5 --outfmt 6 \
qseqid full_qseq evalue staxids sscinames sskingdoms skingdoms sphylums -b 20 -c 1 --compress 0;wait
done
