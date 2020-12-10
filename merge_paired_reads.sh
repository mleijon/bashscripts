#!/usr/bin/env bash


#for f in $dirs;do cd $f;./*.sh; cd ../..;done

FILES="/storage/micke/ku-kc-ranch/01-Cleaned/*R1.fastq.gz"
for f in $FILES; do
cat ${f/R1/R2} >> $f; wait
rm ${f/R1/R2}; wait
rename 's/_R1.fastq/.fastq/' $f; wait
done
