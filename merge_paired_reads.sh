#!/usr/bin/env bash


#for f in $dirs;do cd $f;./*.sh; cd ../..;done

FILES="/storage/micke/SZ-bats/01-Cleaned/*R1.fastq.gz"
for f in $FILES; do
cat ${f/R1/R2} >> $f; wait
cat ${f/R1/SE} >> $f; wait
rm ${f/R1/R2}; wait
rm ${f/R1/SE}; wait
gunzip $f; wait
rename _R1.fastq .fastq ${f/.gz/}; wait
done


for i in a b; do
 for j in a b; do
 echo $j;
 done;
done