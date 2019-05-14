#!/bin/bash

if [[ -z $1 ]]; then
	echo "Use: 'splitfasta.sh gz' or 'splitfasta.sh fastq'"
	exit 0
fi 
FILES=$PWD/*.$1
for f in ${FILES}; do
	if [[ ! -f  "$f" ]]; then
		echo "No $1-files in the current directory!"
		exit 0
	fi
done
nr=0
nr_of_files=$(ls *.$1|wc -l)
for f in ${FILES}; do
    nr=$((nr+1))
	echo -en 'Processing file: '${nr}'/'${nr_of_files}\\r
	if [[ "$1" == "gz" ]]; then
		gunzip -c ${f} > ${f/.gz/}; wait
		u=${f/.gz/}
	else
		u=${f}
	fi
	b=$(dirname "$f")/vrl_blastn/$(basename "$f")
	b=${b/fastq/blastn}
	if [[ "$1" == "gz" ]]; then
	    gunzip -c ${b} > ${b/.gz/}; wait
	    b=${b/.gz/}
	fi
	split_fasta.py -b ${b} -f ${u};wait
	basename=${u##*/}
	hitname=$PWD'/hits_'${basename/fastq/fa}
	mv ${hitname} $PWD/vrlhits
	gzip $PWD/vrlhits/$(basename "$hitname");wait
	rm $PWD'/nohits_'${basename/fastq/fa} $PWD'/'${basename/fastq/fa};wait
	rm ${b} ${u}; wait
done
