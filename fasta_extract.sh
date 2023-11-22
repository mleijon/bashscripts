#!/bin/bash

clear
if [ $2 -ge $(grep -c '>' < $1) ]; then
   echo "Nr of sequence records in fasta file are less or equal to $2. No splitting."
   exit 0
fi
csplit -zs $1 /"$(grep -m$(( $2 + 1 )) '>' $1 |tail -n1)"/
mv xx00 $(awk -v file=$1 -v nr=$2 'BEGIN{sub(/\.fasta/, "_" nr ".fasta" , file); print file}') 
rm xx01
