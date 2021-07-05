#!/bin/bash

nr_lines=$(gunzip -c $1|wc -l)
nr_reads=$((nr_lines/4))
file_size=$(wc -c $1)
file_size=${file_size% *}
max_size=${2:-10000000}
if [ ! $((file_size%max_size)) = '0' ]; then
  nr_files=$((file_size/max_size + 1))
else
  nr_files=$((file_size/max_size))
fi
split_size=$((4*(nr_reads/nr_files)))
gunzip -c $1|split -d --additional-suffix='.fastq' --suffix-length=3 -l \
$split_size - ${1%'.fastq.gz'}.'split'
gzip *.fastq


