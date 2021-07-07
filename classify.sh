#!/bin/bash

nr_lines=$(gunzip -c $1|wc -l)
nr_reads=$((nr_lines/4))
file_size=$(wc -c $1)
file_size=${file_size% *}
max_size=${2:-10000000}
case ${2: -1} in
  G|g)
    max_size=$((1000000000*${2:0:-1}))
    ;;
  M|m)
    max_size=$((1000000*${2:0:-1}))
    ;;
  [0-9])
    ;;
  *)
    echo "Input error: $2"
    exit
    ;;
esac
if [ $file_size -le $max_size ]; then
  echo "No splitting required. The file size of ${1##*/}\
 is smaller than max_size."
  exit
fi
if [ ! $((file_size%max_size)) = '0' ]; then
  nr_files=$((file_size/max_size + 1))
else
  nr_files=$((file_size/max_size))
fi
split_size=$((4*(nr_reads/nr_files + nr_reads%nr_files)))
as='_'${1#*_}
gunzip -c $1|split -l $split_size --additional-suffix=${as%.gz} - ${1%%_*}\
'xxx'
rev=${1/R1/R2}
as='_'${rev#*_}
gunzip -c ${1/R1/R2}|split -l $split_size --additional-suffix=${as%.gz} - \
${rev%%_*}'xxx'
gzip *.fastq



