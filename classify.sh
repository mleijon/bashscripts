#!/bin/bash

if [ ! -d "$1" ]; then
  echo "Please enter the full path to the directory containing the sequencing\
   data as argument."
  exit 0
else
  FILES=$1'/*R1*.gz'
fi
dir=$(dirname "$FILES")
for f in $FILES; do
  base=$(basename "$f")
  mkdir ${base%%_*}
  nr_lines=$(gunzip -c $f|wc -l)
  nr_reads=$((nr_lines/4))
  file_size=$(wc -c $f)
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
    echo "No splitting required. The file size of ${f##*/}\
   is smaller than max_size."
    exit
  fi
  if [ ! $((file_size%max_size)) = '0' ]; then
    nr_files=$((file_size/max_size + 1))
  else
    nr_files=$((file_size/max_size))
  fi
  split_size=$((4*(nr_reads/nr_files + nr_reads%nr_files)))
  as='_'${f#*_}
  gunzip -c $f|split -l $split_size --additional-suffix=${as%.gz} - \
  ${f%%_*}'xxx'
  rev=${f/R1/R2}
  as='_'${rev#*_}
  gunzip -c ${f/R1/R2}|split -l $split_size --additional-suffix=${as%.gz} - \
  ${rev%%_*}'xxx'
  gzip *.fastq
  mv *xxx* ${base%%_*}
done



