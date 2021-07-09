#!/bin/bash

shopt -s nocaseglob

#######INPUT DATA CONTROL########
max_size=${2:-50000000}
if [ ! -d "$1" ]; then
  echo "Please enter the full path to the directory containing the sequencing\
data as argument."
  exit 0
else
  FILES=($1/*R1*.gz)
  if [ ! -f ${FILES[0]} ]; then
    echo "No sequence files in: $1"
    exit
  fi
  dir=$(dirname "$FILES")
fi
if [[ "${2:0:-1}" =~ ^[0-9]+$ ]]; then
  :
elif [ -z $2 ]; then
  :
else
  echo "Input error: $2"
  exit
fi
if [ ! -z $2 ]; then
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
fi
##################################

for f in ${FILES[@]}; do
  base=$(basename "$f")
  ext=${f##*.}
  nr_lines=$(gunzip -c $f|wc -l)
  nr_reads=$((nr_lines/4))
  file_size_fwd=$(wc -c $f)
  file_size_fwd=${file_size_fwd% *}
  if [ ! -f ${f/R1/R2} ]; then
    echo "No reverse file: ${f/R1/R2}. Skipping this sample."
    continue
  fi
  file_size_rev=$(wc -c ${f/R1/R2})
  file_size_rev=${file_size_rev% *}
  if [ $file_size_fwd -lt $file_size_rev ]; then
    file_size=$file_size_rev
  else
    file_size=$file_size_fwd
  fi

  if [ $file_size -le $max_size ]; then
    echo "No splitting required for: ${f/R1/R*}. Copying files..."
    if [ -d $dir/${base%%_*} ]; then
      rm -r $dir/${base%%_*}
    fi
    mkdir $dir/${base%%_*}
    cp $f $dir/${base%%_*}
    cp ${f/R1/R2} $dir/${base%%_*}
    continue
  fi

  if [ ! $((file_size%max_size)) -eq 0 ]; then
    nr_files=$((file_size/max_size + 1))
  else
    nr_files=$((file_size/max_size))
  fi
  if [ $nr_files -gt 100 ]; then
    echo "Too many splits: $nr_files (maximum: 100)"
    exit
  fi
  split_size=$((4*(nr_reads/nr_files + nr_reads%nr_files)))
  as='_'${f#*_}
  gunzip -c $f|split -l $split_size --additional-suffix=${as%.$ext} - \
  ${f%%_*}'xxx'
  rev=${f/R1/R2}
  as='_'${rev#*_}
  echo "splitting ${f/R1/R*}..."
  gunzip -c ${f/R1/R2}|split -l $split_size --additional-suffix=${as%.$ext} - \
  ${rev%%_*}'xxx'
  if [ -d $dir/${base%%_*} ]; then
    rm -r $dir/${base%%_*}
  fi
  mkdir $dir/${base%%_*}
  gzip $dir/*.fastq
  mv $dir/*xxx* $dir/${base%%_*}
  wait;/home/micke/PycharmProjects/bashscripts/rundiam_lf.sh $dir/${base%%_*}
done



