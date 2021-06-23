#!/bin/bash

FILES="/media/micke/7B9D13A66FF9DFB8/parvo/*R1*.gz"
base=${FILES%/*}
for f in $FILES; do
sample_name=${f##*/}
sample_name=${sample_name%%_*}
out_dir=$base/$sample_name
if [ -d "$out_dir" ]
then
  echo "Directory \"$base/$sample_name\" exists."
  exit
else
  mkdir "$out_dir"
  fi_fp=$out_dir/$sample_name'_fp.fastq.gz'
  fi_rp=$out_dir/$sample_name'_rp.fastq.gz'
  fi_fu=$out_dir/$sample_name'_fu.fastq.gz'
  fi_ru=$out_dir/$sample_name'_ru.fastq.gz'
fi
echo "Processing: $sample_name"
wait;trimmomatic PE -threads 8 -phred33 -quiet "$f" "${f/R1/R2}" $fi_fp $fi_fu \
$fi_rp $fi_ru SLIDINGWINDOW:4:15 MINLEN:36
done