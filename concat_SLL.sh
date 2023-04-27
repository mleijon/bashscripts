#!/bin/bash

dir="/media/micke/7B9D13A66FF9DFB8/SA3/inte-klara/*"
for f in $dir; do
   if [ -d "$f" ]; then
      for g in "$f"/*_R1*; do
          f_name="${g##*/}"
          f_name="${f_name#*_}"
          f_name="${f_name%%_*}"
          r1_name="$f"/"$f_name"_1.fastq.gz
          r2_name=${r1_name/_1.fastq/_2.fastq}
          cat "$g" >> "$r1_name"
          cat "${g/R1/R2}" >> "$r2_name"
      done
   fi
done
