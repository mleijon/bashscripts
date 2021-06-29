#!/bin/bash

usearch='usearch11.0.667_i86linux64'

if [ ! -d "$1" ]; then
  echo "Please enter the full path to the directory containing the sequencing\
   data as argument."
  exit 0
else
  FILES=$1'/*R1*.gz'
fi
base=${FILES%/*}
base=$(dirname "$base")'/'$(basename "$base")
for f in $FILES; do
sample_name=${f##*/}
sample_name=${sample_name%%_*}
out_dir=$base/$sample_name
if [ -d "$out_dir" ]; then
  echo "Directory \"$base/$sample_name\" removed."
  rm -r "${base:?}/${sample_name:?}"
fi
mkdir "$out_dir"
fi_out=$out_dir/$sample_name'.fastq'${EXT}
echo "Trimming [Trimmomatic]: $sample_name"
#TRIMMOMATIC
wait;trimmomatic PE -threads 8 -phred33 -quiet -basein "$f" -baseout \
""$fi_out"" -summary "$out_dir"/'summary.log'  SLIDINGWINDOW:4:15 MINLEN:75
cat "$out_dir/$sample_name"* >> "$out_dir"/'__'"$sample_name"'.fastq'
rm "$out_dir"/"$sample_name"*
echo "Dereplicating [Usearch]: $sample_name"
#USEARCH
wait;$usearch -fastx_uniques \
"$out_dir"/'__'"$sample_name"'.fastq' \
-fastaout "$out_dir/$sample_name"'_uq.fasta' -sizeout -relabel Uniq \
-strand both &>/dev/null
rm "$out_dir"/'__'"$sample_name"'.fastq'
wait;gzip "$out_dir/$sample_name"'_uq.fasta'
diamin=$out_dir/$sample_name'_uq.fasta.gz'
echo "Blasting [Diamond]: $sample_name"
#DIAMOND
wait;diamond blastx -d /media/micke/7B9D13A66FF9DFB8/test/test -q "$diamin" \
-o "$out_dir/$sample_name"'.daa' --max-target-seqs 5 --evalue 1E-5 \
--outfmt 102 -b1 -c1 --compress 0 &>/dev/null
echo "Classifying [daa2spec.py]: $sample_name"
#DAA2SPEC.PY
wait;daa2spec.py -f "$out_dir/$sample_name"'.daa' -s -v &>/dev/null
gzip "$out_dir/$sample_name"'.daa'
done
