#!/usr/bin/env bash
dirs1="/ngs/mikael/S_Aureus/KA/Master_Thesis_Emily_Atkins/run2/1*"
dirs2="/ngs/mikael/S_Aureus/KA/Master_Thesis_Emily_Atkins/run2/Ma*"
for d in  $dirs1 ; do
#wait; ariba run /home/mikael/ariba/card "${d}"/*_1.fastq.gz "${d}"/*_2.fastq.gz "${d}"/card
wait; ariba run /ngs/mikael/S_Aureus/KA/Master_Thesis_Emily_Atkins/run2/VF.out "${d}"/*_1.fastq.gz "${d}"/*_2.fastq.gz "${d}"/results-VF
done
for d in  $dirs2 ; do
#wait; ariba run /home/mikael/ariba/card "${d}"/*_1.fastq.gz "${d}"/*_2.fastq.gz "${d}"/card
wait; ariba run /ngs/mikael/S_Aureus/KA/Master_Thesis_Emily_Atkins/run2/VF.out "${d}"/*_1.fastq.gz "${d}"/*_2.fastq.gz "${d}"/results-VF
done

