#!/usr/bin/env bash

vsearch --fastq_join *R1_trunc.fastq --reverse *R2_trunc.fastq --fastqout /dev/stdout|vsearch --derep_fulllength - --output testjoin_derep.fastq --sizeout --log log.txt

