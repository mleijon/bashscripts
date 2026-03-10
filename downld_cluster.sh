#!/bin/bash
set -u

# This loop will download all 86 volumes (00 to 86) and their MD5 checks
for i in $(seq -f "%02g" 0 86); do
    wget -c -nv "https://ftp.ncbi.nlm.nih.gov/blast/db/experimental/nr_cluster_seq.$i.tar.gz"
    wget -c -nv "https://ftp.ncbi.nlm.nih.gov/blast/db/experimental/nr_cluster_seq.$i.tar.gz.md5"
done
