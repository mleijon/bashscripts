#!/bin/bash

rm -f data.fa
temp_file=$(mktemp)
for group in $(seq $1); do
  curl "https://data.orthodb.org/current/fasta?id=${group}at$2&species=$2" -L|sed /^{/d >> "$temp_file"
done
sed /^$/d <"$temp_file"> data.fa
