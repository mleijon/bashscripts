#!/bin/bash

FILES="/storage/micke/kc_ranch/run2/diamond/run2_class/*eukaryota.fa"
dir=$(dirname "$FILES")
rm -rf "${dir:?}"/"$1"
mkdir "$dir"/"$1"
for f in $FILES; do
  name=$(basename "$f")
  name=${name%%_*}
  grep -i -A1 "$1" "$f"|grep -v '--' - > "$dir"/"$1"/"$name"_"$1".fa
done
