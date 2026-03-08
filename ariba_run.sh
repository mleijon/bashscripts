#!/usr/bin/env bash
set -euo pipefail

if [[ "${CONDA_DEFAULT_ENV:-}" != "ariba" ]]; then
    echo "Error: Conda environment 'ariba' is not active." >&2
    exit 1
fi

base_dir="/ngs/mikael/S_Aureus/KA/Master_Thesis_Emily_Atkins/run2"

for d in "$base_dir"/1* "$base_dir"/Ma*; do
  if [ -d "$d" ]; then
    ariba run "$base_dir"/VF.out "${d}"/*_1.fastq.gz "${d}"/*_2.fastq.gz "${d}"/results-VF
  fi
done
