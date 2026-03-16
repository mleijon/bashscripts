#!/bin/bash
# Hardening: -u (error on unset vars), -e (exit on error), -o pipefail (catch pipe errors)
set -ueo pipefail

# --- CONFIGURATION ---
SAMPLES=("S1" "S2" "S3")
THREADS=16
MMSEQS_ENV="shared_viruses_env"

echo "================================================"
echo "STARTING MULTI-SAMPLE VIRAL PIPELINE"
echo "================================================"

# 1. Run ngsclass.sh (Assumes you are in the ngsclass env already)
for SAMPLE in "${SAMPLES[@]}"; do
    echo ">>> [PRIMARY] Processing Sample: $SAMPLE"

    # Replace with your actual flags for ngsclass.sh
    ./ngsclass.sh -s "$SAMPLE" -t "$THREADS"
done

echo "------------------------------------------------"
echo ">>> [COMPARATIVE] Starting Shared Virus Analysis"
echo ">>> Using Environment: $MMSEQS_ENV"
echo "------------------------------------------------"

# 2. Run the second script inside the dedicated mmseqs environment
# 'conda run' is the most robust way to do this inside a shell script
conda run -n "$MMSEQS_ENV" ./shared_viruses.sh -t "$THREADS"

echo "================================================"
echo "PIPELINE COMPLETE"
echo "================================================"