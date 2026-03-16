#!/bin/bash
set -ueo pipefail

# --- Configuration ---
SAMPLES=("S1" "S2" "S3" "S4" "S5")  # Add your sample IDs here
THREADS=16
GENBANK_DB="/path/to/your/nr_cluster_seq_2026_270.0"

# MMseqs2 Params
MIN_ID=1.0
COV=0.9
MODE=1

echo "================================================"
echo "STARTING MULTI-SAMPLE VIRAL PIPELINE"
echo "================================================"

# 1. Run ngsclass.sh for each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo ">>> Processing Sample: $SAMPLE"

    # Example call to your existing ngsclass script
    # Adjust flags (-i, -p, etc.) to match your actual ngsclass.sh usage
    ./ngsclass.sh -i "${SAMPLE}_R1.fastq.gz" -p "${SAMPLE}_R2.fastq.gz" -t $THREADS -d $GENBANK_DB

    if [ $? -ne 0 ]; then
        echo "Error: ngsclass.sh failed on $SAMPLE. Aborting."
        exit 1
    fi
done

echo "------------------------------------------------"
echo "All samples processed. Starting Shared Virus Analysis."
echo "------------------------------------------------"

# 2. Run the shared virus identification script
# This script will automatically find the diamond_[VERSION] folder
./shared_viruses.sh -i $MIN_ID -c $COV -m $MODE -t $THREADS

echo "================================================"
echo "PIPELINE COMPLETE"
echo "================================================"