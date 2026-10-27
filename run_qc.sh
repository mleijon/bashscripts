#!/bin/bash

# Exit immediately if any command fails
set -e

# Assign input arguments to variables
INPUT_DIR="$1"
OUTPUT_DIR="$2"
# Optional 3rd argument for threads (defaults to 16 if not provided)
THREADS="${3:-16}" 

# Check if the minimum required arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory> [threads]"
    echo "Example: $0 /path/to/raw_reads /path/to/qc_results 32"
    exit 1
fi

# Create the output directory if it doesn't already exist
mkdir -p "$OUTPUT_DIR"

echo "==> Starting FastQC on $INPUT_DIR using $THREADS threads..."
# Run FastQC. The -t flag determines how many files are processed simultaneously.
# Note: Adjust *.fastq.gz to *.fq.gz or *.fastq if your file extensions differ.
fastqc -t "$THREADS" -o "$OUTPUT_DIR" "$INPUT_DIR"/*.fastq.gz

echo "==> Starting MultiQC..."
# MultiQC scans the output directory for FastQC zip files and compiles the report
multiqc -o "$OUTPUT_DIR" "$OUTPUT_DIR"

echo "==> QC complete! You can view the aggregate report at: $OUTPUT_DIR/multiqc_report.html"
