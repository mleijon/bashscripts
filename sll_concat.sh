#!/bin/bash
# Safety: exit on error or undefined variables
set -ueo pipefail

# --- CONFIGURATION ---
# Leave empty to be provided when running, or fill in the paths here after :-
INPUT_DIR="${1:-}"  # Directory containing your sample folders or raw files
OUT_DIR="${2:-}"    # Where you want the merged files to go

# --- USAGE CHECK ---
if [[ -z "$INPUT_DIR" || -z "$OUT_DIR" ]]; then
    echo "Usage: $0 <INPUT_DIR> <OUT_DIR>"
    echo "Example: $0 /path/to/raw_data /path/to/merged_output"
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "🧵 Starting lane concatenation..."

# Iterate through files to identify unique samples
# This assumes the Illumina format: SampleName_S1_L001_R1_001.fastq.gz
for f in "$INPUT_DIR"/*_L001_R1_001.fastq.gz; do
    [ -e "$f" ] || continue

    # Extract the base sample name (everything before the _L001)
    # e.g., "Sample-A_S5_L001_R1_001.fastq.gz" -> "Sample-A_S5"
    base_name=$(basename "$f")
    sample_id="${base_name%_L00*}"

    echo "------------------------------------------------"
    echo "📦 Processing: $sample_id"

    # Define the specific output names you requested
    merged_r1="$OUT_DIR/${sample_id}_L001_R1_001.fastq.gz"
    merged_r2="$OUT_DIR/${sample_id}_L001_R2_001.fastq.gz"

    # Check if files already exist to avoid double-appending if script is re-run
    if [[ -f "$merged_r1" ]]; then
        echo ">> ⚠️ Skipping $sample_id (Output files already exist in $OUT_DIR)"
        continue
    fi

    # Concatenate all lanes (L001, L002, etc.) for this sample ID
    echo ">> Concatenating R1 lanes..."
    cat "$INPUT_DIR/${sample_id}"_*_R1_*.fastq.gz > "$merged_r1"

    echo ">> Concatenating R2 lanes..."
    cat "$INPUT_DIR/${sample_id}"_*_R2_*.fastq.gz > "$merged_r2"

    echo "✅ Created: $(basename "$merged_r1")"
    echo "✅ Created: $(basename "$merged_r2")"
done

echo "------------------------------------------------"
echo "Done! All samples merged into L001 format in: $OUT_DIR"