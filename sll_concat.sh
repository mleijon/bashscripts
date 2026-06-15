#!/bin/bash
# Safety: exit on error or undefined variables
set -ueo pipefail

# --- CONFIGURATION ---
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

# Track processed samples to avoid running the same sample multiple times
declare -A processed_samples

# Iterate through files using any lane (L001, L002, etc.) to discover samples
for f in "$INPUT_DIR"/*_L00*_R1_001.fastq.gz; do
    [ -e "$f" ] || continue

    # Extract the base sample name (everything before the _L00x)
    base_name=$(basename "$f")
    # Using regex matching to clip off at _L00
    sample_id="${base_name%_L00*}"

    # Skip if we already processed this sample via another lane file
    if [[ -n "${processed_samples[$sample_id]:-}" ]]; then
        continue
    fi
    processed_samples[$sample_id]=1

    echo "------------------------------------------------"
    echo "📦 Processing: $sample_id"

    # Define the output names (hardcoded to L001 as per your original logic)
    merged_r1="$OUT_DIR/${sample_id}_L001_R1_001.fastq.gz"
    merged_r2="$OUT_DIR/${sample_id}_L001_R2_001.fastq.gz"

    # Check if files already exist
    if [[ -f "$merged_r1" ]]; then
        echo ">> ⚠️ Skipping $sample_id (Output files already exist in $OUT_DIR)"
        continue
    fi

    # Concatenate all lanes, explicitly sorting them so L002 comes before L003, etc.
    echo ">> Concatenating R1 lanes..."
    # Using 'ls | sort' inside an array expansion ensures the lanes are merged in numerical order
    cat $(ls "$INPUT_DIR/${sample_id}"_*_R1_*.fastq.gz | sort) > "$merged_r1"

    echo ">> Concatenating R2 lanes..."
    cat $(ls "$INPUT_DIR/${sample_id}"_*_R2_*.fastq.gz | sort) > "$merged_r2"

    echo "✅ Created: $(basename "$merged_r1")"
    echo "✅ Created: $(basename "$merged_r2")"
done
