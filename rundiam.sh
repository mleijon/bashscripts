#!/bin/bash
# Safety: exit on error or undefined variables
set -ueo pipefail

# --- CONFIGURATION (Fill these in or use command-line arguments) ---
# Example: DB_PATH="/path/to/diamond/nr"
DB_PATH="${1:-}"        # Takes the first argument from command line
INPUT_PATTERN="${2:-}"  # Takes the second argument (e.g., "/data/*/contigs.fasta")
OUT_DIR="${3:-}"        # Takes the third argument

# --- USAGE CHECK ---
if [[ -z "$DB_PATH" || -z "$INPUT_PATTERN" || -z "$OUT_DIR" ]]; then
    echo "Usage: $0 <DB_PATH> <INPUT_PATTERN> <OUT_DIR>"
    echo "Example: $0 /db/nr '/spades/*/contigs.fasta' ./diamond_out"
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "🚀 Starting Diamond Blastx Pipeline..."

# Use quotes around $INPUT_PATTERN to handle paths correctly
for f in $INPUT_PATTERN; do
    # Avoid loop errors if no files match
    [ -e "$f" ] || continue

    # Extract sample name (Logic adapted from your original script)
    # This assumes 'spades/' is in the path; adjust if your folder structure changes.
    outname=$(basename "$(dirname "$f")")
    output_file="$OUT_DIR/${outname}.daa"

    # --- IDEMPOTENCY CHECK ---
    # Skips files that have already been processed
    if [[ -f "$output_file" ]]; then
        echo ">> ✅ Skipping (Already Finished): $outname"
        continue
    fi

    echo "------------------------------------------------"
    echo "💎 Blasting: $outname"
    echo "------------------------------------------------"

    # Run Diamond
    diamond blastx \
        -d "$DB_PATH" \
        -q "$f" \
        -o "$output_file" \
        --max-target-seqs 1 \
        --evalue 1E-5 \
        --outfmt 6 qseqid full_qseq evalue staxids sscinames sskingdoms skingdoms sphylums \
        -b 20 -c 1 --compress 0
done

echo "------------------------------------------------"
echo "All files processed! Output is in $OUT_DIR"