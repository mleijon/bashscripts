#!/bin/bash
# Safety: exit on error or undefined variables
set -ueo pipefail

# --- CONFIGURATION ---
THREADS=72              # Adjusted for your high-performance machine
TRIMMED_DIR="./trimmed"  # Where your cleaned FASTQ files live
ASM_DIR="./megahit"      # Where the assemblies will be saved

mkdir -p "$ASM_DIR"

# ==============================================================================
# 2. ASSEMBLY (MEGAHIT)
# ==============================================================================
echo "🚀 Starting MEGAHIT Metagenomic Assembly Pipeline..."

# Loop through paired-end forward reads
for f1p in "$TRIMMED_DIR"/*_1P.fastq.gz; do
    # Check if files exist to avoid loop errors
    [ -e "$f1p" ] || continue
    
    # Extract sample name (e.g., SampleA_1P.fastq.gz -> SampleA)
    sample_name=$(basename "$f1p" _1P.fastq.gz)
    asm_out="$ASM_DIR/$sample_name"

    # --- IDEMPOTENCY CHECK ---
    # If the final assembly already exists, skip it.
    if [[ -f "$asm_out/final.contigs.fa" ]]; then
        echo ">> ✅ Skipping Assembly (Already Finished): $sample_name"
        continue
    fi

    echo "------------------------------------------------"
    echo "🧬 Assembling: $sample_name"
    echo "------------------------------------------------"

    # Run MEGAHIT
    # -1 / -2: Paired-end reads
    # -r: Unpaired/orphan reads (from trimming)
    # --continue: Resume if the run was previously interrupted
    megahit -o "$asm_out" \
        -1 "$f1p" \
        -2 "${f1p/_1P/_2P}" \
        -r "${f1p/_1P/_1U},${f1p/_1P/_2U}" \
        -t "$THREADS" \
        --continue

    # --- SUCCESS CHECK ---
    if [[ -f "$asm_out/final.contigs.fa" ]]; then
        echo ">> 🏆 Assembly Success: $sample_name"
    else
        echo ">> ❌ Error: Assembly failed for $sample_name. Check logs in $asm_out"
    fi
done

echo "------------------------------------------------"
echo "All assemblies processed! Output is in $ASM_DIR"
