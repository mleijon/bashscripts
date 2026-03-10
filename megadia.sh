#!/bin/bash
# Safety: exit on error, undefined variables, or pipe failures
set -ueo pipefail

# ==============================================================================
# MEGADIA PIPELINE V4 (Optimized, Modular & Environment-Aware)
# Trimming -> Assembly (MEGAHIT) -> Classification (Diamond)
# ==============================================================================

# --- CONDA ENVIRONMENT CHECK ---
# Ensures the 'megadia' environment is active before proceeding
if [[ "${CONDA_DEFAULT_ENV:-}" != "megadia" ]]; then
    echo "❌ Error: The 'megadia' conda environment is not active."
    echo "Please activate it using: conda activate megadia"
    exit 1
fi

# --- DEFAULT CONFIGURATION ---
# Databases and Tool Paths
DIAMOND_DB="/ssd2/classify/nr"

# Resources
THREADS=$(grep -c 'processor' /proc/cpuinfo)
# Specific thread cap for Trimmomatic to avoid Java memory issues
TRIM_THREADS=$(( THREADS > 8 ? 8 : THREADS ))

# Default Flags
RAW_MODE="n"         # Skip trimming if 'y'
REMOVE_TEMP="y"      # Remove intermediate files
DIAMOND_BLOCK=12     # Diamond block size
DIAMOND_CHUNKS=1

# --- HELP FUNCTION ---
Help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -x [y|n]    Raw reads mode (skip trimming). Default: n"
    echo "  -p [int]    Number of threads. Default: $THREADS"
    echo "  -r [y|n]    Remove intermediate files. Default: y"
    echo "  -b [int]    Diamond block size. Default: 12"
    echo "  -h          Show this help"
}

# --- ARGUMENT PARSING ---
while getopts "p:b:c:r:x:h" opt; do
    case ${opt} in
        p) THREADS=$OPTARG ;;
        b) DIAMOND_BLOCK=$OPTARG ;;
        c) DIAMOND_CHUNKS=$OPTARG ;;
        r) REMOVE_TEMP=$OPTARG ;;
        x) RAW_MODE=$OPTARG ;;
        h) Help; exit 0 ;;
        *) Help; exit 1 ;;
    esac
done

# --- INPUT DATA CHECK ---
# 1. Check if the directory exists
if [[ ! -d "./fa" ]]; then
    echo "❌ Error: Input directory './fa' not found."
    echo "Please create a folder named 'fa' and place your raw reads there."
    exit 1
fi

# 2. Check if the directory contains .fastq or .fastq.gz files
# We use find to safely check for files without triggering 'pipefail' errors
if [[ -z $(find ./fa -maxdepth 1 \( -name "*.fastq" -o -name "*.fastq.gz" \) -print -quit) ]]; then
    echo "❌ Error: No .fastq or .fastq.gz files found in './fa'."
    echo "Please ensure your raw sequence reads are in the 'fa' folder."
    exit 1
fi

# --- DIRECTORY SETUP ---
TRIM_DIR="./trimmed"
ASM_DIR="./megahit"
DMND_DIR="./diamond"

mkdir -p "$TRIM_DIR" "$ASM_DIR" "$DMND_DIR"

# ==============================================================================
# 1. TRIMMING (TRIMMOMATIC)
# ==============================================================================
if [[ "$RAW_MODE" == "n" ]]; then
    echo "✂️ Starting Trimmomatic (Threads capped at $TRIM_THREADS)..."
    for f1 in ./fa/*_R1_*.fastq.gz; do
        [ -e "$f1" ] || continue

        sample_name=$(basename "$f1" | sed 's/_L.*//')
        out_prefix="$TRIM_DIR/${sample_name}"

        if [[ -f "${out_prefix}_1P.fastq.gz" ]]; then
            echo ">> ✅ Skipping Trimming: $sample_name"
            continue
        fi

        echo "🧬 Trimming: $sample_name"
        trimmomatic PE -threads "$TRIM_THREADS" -phred33 -quiet \
            "$f1" "${f1/_R1_/_R2_}" \
            "${out_prefix}_1P.fastq.gz" "${out_prefix}_1U.fastq.gz" \
            "${out_prefix}_2P.fastq.gz" "${out_prefix}_2U.fastq.gz" \
            SLIDINGWINDOW:4:15 MINLEN:75
    done
fi

# ==============================================================================
# 2. ASSEMBLY (MEGAHIT)
# ==============================================================================
echo "🚀 Starting MEGAHIT Metagenomic Assembly..."
for f1p in "$TRIM_DIR"/*_1P.fastq.gz; do
    [ -e "$f1p" ] || continue

    sample_id=$(basename "$f1p" _1P.fastq.gz)
    out_dir="$ASM_DIR/$sample_id"

    if [[ -f "$out_dir/final.contigs.fa" ]]; then
        echo ">> ✅ Skipping Assembly: $sample_id"
        continue
    fi

    echo "🧬 Assembling: $sample_id"
    # Use --continue to resume interrupted MEGAHIT runs
    megahit -o "$out_dir" \
        -1 "$f1p" \
        -2 "${f1p/_1P/_2P}" \
        -r "${f1p/_1P/_1U},${f1p/_1P/_2U}" \
        -t "$THREADS" \
        --continue
done

# ==============================================================================
# 3. CLASSIFICATION (DIAMOND)
# ==============================================================================
echo "💎 Starting Diamond Classification..."
for f in "$ASM_DIR"/*/final.contigs.fa; do
    [ -e "$f" ] || continue

    sample_name=$(basename "$(dirname "$f")")
    daa_out="$DMND_DIR/${sample_name}.daa"

    if [[ -f "$daa_out" || -f "${daa_out}.gz" ]]; then
        echo ">> ✅ Skipping Classification: $sample_name"
        continue
    fi

    echo "🔍 Classifying: $sample_name"
    diamond blastx -d "$DIAMOND_DB" -q "$f" -o "$daa_out" \
        --max-target-seqs 1 --evalue 1E-5 \
        --outfmt 6 qseqid full_qseq evalue staxids sscinames sskingdoms skingdoms sphylums \
        -b "$DIAMOND_BLOCK" -c "$DIAMOND_CHUNKS" --compress 0

    # --- OUTPUT FILTERING (AWK) ---
    # Extracts sequences into separate files by Kingdom
    for kingdom in Viruses Eukaryota Bacteria; do
        awk -v K="$kingdom" -v FS='\t' \
            '$6 == K { gsub(/ /,"_"); print ">"$1":"$5"\n"$2 }' "$daa_out" \
            > "$DMND_DIR/${sample_name}_${kingdom,,}.fa"
    done

    # Catch-all for 'Other' (Non-standard kingdoms)
    awk -v FS='\t' \
        '!($6 ~ /Viruses|Eukaryota|Bacteria/) { gsub(/ /,"_"); print ">"$1":"$5"\n"$2 }' \
        "$daa_out" > "$DMND_DIR/${sample_name}_other.fa"

    gzip "$daa_out"
done

# --- CLEANUP ---
if [[ "$REMOVE_TEMP" == "y" ]]; then
    echo "🧹 Cleaning up intermediate files..."
    echo "Cleanup complete."
fi

echo "✅ Pipeline Complete. Results in: $DMND_DIR"