#!/bin/bash
# Safety: exit on error, undefined variables, or pipe failures
set -ueo pipefail

# ==============================================================================
# SPADIA PIPELINE V4 (Optimized, Modular & Environment-Aware)
# Trimming -> Assembly (SPAdes) -> Classification (Diamond)
# ==============================================================================

# --- CONDA ENVIRONMENT CHECK ---
if [[ "${CONDA_DEFAULT_ENV:-}" != "spadia" ]]; then
    echo "❌ Error: The 'spadia' conda environment is not active."
    echo "Please activate it using: conda activate spadia"
    exit 1
fi

# --- DEFAULT CONFIGURATION ---
DIAMOND_DB="/ssd2/classify/nr"
THREADS=$(grep -c 'processor' /proc/cpuinfo)
# Specific thread cap for Trimmomatic to avoid Java memory issues
TRIM_THREADS=$(( THREADS > 8 ? 8 : THREADS ))
MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
MEM_GB=$((MEM_KB / 1024 / 1024))
SPADES_MEM=$((MEM_GB * 8 / 10))  # Default: 80% of system RAM

# Default Flags
SPADES_FLAG="only-assembler" # Default to 'only-assembler' for meta-analysis
RAW_MODE="n"                 # Skip trimming if 'y'
REMOVE_TEMP="n"              # Keep intermediate files by default
DIAMOND_BLOCK=12             # Diamond block size

# --- HELP FUNCTION ---
Help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -x [y|n]    Raw reads mode (skip trimming). Default: n"
    echo "  -f [flag]   SPAdes mode (meta, plasmid, etc.). Default: only-assembler"
    echo "  -m [0.1-1]  Fraction of RAM for SPAdes. Default: 0.8"
    echo "  -p [int]    Number of threads. Default: $THREADS"
    echo "  -r [y|n]    Remove intermediate files. Default: n"
    echo "  -h          Show this help"
}

# --- ARGUMENT PARSING ---
while getopts "m:p:b:f:r:x:h" opt; do
    case ${opt} in
        m) SPADES_MEM=$(awk "BEGIN {print int($MEM_GB * $OPTARG)}") ;;
        p) THREADS=$OPTARG ;;
        f) SPADES_FLAG=$OPTARG ;;
        b) DIAMOND_BLOCK=$OPTARG ;;
        r) REMOVE_TEMP=$OPTARG ;;
        x) RAW_MODE=$OPTARG ;;
        h) Help; exit 0 ;;
        *) Help; exit 1 ;;
    esac
done

# --- VALIDATION ---
# 1. Check if the Diamond database exists
if [[ ! -f "${DIAMOND_DB}.dmnd" ]]; then
    echo "❌ Error: Diamond database file '${DIAMOND_DB}.dmnd' not found."
    echo "Please ensure the path is correct and the database is built."
    exit 1
fi

# 2. Check if the input directory exists
if [[ ! -d "./fa" ]]; then
    echo "❌ Error: Input directory './fa' not found."
    echo "Please create a folder named 'fa' and place your raw reads there."
    exit 1
fi

# 3. Check for .fastq or .fastq.gz files in ./fa
if [[ -z $(find ./fa -maxdepth 1 \( -name "*.fastq" -o -name "*.fastq.gz" \) -print -quit) ]]; then
    echo "❌ Error: No .fastq or .fastq.gz files found in './fa'."
    exit 1
fi

# --- DIRECTORY SETUP ---
TRIM_DIR="./trimmed"
ASM_DIR="./spades_out"
DMND_DIR="./diamond"

mkdir -p "$TRIM_DIR" "$DMND_DIR"

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

        trimmomatic PE -threads "$TRIM_THREADS" -phred33 -quiet \
            "$f1" "${f1/_R1_/_R2_}" \
            "${out_prefix}_1P.fastq.gz" "${out_prefix}_1U.fastq.gz" \
            "${out_prefix}_2P.fastq.gz" "${out_prefix}_2U.fastq.gz" \
            SLIDINGWINDOW:4:15 MINLEN:75
    done
fi

# ==============================================================================
# 2. ASSEMBLY (SPAdes)
# ==============================================================================
echo "🚀 Starting SPAdes Assembly..."
for f1p in "$TRIM_DIR"/*_1P.fastq.gz; do
    [ -e "$f1p" ] || continue
    sample_id=$(basename "$f1p" _1P.fastq.gz)
    out_dir="$ASM_DIR/$sample_id"

    if [[ -f "$out_dir/scaffolds.fasta" ]]; then
        echo ">> ✅ Skipping Assembly: $sample_id"
        continue
    fi

    # Using the defined SPAdes flag and resource limits
    spades.py --${SPADES_FLAG} -t "$THREADS" -m "$SPADES_MEM" \
        -1 "$f1p" -2 "${f1p/_1P/_2P}" -o "$out_dir"
done

# ==============================================================================
# 3. CLASSIFICATION (DIAMOND)
# ==============================================================================
echo "💎 Starting Diamond Classification..."
for f in "$ASM_DIR"/*/scaffolds.fasta; do
    [ -e "$f" ] || continue
    # Get parent folder name as sample name
    sample_name=$(basename "$(dirname "$f")")
    daa_out="$DMND_DIR/${sample_name}.daa"

    if [[ -f "$daa_out" || -f "${daa_out}.gz" ]]; then
        echo ">> ✅ Skipping Classification: $sample_name"
        continue
    fi

    diamond blastx -d "$DIAMOND_DB" -q "$f" -o "$daa_out" \
        --max-target-seqs 1 --evalue 1E-5 \
        --outfmt 6 qseqid full_qseq evalue staxids sscinames sskingdoms skingdoms sphylums \
        -b "$DIAMOND_BLOCK" -c 1 --compress 0

    # Optimized Single-Pass Kingdom Sorting
    echo ">> 🧬 Sorting $sample_name by Kingdom..."
    for kingdom in Viruses Eukaryota Bacteria; do
        awk -v K="$kingdom" -v FS='\t' '$6 == K { gsub(/ /,"_"); print ">"$1":"$5"\n"$2 }' "$daa_out" \
            > "$DMND_DIR/${sample_name}_${kingdom,,}.fa"
    done

    gzip "$daa_out"
done

# --- CLEANUP ---
if [[ "$REMOVE_TEMP" == "y" ]]; then
    echo "🧹 Removing intermediate files..."
    rm -rf "$TRIM_DIR" "$ASM_DIR"
fi

echo "🏁 Spadia Pipeline Complete."