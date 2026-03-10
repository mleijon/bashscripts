#!/bin/bash
# Safety: exit on error, undefined variables, or pipe failures
set -ueo pipefail

# ==============================================================================
# SPADIA PIPELINE V3 (Optimized & Modular)
# Trimming (Trimmomatic) -> Assembly (SPAdes/USEARCH) -> Classification (Diamond)
# ==============================================================================

# --- DEFAULT CONFIGURATION ---
# Databases and Tool Paths
DIAMOND_DB="/ssd2/classify/nr"
USEARCH_BIN="/ssd2/classify/usearch11.0.667_i86linux64"

# Resources
THREADS=$(grep -c 'processor' /proc/cpuinfo)
# Specific thread cap for Trimmomatic to avoid Java memory issues
TRIM_THREADS=$THREADS
if [ "$TRIM_THREADS" -gt 8 ]; then TRIM_THREADS=8; fi
MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
MEM_GB=$((MEM_KB / 1024 / 1024))
SPADES_MEM=$((MEM_GB * 8 / 10))  # Default: 80% of system RAM

# Default Flags
SPADES_FLAG=""       # Standard mode
SEQ_TYPE="c"         # 'c' for contigs, 'r' for reads
RAW_MODE="n"         # Skip trimming if 'y'
REMOVE_TEMP="y"      # Remove intermediate files
DIAMOND_BLOCK=12     # Diamond block size
DIAMOND_CHUNKS=1

# --- HELP FUNCTION ---
Help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -s [c|r]    Sequence type: contigs (c) or reads (r). Default: c"
    echo "  -x [y|n]    Raw reads mode (skip trimming/de-rep). Default: n"
    echo "  -f [flag]   SPAdes mode (meta, plasmid, etc.). Default: standard"
    echo "  -m [0.1-1]  Fraction of RAM for SPAdes. Default: 0.8"
    echo "  -p [int]    Number of threads. Default: $THREADS"
    echo "  -r [y|n]    Remove intermediate files. Default: y"
    echo "  -h          Show this help"
}

# --- ARGUMENT PARSING ---
while getopts "m:p:b:c:f:s:r:x:h" opt; do
    case ${opt} in
        m) SPADES_MEM=$(awk "BEGIN {print int($MEM_GB * $OPTARG)}") ;;
        p) THREADS=$OPTARG ;;
        f) SPADES_FLAG=$OPTARG ;;
        b) DIAMOND_BLOCK=$OPTARG ;;
        c) DIAMOND_CHUNKS=$OPTARG ;;
        s) SEQ_TYPE=$OPTARG ;;
        r) REMOVE_TEMP=$OPTARG ;;
        x) RAW_MODE=$OPTARG ;;
        h) Help; exit 0 ;;
        *) Help; exit 1 ;;
    esac
done

# --- DIRECTORY SETUP ---
TRIM_DIR="./trimmed"
ASM_DIR="./${SPADES_FLAG:-standard}_spades"
MERGED_DIR="./merged"
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

        echo "🧬 Trimming: $sample_name"
        trimmomatic PE -threads "$TRIM_THREADS" -phred33 -quiet \
            "$f1" "${f1/_R1_/_R2_}" \
            "${out_prefix}_1P.fastq.gz" "${out_prefix}_1U.fastq.gz" \
            "${out_prefix}_2P.fastq.gz" "${out_prefix}_2U.fastq.gz" \
            SLIDINGWINDOW:4:15 MINLEN:75
    done
fi

# ==============================================================================
# 2. ASSEMBLY OR DEREPLICATION
# ==============================================================================
if [[ "$SEQ_TYPE" == "c" ]]; then
    # --- SPAdes ASSEMBLY ---
    echo "🚀 Starting SPAdes Assembly..."
    for f1p in "$TRIM_DIR"/*_1P.fastq.gz; do
        [ -e "$f1p" ] || continue

        sample_id=$(basename "$f1p" _1P.fastq.gz)
        out_dir="$ASM_DIR/$sample_id"

        if [[ -f "$out_dir/contigs.fasta" ]]; then
            echo ">> ✅ Skipping Assembly: $sample_id"
            continue
        fi

        echo "🧬 Assembling: $sample_id (Mode: ${SPADES_FLAG:-standard})"

        # Build SPAdes command dynamically
        cmd="spades.py -o $out_dir -t $THREADS -m $SPADES_MEM"
        [[ -n "$SPADES_FLAG" ]] && cmd="$cmd --$SPADES_FLAG"

        if [[ "$SPADES_FLAG" == "meta" ]]; then
            $cmd -1 "$f1p" -2 "${f1p/_1P/_2P}"
        else
            $cmd -1 "$f1p" -2 "${f1p/_1P/_2P}" -s "${f1p/_1P/_1U}" -s "${f1p/_1P/_2U}"
        fi
    done
    FILES_TO_CLASS="$ASM_DIR/*/contigs.fasta"

else
    # --- USEARCH DEREPLICATION (READS MODE) ---
    echo "🧬 Processing Reads (De-replication)..."
    mkdir -p "$MERGED_DIR"
    for f1p in "$TRIM_DIR"/*_1P.fastq.gz; do
        [ -e "$f1p" ] || continue

        sample_id=$(basename "$f1p" _1P.fastq.gz)
        merged_fq="$MERGED_DIR/${sample_id}.fastq"
        uq_fasta="$MERGED_DIR/${sample_id}_uq.fasta"

        if [[ -f "${uq_fasta}.gz" ]]; then
            echo ">> ✅ Skipping De-replication: $sample_id"
            continue
        fi

        # Merge all trimmed reads (Paired and Unpaired)
        zcat "$f1p" "${f1p/_1P/_2P}" "${f1p/_1P/_1U}" "${f1p/_1P/_2U}" > "$merged_fq"

        "$USEARCH_BIN" -fastx_uniques "$merged_fq" -fastaout "$uq_fasta" \
            -sizeout -relabel Uniq -strand both

        gzip "$uq_fasta"
        rm "$merged_fq"
    done
    FILES_TO_CLASS="$MERGED_DIR/*.fasta.gz"
fi

# ==============================================================================
# 3. CLASSIFICATION (DIAMOND)
# ==============================================================================
echo "💎 Starting Diamond Classification..."
for f in $FILES_TO_CLASS; do
    [ -e "$f" ] || continue

    # Determine output name based on path
    if [[ "$SEQ_TYPE" == "c" ]]; then
        sample_name=$(basename "$(dirname "$f")")
    else
        sample_name=$(basename "$f" .fasta.gz)
    fi

    daa_out="$DMND_DIR/${sample_name}.daa"

    if [[ -f "$daa_out" ]]; then
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
    rm -rf "$TRIM_DIR" "$MERGED_DIR" "$ASM_DIR" 2>/dev/null
fi

echo "✅ Pipeline Complete. Results in: $DMND_DIR"