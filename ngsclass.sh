#!/bin/bash
set -ueo pipefail

# ==============================================================================
# NGSCLASS.SH (v1.2)
# Modular NGS Classification Pipeline
# ==============================================================================

# --- CONDA ENVIRONMENT CHECK ---
if [[ "${CONDA_DEFAULT_ENV:-}" != "ngsclass" ]]; then
    echo "❌ Error: The 'ngsclass' conda environment is not active."
    echo "Please activate it using: mamba activate ngsclass"
    exit 1
fi

# --- 1. Global Configuration ---
DIAMOND_DB="/mnt/micke_ssd/resources/nr_cluster_seq_2026_270.0"
USEARCH="usearch11.0.667_i86linux64"

# Resources & Defaults
THREADS=$(grep -c 'processor' /proc/cpuinfo)
TRIM_THREADS=8
DIAMOND_BLOCK=12
DIAMOND_CHUNKS=1
RAW_MODE="n"        # -x y: Derep mode
ASM_MODE="m"        # -a m: Megahit, -a s: SPAdes
SPADES_FLAGS=""     # -f: SPAdes specific mode
REMOVE_TEMP="n"

# Versioning for Diamond directory
GB_VER="${DIAMOND_DB##*_}"
DMND_DIR="./diamond_${GB_VER}"

export TRIMMOMATIC_JAVA_OPTS="-Xmx128G -XX:ParallelGCThreads=2"

# --- 2. Help & Argument Parsing ---
Help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -x [y|n]    Dereplication mode (skip assembly). Default: $RAW_MODE"
    echo "  -a [m|s]    Assembly tool: (m)egahit or (s)pades. Default: $ASM_MODE"
    echo "  -f [str]    SPAdes mode flags. Options include:"
    echo "                --isolate       (High-coverage isolates)"
    echo "                --sc            (Single-cell MDA data)"
    echo "                --meta          (Metagenomic data)"
    echo "                --rna           (Transcriptomic data)"
    echo "                --plasmid       (Plasmid assembly from WGS)"
    echo "                --metaviral     (Viral assembly from metagenomes)"
    echo "                --metaplasmid   (Plasmid assembly from metagenomes)"
    echo "                --rnaviral      (Viral RNA-Seq data)"
    echo "                --bio           (Biosynthetic gene clusters)"
    echo "                --corona        (SARS-CoV-2 data)"
    echo "                --sewage        (Wastewater metagenomics)"
    echo "              Default: $SPADES_FLAGS"
    echo "  -p [int]    Threads. Default: $THREADS"
    echo "  -b [int]    Diamond block size. Default: $DIAMOND_BLOCK"
    echo "  -c [int]    Diamond chunks. Default: $DIAMOND_CHUNKS"
    echo "  -r [y|n]    Remove intermediate files. Default: $REMOVE_TEMP"
    echo "  -h          Show this help"
}

while getopts "x:a:f:p:b:c:r:h" opt; do
    case ${opt} in
        x) RAW_MODE=$OPTARG ;;
        a) ASM_MODE=$OPTARG ;;
        f) SPADES_FLAGS=$OPTARG ;;
        p) THREADS=$OPTARG ;;
        b) DIAMOND_BLOCK=$OPTARG ;;
        c) DIAMOND_CHUNKS=$OPTARG ;;
        r) REMOVE_TEMP=$OPTARG ;;
        h) Help; exit 0 ;;
        *) Help; exit 1 ;;
    esac
done

# --- 3. Directory Setup ---
TRIM_DIR="./trimmed"
ASM_DIR="./assembly"
DERE_DIR="./dereplicated"
mkdir -p "$ASM_DIR" "$DMND_DIR" "./logs"
if [[ "$RAW_MODE" == "n" ]]; then
    mkdir -p "$TRIM_DIR"
else
    mkdir -p "$DERE_DIR"
fi

#

# ==============================================================================
# SECTION 1: TRIMMING OR DEREPLICATION
# ==============================================================================
if [[ "$RAW_MODE" == "n" ]]; then
    echo "✂️ Mode: Trimming (Preparing for Assembly)"
    for f1 in ./fa/*_R1_*.fastq.gz; do
        [ -e "$f1" ] || continue
        sample_name=$(basename "$f1" | sed 's/_L.*//')
        out_p="$TRIM_DIR/${sample_name}"
        if [[ ! -f "${out_p}_1P.fastq.gz" ]]; then
            trimmomatic PE -threads "$TRIM_THREADS" -phred33 -quiet \
                -summary "./logs/${sample_name}"_trimmed.log "$f1" "${f1/_R1_/_R2_}" \
                "${out_p}_1P.fastq.gz" "${out_p}_1U.fastq.gz" \
                "${out_p}_2P.fastq.gz" "${out_p}_2U.fastq.gz" \
                SLIDINGWINDOW:4:15 MINLEN:75
        fi
    done
else
    echo "🌀 Mode: Raw Read Dereplication"
    for f1 in ./fa/*_R1_*.fastq.gz; do
        [ -e "$f1" ] || continue
        sample_name=$(basename "$f1" | sed 's/_L.*//')
        derep_f="$DERE_DIR/${sample_name}_derep.fa"
        if [[ ! -f "$derep_f" ]]; then
           "$USEARCH" -fastx_uniques <(zcat "$f1" "${f1/_R1_/_R2_}") \
          -fastaout "$derep_f" \
          -sizeout -relabel "${sample_name}_" \
          -threads "$THREADS"
        fi
    done
fi

# ==============================================================================
# SECTION 2: ASSEMBLY
# ==============================================================================
if [[ "$RAW_MODE" == "n" ]]; then
    for f1p in "$TRIM_DIR"/*_1P.fastq.gz; do
        [ -e "$f1p" ] || continue
        sample_id=$(basename "$f1p" _1P.fastq.gz)
        out_dir="$ASM_DIR/$sample_id"

        if [[ "$ASM_MODE" == "m" ]]; then
            if [[ ! -f "$out_dir/final.contigs.fa" ]]; then
                echo "🧬 Megahit Assembling: $sample_id"
                megahit -o "$out_dir" -1 "$f1p" -2 "${f1p/_1P/_2P}" -t "$THREADS" --continue \
                2> "./logs/$sample_id.megahit.log"
            fi
        else
            if [[ ! -f "$out_dir/contigs.fasta" ]]; then
                echo "🧬 SPAdes Assembling ($SPADES_FLAGS): $sample_id"
                spades.py $SPADES_FLAGS -t "$THREADS" -m 450 \
                    -1 "$f1p" -2 "${f1p/_1P/_2P}" -o "$out_dir" \
                    2>&1 > "./logs/$sample_id.spades.log"
            fi
        fi
    done
fi

# ==============================================================================
# SECTION 3: CLASSIFICATION (Diamond)
# ==============================================================================
echo "💎 Starting Diamond Classification..."

if [[ "$RAW_MODE" == "y" ]]; then
    INPUTS=("$DERE_DIR"/*_derep.fa)
else
    [[ "$ASM_MODE" == "m" ]] && contig_file="final.contigs.fa" || contig_file="contigs.fasta"
    INPUTS=("$ASM_DIR"/*/"$contig_file")
fi

for f in "${INPUTS[@]}"; do
    [ -e "$f" ] || continue
    sample_name=$( [[ "$RAW_MODE" == "y" ]] && basename "$f" _derep.fa || basename "$(dirname "$f")" )
    tsv_out="$DMND_DIR/${sample_name}.tsv"

    if [[ ! -f "$tsv_out" && ! -f "${tsv_out}.gz" ]]; then
        echo "🔍 Classifying: $sample_name"
        # Updated --outfmt to include sphylums for dmnd2class.sh compatibility
        diamond blastx -d "$DIAMOND_DB" -q "$f" -o "$tsv_out" \
            --max-target-seqs 1 --evalue 1E-5 -b "$DIAMOND_BLOCK" -c "$DIAMOND_CHUNKS" \
            --outfmt 6 qseqid full_qseq evalue staxids sscinames sskingdoms skingdoms sphylums \
            2> "./logs/$sample_name.diamond.log"

        for K in Viruses Eukaryota Bacteria; do
            awk -v K="$K" -v FS='\t' '$6 == K { gsub(/ /,"_"); print ">"$1":"$5"\n"$2 }' "$tsv_out" \
            > "$DMND_DIR/${sample_name}_${K,,}.fa"
        done
        gzip "$tsv_out"
    fi
done

# ==============================================================================
# SECTION 4: CLEANUP
# ==============================================================================
if [[ "$REMOVE_TEMP" == "y" ]]; then
    echo "🧹 Cleanup active..."
    [[ "$RAW_MODE" == "n" ]] && rm -rf "$TRIM_DIR" "$ASM_DIR" || rm -rf "$DERE_DIR"
fi

echo "✅ ngsclass.sh Complete."