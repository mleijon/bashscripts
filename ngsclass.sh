#!/bin/bash
set -ueo pipefail

# ==============================================================================
# NGSCLASS.SH (v1.8.1 - Fixed Other Kingdom Extraction)
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

THREADS=$(grep -c 'processor' /proc/cpuinfo)
TRIM_THREADS=8
DIAMOND_BLOCK=12
DIAMOND_CHUNKS=1
RAW_MODE="n"
ASM_MODE="m"
SPADES_FLAGS=""
REMOVE_TEMP="n"

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
    echo "                --isolate, --sc, --meta, --rna, --plasmid,"
    echo "                --metaviral, --metaplasmid, --rnaviral, --bio,"
    echo "                --corona, --sewage"
    echo "              Default: $SPADES_FLAGS"
    echo "  -p [int]    Threads. Default: $THREADS"
    echo "  -b [int]    Diamond block size. Default: $DIAMOND_BLOCK"
    echo "  -c [int]    Diamond chunks. Default: $DIAMOND_CHUNKS"
    echo "  -r [y|n]    Remove intermediate files (keeps coverage results). Default: $REMOVE_TEMP"
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

# --- 3. Input Validation & Directory Setup ---
# This pattern matches _R1.fastq, _R1_001.fastq, and _R1_S1.fastq
INPUT_PATTERN="./fa/*_R1*fastq*"

# Use an array to check for matches safely with set -u
# We temporarily enable nullglob so the array is empty if no files match
shopt -s nullglob
files=($INPUT_PATTERN)
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
    echo "❌ Error: No input files found in ./fa/ matching: $INPUT_PATTERN"
    echo "Ensure your files are in the 'fa' folder and contain '_R1'."
    exit 1
fi

TRIM_DIR="./trimmed"; ASM_DIR="./assembly"; DERE_DIR="./dereplicated"; COV_DIR="./coverage"
mkdir -p "$ASM_DIR" "$DMND_DIR" "$COV_DIR" "./logs"
[[ "$RAW_MODE" == "n" ]] && mkdir -p "$TRIM_DIR" || mkdir -p "$DERE_DIR"

# ==============================================================================
# SECTION 1: TRIMMING OR DEREPLICATION
# ==============================================================================
if [[ "$RAW_MODE" == "n" ]]; then
    echo "✂️ Mode: Trimming (Preparing for Assembly)"
    for f1 in $INPUT_PATTERN; do
        [ -e "$f1" ] || continue
        # Robustly swap R1 for R2 regardless of following characters
        f2="${f1/_R1/_R2}"
        sample_name=$(basename "$f1" | sed 's/_L.*//')
        out_p="$TRIM_DIR/${sample_name}"
        if [[ ! -f "${out_p}_1P.fastq.gz" ]]; then
            echo "   Processing: $sample_name"
            trimmomatic PE -threads "$TRIM_THREADS" -phred33 -quiet \
                -summary "./logs/${sample_name}"_trimmed.log "$f1" "$f2" \
                "${out_p}_1P.fastq.gz" "${out_p}_1U.fastq.gz" \
                "${out_p}_2P.fastq.gz" "${out_p}_2U.fastq.gz" \
                SLIDINGWINDOW:4:15 MINLEN:75
        fi
    done
else
    echo "🌀 Mode: Raw Read Dereplication (USEARCH)"
    for f1 in $INPUT_PATTERN; do
        [ -e "$f1" ] || continue
        sample_name=$(basename "$f1" | sed 's/_L.*//')
        derep_f="$DERE_DIR/${sample_name}_derep.fa"
        if [[ ! -f "$derep_f" ]]; then
            echo "   Dereplicating: $sample_name"
            temp_fq="$DERE_DIR/${sample_name}_temp.fastq"
            gzip -dc -f "$f1" "${f1/_R1/_R2}" > "$temp_fq"
            "$USEARCH" -fastx_uniques "$temp_fq" -fastaout "$derep_f" \
                -sizeout -relabel "${sample_name}_" -threads "$THREADS" \
                &> "./logs/$sample_name.usearch.log"
            rm "$temp_fq"
        fi
    done
fi

# ==============================================================================
# SECTION 2: ASSEMBLY & COVERAGE
# ==============================================================================
if [[ "$RAW_MODE" == "n" ]]; then
    for f1p in "$TRIM_DIR"/*_1P.fastq.gz; do
        [ -e "$f1p" ] || continue
        sample_id=$(basename "$f1p" _1P.fastq.gz)
        f2p="${f1p/_1P.fastq.gz/_2P.fastq.gz}"

        # Define Singletons
        f1u="${f1p/_1P.fastq.gz/_1U.fastq.gz}"
        f2u="${f1p/_1P.fastq.gz/_2U.fastq.gz}"

        out_dir="$ASM_DIR/$sample_id"
        sample_cov_dir="$COV_DIR/$sample_id"

        if [[ "$ASM_MODE" == "m" ]]; then
            target_fa="$out_dir/final.contigs.fa"
            if [[ ! -f "$target_fa" ]]; then
                echo "🧬 Megahit Assembling: $sample_id"
                megahit -o "$out_dir" -1 "$f1p" -2 "$f2p" -r "$f1u,$f2u" -t "$THREADS" --continue \
                &> "./logs/$sample_id.megahit.log"
                if [[ -f "$target_fa" ]]; then
                    echo "📏 Sorting contigs: $sample_id"
                    seqkit sort --by-length --reverse "$target_fa" -o "${target_fa}.tmp"
                    mv "${target_fa}.tmp" "$target_fa"
                fi
            fi
        else
            target_fa="$out_dir/contigs.fasta"
            if [[ ! -f "$target_fa" ]]; then
                echo "🧬 SPAdes Assembling ($SPADES_FLAGS): $sample_id"
                spades.py $SPADES_FLAGS -t "$THREADS" -m 450 -1 "$f1p" -2 "$f2p" -s "$f1u" -s "$f2u" \
                -o "$out_dir" &> "./logs/$sample_id.spades.log"
            fi
        fi

        if [[ -f "$target_fa" && ! -f "$sample_cov_dir/coverage.txt" ]]; then
            echo "🎯 Calculating Coverage: $sample_id"
            mkdir -p "$sample_cov_dir"
            bowtie2-build --threads "$THREADS" "$target_fa" "$out_dir/bt2_idx" > /dev/null 2>&1
            bowtie2 -x "$out_dir/bt2_idx" -1 "$f1p" -2 "$f2p" -U "$f1u,$f2u" -p "$THREADS" \
            --very-fast-local --no-unal 2> "./logs/$sample_id.bowtie2.log" \
                | samtools sort -@ "$THREADS" -o "$sample_cov_dir/mapped.bam"
            samtools coverage "$sample_cov_dir/mapped.bam" -o "$sample_cov_dir/coverage.txt"
            samtools index "$sample_cov_dir/mapped.bam"
            rm "$out_dir"/bt2_idx.*
        fi
    done
fi

# ==============================================================================
# SECTION 3: CLASSIFICATION & ENRICHMENT
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
    cov_file="$COV_DIR/${sample_name}/coverage.txt"

    if [[ ! -f "$tsv_out" && ! -f "${tsv_out}.gz" ]]; then
        echo "🔍 Classifying: $sample_name"
        diamond blastx -d "$DIAMOND_DB" -q "$f" -o "$tsv_out" \
            --max-target-seqs 1 --evalue 1E-5 -b "$DIAMOND_BLOCK" -c "$DIAMOND_CHUNKS" \
            --outfmt 6 qseqid full_qseq evalue staxids sscinames sskingdoms skingdoms sphylums \
            &> "./logs/$sample_name.diamond.log"

        # --- ENRICH TSV WITH COVERAGE ---
        if [[ -f "$cov_file" ]]; then
            echo "➕ Adding coverage to TSV: $sample_name"
            # Map $1 (ID) to $7 (meandepth) from coverage.txt, then update $1 in Diamond TSV
            awk -v FS='\t' -v OFS='\t' '
                NR==FNR { if($1 != "#rname") cov[$1]=$7; next }
                { if($1 in cov) $1 = $1 ":cov_" cov[$1]; print }
            ' "$cov_file" "$tsv_out" > "${tsv_out}.tmp" && mv "${tsv_out}.tmp" "$tsv_out"
        fi

        # Extract Fastas (Now automatically includes :cov_XX in headers)
        for K in Viruses Eukaryota Bacteria; do
            awk -v K="$K" -v FS='\t' '$6 == K { gsub(/ /,"_"); print ">"$1":"$5"\n"$2 }' "$tsv_out" \
            > "$DMND_DIR/${sample_name}_${K,,}.fa"
        done

        # Catch-all for 'Other' (Non-standard kingdoms)
        awk -v FS='\t' '!($6 ~ /Viruses|Eukaryota|Bacteria/) { gsub(/ /,"_"); print ">"$1":"$5"\n"$2 }' "$tsv_out" \
        > "$DMND_DIR/${sample_name}_other.fa"

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