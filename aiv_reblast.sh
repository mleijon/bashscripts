#!/bin/bash
# Safety: exit on error or undefined variables
set -ueo pipefail

DB_NAME="/mnt/micke_ssd/resources/VRL_270.0"

# NCBI API Key (Ensure this is valid)
export NCBI_API_KEY="11a4df7a92bf7aef31e2ba1f4570300f6009"

if [[ -z "${1:-}" || -z "${2:-}" ]]; then
    echo "Usage: $0 <input.fasta> <output_directory>"
    exit 1
fi

RAW_INPUT="$1"
OUT_DIR="$2"
mkdir -p "$OUT_DIR"

# --- STEP 1: Create a cleaned temporary FASTA ---
# Truncates headers at the first space/semicolon for tool compatibility
INPUT_FILE=$(mktemp)
# First clean the headers, then filter for Influenza records
sed 's/[ ;].*//' "$RAW_INPUT" | grep -iA1 "influenza" | grep -v "^--" > "$INPUT_FILE"

BASE_NAME=$(basename "${RAW_INPUT%.*}")
MASTER_CSV="${OUT_DIR}/${BASE_NAME}_master_summary.csv"

blast_results=$(mktemp)
seen_list=$(mktemp)
# Cleanup temp files on exit
trap 'rm -f "$blast_results" "$seen_list" "$INPUT_FILE"' EXIT

echo "Sending batch request to NCBI (Remote)..."

# --- STEP 2: Run Batch BLAST (using cleaned input) ---
# Filtered for viruses (txid10239) and remote execution
blastn -query "$INPUT_FILE" -db "$DB_NAME" -max_target_seqs 5 -max_hsps 1 \
        -taxids 2955291,11320 \
        -outfmt "6 qseqid sacc sstrand bitscore stitle" > "$blast_results"

# --- STEP 3: Process and Sort results ---
while IFS=$'\t' read -r qseqid sacc strand score stitle; do

    # Skip if we already processed the top hit for this ID
    if grep -q "^$qseqid$" "$seen_list"; then
        continue
    fi

    if [ -n "$sacc" ]; then
        query_id="$qseqid"

        # --- SMART SEGMENT EXTRACTION ---
        # Priority 1: Check for "segment X"
        raw_seg=$(echo "$stitle" | grep -oiP "segment[:\s]*\K\w+")

        if [ -n "$raw_seg" ]; then
            seg_suffix="s${raw_seg,,}"
        # Priority 2: Map protein names to segments (Standard Influenza A)
        elif echo "$stitle" | grep -qi "PB2"; then seg_suffix="s1"
        elif echo "$stitle" | grep -qi "PB1"; then seg_suffix="s2"
        elif echo "$stitle" | grep -qi "PA";  then seg_suffix="s3"
        elif echo "$stitle" | grep -qiE "HA|hemagglutinin"; then seg_suffix="s4"
        elif echo "$stitle" | grep -qiE "NP|nucleoprotein"; then seg_suffix="s5"
        elif echo "$stitle" | grep -qiE "NA|neuraminidase"; then seg_suffix="s6"
        elif echo "$stitle" | grep -qiE "M1|M2|MP|matrix";   then seg_suffix="s7"
        elif echo "$stitle" | grep -qiE "NS1|NS2|NS";       then seg_suffix="s8"
        else
            seg_suffix="unknown"
        fi

        SEG_FASTA="${OUT_DIR}/${BASE_NAME}_${seg_suffix}.fa"
        SEG_CSV="${OUT_DIR}/${BASE_NAME}_${seg_suffix}.csv.tmp"

        # --- Reorient sequences to Plus strand and Append ---
        if [ "$strand" == "plus" ]; then
            seqkit faidx "$INPUT_FILE" "$qseqid" | seqkit replace -p ".*" -r "${query_id} | ${stitle}" | seqkit seq -w 0 >> "$SEG_FASTA"
        else
            # Reverse complement if on minus strand
            seqkit faidx "$INPUT_FILE" "$qseqid" | seqkit replace -p ".*" -r "${query_id} | ${stitle}" | seqkit seq -t DNA -r -p -w 0 >> "$SEG_FASTA"
        fi

        echo "$seg_suffix,$query_id,$sacc,$strand,$score,\"$stitle\"" >> "$SEG_CSV"
        echo "$qseqid" >> "$seen_list"
    fi

done < <(sort -t$'\t' -k1,1 -k4,4rn "$blast_results")

# --- STEP 4: Finalizing Summary and Sorting Output ---
echo "------------------------------------------------"
echo "Finalizing and Sorting files..."
echo "Segment_Group,Megahit_ID,Accession,Strand,Bitscore,Full_Title" > "$MASTER_CSV"

# 1. Process FASTA files
for fa in "${OUT_DIR}"/*.fa; do
    if [ -f "$fa" ]; then
        seqkit sort -l -r -w 0 "$fa" -o "${fa}.tmp" && mv "${fa}.tmp" "$fa"
        count=$(grep -c ">" "$fa")
        echo "✅ $(basename "$fa"): $count sequences"
    fi
done

# 2. Concatenate all .tmp CSV files and then remove them
# Using find ensures we catch all files created in Step 3
if ls "${OUT_DIR}"/*.csv.tmp >/dev/null 2>&1; then
    cat "${OUT_DIR}"/*.csv.tmp >> "$MASTER_CSV"
    rm "${OUT_DIR}"/*.csv.tmp
fi

echo "Done! Check: $MASTER_CSV"
