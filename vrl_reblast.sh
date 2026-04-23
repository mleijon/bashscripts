#!/bin/bash

# --- Path Configuration ---
RESOURCE_DIR="/mnt/micke_ssd/resources"
export BLASTDB="$RESOURCE_DIR"
DB_NAME="VRL_latest"
VERT_TAXIDS="$RESOURCE_DIR/vertebrate_virus_taxids.txt"
MAMM_TAXIDS="$RESOURCE_DIR/mammal_virus_taxids.txt"
EXCLUDED_FILE="$RESOURCE_DIR/excluded.txt"

usage() {
    echo "Usage: $0 -q <query.fasta> [-o <prefix>] [-h <v|m>] [-e]"
    echo ""
    echo "Options:"
    echo "  -q    Query FASTA file"
    echo "  -o    Optional: Manual output prefix (default: extracted from -q)"
    echo "  -h    Host restriction: 'v' (vertebrates) or 'm' (mammals)"
    echo "  -e    Enable case-insensitive exclusion filter (uses $EXCLUDED_FILE)"
    echo ""
    exit 1
}

# Initialize variables
HOST_FLAG=""
TAX_FILTER=""
OUT_PREFIX=""
USE_EXCLUSION=false
MAX_TARGETS=20
SUFFIX="rbl"

while getopts "q:o:h:e" opt; do
    case $opt in
        q) QUERY="$OPTARG" ;;
        o) OUT_PREFIX="$OPTARG" ;;
        h) HOST_FLAG="$OPTARG" ;;
        e) USE_EXCLUSION=true ;;
        *) usage ;;
    esac
done

if [[ -z "$QUERY" ]]; then
    usage
fi

# --- 1. Automatic Prefix Extraction ---
if [[ -z "$OUT_PREFIX" ]]; then
    QUERY_FILENAME=$(basename "$QUERY")
    OUT_PREFIX="${QUERY_FILENAME%%_R1*}"
    if [[ "$OUT_PREFIX" == "$QUERY_FILENAME" ]]; then
        OUT_PREFIX="${QUERY_FILENAME%.*}"
    fi
fi

# --- 2. Host Filter Setup ---
if [[ -n "$HOST_FLAG" ]]; then
    case "$HOST_FLAG" in
        v) TAX_FILTER="-taxidlist $VERT_TAXIDS"; SUFFIX="vertebrate_viruses_rbl" ;;
        m) TAX_FILTER="-taxidlist $MAMM_TAXIDS"; SUFFIX="mammal_viruses_rbl" ;;
        *) usage ;;
    esac
fi

if [ "$USE_EXCLUSION" = true ]; then
    MAX_TARGETS=50
    SUFFIX="${SUFFIX}_cleaned"
fi

DB_PATH="$RESOURCE_DIR/$DB_NAME"
TSV_OUT="${OUT_PREFIX}_${SUFFIX}.tsv"
FASTA_OUT="${OUT_PREFIX}_${SUFFIX}.fa"

# --- 3. Run BLAST ---
echo "Running BLAST for $OUT_PREFIX against $DB_NAME..."
blastn -query "$QUERY" \
       -db "$DB_PATH" \
       -out "$TSV_OUT" \
       $TAX_FILTER \
       -max_target_seqs "$MAX_TARGETS" \
       -max_hsps 1 \
       -outfmt "6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand staxids sscinames stitle" \
       -num_threads $(nproc)

# --- 4. Optional Case-Insensitive Exclusion Filter ---
if [ "$USE_EXCLUSION" = true ]; then
    if [[ -f "$EXCLUDED_FILE" && -s "$EXCLUDED_FILE" ]]; then
        echo "--> Applying exclusion filter..."
        grep -v -i -F -f "$EXCLUDED_FILE" "$TSV_OUT" > "${TSV_OUT}.tmp" && mv "${TSV_OUT}.tmp" "$TSV_OUT"
    else
        echo "Warning: -e set but $EXCLUDED_FILE missing or empty."
    fi
fi

# --- 5. Process Unique Hits & Sort by Length ---
if [[ ! -s "$TSV_OUT" ]]; then
    echo "No valid hits found for $OUT_PREFIX."
    exit 0
fi

echo "Filtering for top hits and extracting sequences..."
# Create a temporary file for unsorted results
TMP_FASTA=$(mktemp)

# Sort by bitscore and take the first hit per query ID
sort -k1,1 -k12,12rn "$TSV_OUT" | awk -F$'\t' '!seen[$1]++' | while IFS=$'\t' read -r QSEQID SACC PIDENT LENGTH MISMATCH GAPOPEN QSTART QEND SSTART SEND EVALUE BITSCORE SSTRAND STAXIDS SSCINAMES STITLE; do

    if [ "$SSTRAND" == "plus" ]; then
        seqkit faidx "$QUERY" "$QSEQID" | \
            seqkit replace -p "(.*)" -r "\$1 hit:${SACC} ${STITLE}" | \
            seqkit seq -w 0 >> "$TMP_FASTA"
    else
        seqkit faidx "$QUERY" "$QSEQID" | \
            seqkit replace -p "(.*)" -r "\$1 hit:${SACC} ${STITLE}" | \
            seqkit seq -t DNA -r -p -w 0 >> "$TMP_FASTA"
    fi

done

# Now sort the collected FASTA sequences by length (longest first)
if [[ -s "$TMP_FASTA" ]]; then
    echo "Ordering output by sequence length (longest first)..."
    seqkit sort -l -r "$TMP_FASTA" -w 0 > "$FASTA_OUT"
    rm "$TMP_FASTA"
fi

TOTAL_HITS=$(grep -c "^>" "$FASTA_OUT" || true)
echo "-------------------------------------------------------"
echo "Done. Extracted $TOTAL_HITS sequences to $FASTA_OUT"
echo "Sorted by: Length (descending)"