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
    echo "  -e    Enable exclusion filter (uses $EXCLUDED_FILE)"
    echo ""
    exit 1
}

# Initialize variables
HOST_FLAG=""
TAX_FILTER=""
OUT_PREFIX=""
USE_EXCLUSION=false
MAX_TARGETS=5
SUFFIX="viruses_rbl"

# Add 'e' to getopts
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

# If exclusion is enabled, increase targets and update suffix
if [ "$USE_EXCLUSION" = true ]; then
    MAX_TARGETS=50
    SUFFIX="${SUFFIX}_cleaned"
fi

DB_PATH="$RESOURCE_DIR/$DB_NAME"
TSV_OUT="${OUT_PREFIX}_${SUFFIX}.tsv"
FASTA_OUT="${OUT_PREFIX}_${SUFFIX}.fasta"

# --- 3. Run BLAST ---
echo "Running BLAST for Sample: $OUT_PREFIX using $DB_NAME"
# We include 'stitle' (col 16) regardless so the TSV is always descriptive
blastn -query "$QUERY" \
       -db "$DB_PATH" \
       -out "$TSV_OUT" \
       $TAX_FILTER \
       -max_target_seqs "$MAX_TARGETS" \
       -max_hsps 1 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand staxids sscinames stitle" \
       -num_threads $(nproc)

# --- 4. Optional Exclusion Filter ---
if [ "$USE_EXCLUSION" = true ]; then
    if [[ -f "$EXCLUDED_FILE" && -s "$EXCLUDED_FILE" ]]; then
        echo "--> Applying exclusion filter from $EXCLUDED_FILE..."
        grep -v -F -f "$EXCLUDED_FILE" "$TSV_OUT" > "${TSV_OUT}.tmp" && mv "${TSV_OUT}.tmp" "$TSV_OUT"
    else
        echo "Warning: -e set but $EXCLUDED_FILE not found or empty. Skipping filter."
    fi
fi

# --- 5. Process Top Hits ---
if [[ ! -s "$TSV_OUT" ]]; then
    echo "No valid hits found for $OUT_PREFIX."
    exit 0
fi

TOP_HIT=$(head -n 1 "$TSV_OUT")
SSEQID=$(echo "$TOP_HIT" | awk '{print $2}')
SSTRAND=$(echo "$TOP_HIT" | awk '{print $13}')
STITLE=$(echo "$TOP_HIT" | awk -F'\t' '{print $16}')

echo "--> Top Valid Hit: $SSEQID"
echo "--> Definition: $STITLE"

if [[ "$SSTRAND" == "minus" ]]; then
    echo "--> Orientation: Reverse. Extracting RC..."
    blastdbcmd -db "$DB_PATH" -entry "$SSEQID" -strand minus > "$FASTA_OUT"
else
    echo "--> Orientation: Forward. Extracting..."
    blastdbcmd -db "$DB_PATH" -entry "$SSEQID" -strand plus > "$FASTA_OUT"
fi

echo "Done. Results: $TSV_OUT and $FASTA_OUT"