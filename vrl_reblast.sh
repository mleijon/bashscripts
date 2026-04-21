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
MAX_TARGETS=5
SUFFIX="viruses_rbl"

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
FASTA_OUT="${OUT_PREFIX}_${SUFFIX}.fasta"

# --- 3. Run BLAST ---
echo "Running BLAST for $OUT_PREFIX against $DB_NAME..."
blastn -query "$QUERY" \
       -db "$DB_PATH" \
       -out "$TSV_OUT" \
       $TAX_FILTER \
       -max_target_seqs "$MAX_TARGETS" \
       -max_hsps 1 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand staxids sscinames stitle" \
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

# --- 5. Process ALL Top Hits ---
if [[ ! -s "$TSV_OUT" ]]; then
    echo "No valid hits found for $OUT_PREFIX."
    exit 0
fi

echo "Processing top hits and removing line breaks with seqkit..."
# Clear/initialize output file
> "$FASTA_OUT"

# Select only the first hit for each unique query ID
awk '!seen[$1]++' "$TSV_OUT" | while read -r line; do
    QSEQID=$(echo "$line" | awk '{print $1}')
    SSEQID=$(echo "$line" | awk '{print $2}')
    SSTRAND=$(echo "$line" | awk '{print $13}')

    # Determine orientation flag
    STRAND_VAL="plus"
    [[ "$SSTRAND" == "minus" ]] && STRAND_VAL="minus"

    # Extract, prepend query ID to header, and remove line breaks with seqkit
    blastdbcmd -db "$DB_PATH" -entry "$SSEQID" -strand "$STRAND_VAL" | \
        sed "s/^>/>${QSEQID}_hit_/" | \
        seqkit seq -w 0 >> "$FASTA_OUT"
done

TOTAL_HITS=$(grep -c "^>" "$FASTA_OUT" || true)
echo "-------------------------------------------------------"
echo "Done. Extracted $TOTAL_HITS sequences to $FASTA_OUT"
echo "Format: Single-line sequences (ready for grep -A1)"