#!/bin/bash

# --- Path Configuration ---
RESOURCE_DIR="/mnt/micke_ssd/resources"
VERT_TAXIDS="$RESOURCE_DIR/vertebrate_virus_taxids.txt"
MAMM_TAXIDS="$RESOURCE_DIR/mammal_virus_taxids.txt"

usage() {
    echo "Usage: $0 -q <query.fasta> -d <db_name> -o <output_file> [-h <v|m>]"
    echo ""
    echo "Options:"
    echo "  -q    Query FASTA file"
    echo "  -d    BLAST database name (in $RESOURCE_DIR)"
    echo "  -o    Output prefix (for .tsv and .fasta)"
    echo "  -h    Host restriction: 'v' (vertebrates) or 'm' (mammals)"
    echo ""
    exit 1
}

# Initialize variables
HOST_FLAG=""
TAX_FILTER=""

while getopts "q:d:o:h:" opt; do
    case $opt in
        q) QUERY="$OPTARG" ;;
        d) DB_NAME="$OPTARG" ;;
        o) OUT_PREFIX="$OPTARG" ;;
        h) HOST_FLAG="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$QUERY" || -z "$DB_NAME" || -z "$OUT_PREFIX" ]]; then
    usage
fi

DB_PATH="$RESOURCE_DIR/$DB_NAME"
TSV_OUT="${OUT_PREFIX}.tsv"
FASTA_OUT="${OUT_PREFIX}_top_hit.fasta"

# --- 1. Host Filter Setup ---
if [[ -n "$HOST_FLAG" ]]; then
    case "$HOST_FLAG" in
        v) TAX_FILTER="-taxidlist $VERT_TAXIDS" ;;
        m) TAX_FILTER="-taxidlist $MAMM_TAXIDS" ;;
        *) echo "Error: Invalid host option"; exit 1 ;;
    esac
fi

# --- 2. Run BLAST ---
# Using -max_target_seqs 5 and -max_hsps 1 as requested
# We include 'sstrand' in the output to check orientation
echo "Running BLAST search..."
blastn -query "$QUERY" \
       -db "$DB_PATH" \
       -out "$TSV_OUT" \
       $TAX_FILTER \
       -max_target_seqs 5 \
       -max_hsps 1 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand staxids sscinames" \
       -num_threads $(nproc)

# --- 5. Process Top Hits ---
if [[ ! -s "$TSV_OUT" ]]; then
    echo "No hits found."
    exit 0
fi

echo "Processing top hit and checking orientation..."

# Grab the first line (best hit)
TOP_HIT=$(head -n 1 "$TSV_OUT")

# Extract details from columns
SSEQID=$(echo "$TOP_HIT" | awk '{print $2}')
SSTRAND=$(echo "$TOP_HIT" | awk '{print $13}')
SSTART=$(echo "$TOP_HIT" | awk '{print $9}')
SEND=$(echo "$TOP_HIT" | awk '{print $10}')

echo "Top Hit: $SSEQID ($SSTRAND orientation)"

# Extract the sequence using blastdbcmd
# If sstrand is 'minus', we use -strand minus to get the reverse complement
if [[ "$SSTRAND" == "minus" ]]; then
    echo "--> Orientation is reverse. Extracting Reverse Complement..."
    blastdbcmd -db "$DB_PATH" -entry "$SSEQID" -strand minus > "$FASTA_OUT"
else
    echo "--> Orientation is forward. Extracting as-is..."
    blastdbcmd -db "$DB_PATH" -entry "$SSEQID" -strand plus > "$FASTA_OUT"
fi

echo "-------------------------------------------------------"
echo "Process Complete."
echo "Full table: $TSV_OUT"
echo "Top hit FASTA: $FASTA_OUT"