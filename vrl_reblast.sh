#!/bin/bash

# --- Path Configuration ---
RESOURCE_DIR="/mnt/micke_ssd/resources"
export BLASTDB="$RESOURCE_DIR"
VERT_TAXIDS="$RESOURCE_DIR/vertebrate_virus_taxids.txt"
MAMM_TAXIDS="$RESOURCE_DIR/mammal_virus_taxids.txt"

usage() {
    echo "Usage: $0 -q <query.fasta> -d <db_name> [-o <prefix>] [-h <v|m>]"
    echo ""
    echo "Options:"
    echo "  -q    Query FASTA file (e.g., SampleName_R1_contigs.fasta)"
    echo "  -d    BLAST database name (in $RESOURCE_DIR)"
    echo "  -o    Optional: Manual output prefix (default: extracted from -q)"
    echo "  -h    Host restriction: 'v' (vertebrates) or 'm' (mammals)"
    echo ""
    exit 1
}

# Initialize variables
HOST_FLAG=""
TAX_FILTER=""
OUT_PREFIX=""

while getopts "q:d:o:h:" opt; do
    case $opt in
        q) QUERY="$OPTARG" ;;
        d) DB_NAME="$OPTARG" ;;
        o) OUT_PREFIX="$OPTARG" ;;
        h) HOST_FLAG="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$QUERY" || -z "$DB_NAME" ]]; then
    usage
fi

# --- 1. Automatic Prefix Extraction ---
# If -o is not provided, extract the sample name from the query filename
if [[ -z "$OUT_PREFIX" ]]; then
    QUERY_FILENAME=$(basename "$QUERY")
    # This strips everything from _R1 onwards (e.g. SampleA_R1_vrl.fasta -> SampleA)
    OUT_PREFIX="${QUERY_FILENAME%%_R1*}"

    # Fallback: if _R1 wasn't found, just strip the extension
    if [[ "$OUT_PREFIX" == "$QUERY_FILENAME" ]]; then
        OUT_PREFIX="${QUERY_FILENAME%.*}"
    fi
fi

DB_PATH="$RESOURCE_DIR/$DB_NAME"
TSV_OUT="${OUT_PREFIX}.tsv"
FASTA_OUT="${OUT_PREFIX}_top_hit.fasta"

# --- 2. Host Filter Setup ---
if [[ -n "$HOST_FLAG" ]]; then
    case "$HOST_FLAG" in
        v) TAX_FILTER="-taxidlist $VERT_TAXIDS" ;;
        m) TAX_FILTER="-taxidlist $MAMM_TAXIDS" ;;
        *) echo "Error: Invalid host option"; exit 1 ;;
    esac
fi

# --- 3. Run BLAST ---
echo "Running BLAST for Sample: $OUT_PREFIX"
blastn -query "$QUERY" \
       -db "$DB_PATH" \
       -out "$TSV_OUT" \
       $TAX_FILTER \
       -max_target_seqs 5 \
       -max_hsps 1 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand staxids sscinames" \
       -num_threads $(nproc)

# --- 4. Process Top Hits ---
if [[ ! -s "$TSV_OUT" ]]; then
    echo "No hits found for $OUT_PREFIX."
    exit 0
fi

# Grab the first line (best hit)
TOP_HIT=$(head -n 1 "$TSV_OUT")
SSEQID=$(echo "$TOP_HIT" | awk '{print $2}')
SSTRAND=$(echo "$TOP_HIT" | awk '{print $13}')

# Extract and Orient the sequence
if [[ "$SSTRAND" == "minus" ]]; then
    echo "--> Top Hit: $SSEQID (Reverse). Extracting RC..."
    blastdbcmd -db "$DB_PATH" -entry "$SSEQID" -strand minus > "$FASTA_OUT"
else
    echo "--> Top Hit: $SSEQID (Forward). Extracting..."
    blastdbcmd -db "$DB_PATH" -entry "$SSEQID" -strand plus > "$FASTA_OUT"
fi

echo "Done. Results: $TSV_OUT and $FASTA_OUT"