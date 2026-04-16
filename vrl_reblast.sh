#!/bin/bash
# Safety: exit on error or undefined variables
set -ueo pipefail

THREADS=$(grep -c 'processor' /proc/cpuinfo)

# Define the error handler
error_report() {
    local exit_code=$?
    local line_number=$1
    if [ "$exit_code" -ne 0 ] && [ "$exit_code" -ne 130 ]; then
        echo "❌ ERROR: Script failed at line $line_number with Exit Code: $exit_code"
    fi
}

# Define cleanup function
cleanup() {
    rm -f "$blast_results" "$seen_list" "$INPUT_FILE" GB_Release_Number 2>/dev/null
}

trap 'error_report $LINENO' ERR
trap cleanup EXIT

usage() {
    echo "Usage: $(basename "$0") [-r] <input.fasta> <output_directory>"
    echo "Options:"
    echo "  -r    Use remote NCBI BLAST (requires internet)"
    exit "${1:-0}"
}

REMOTE_MODE=false
while getopts "rh" opt; do
  case $opt in
    r) REMOTE_MODE=true ;;
    h) usage ;;
    \?) exit 1 ;;
  esac
done
shift $((OPTIND -1))

if [[ -z "${1:-}" || -z "${2:-}" ]]; then
    usage 1
fi

# 1. Database Discovery (VRL partition)
RESOURCES_DIR="/mnt/micke_ssd/resources"
LATEST_DB_FILE=$(ls -v "$RESOURCES_DIR"/VRL_[0-9]*.nsq "$RESOURCES_DIR"/VRL_[0-9]*.nal 2>/dev/null | tail -n 1)

if [[ -z "$LATEST_DB_FILE" ]]; then
    if [[ -f "$RESOURCES_DIR/VRL_latest.nal" ]]; then
        LATEST_DB_FILE="$RESOURCES_DIR/VRL_latest.nal"
    else
        echo "❌ ERROR: No local VRL BLAST databases found in $RESOURCES_DIR"
        exit 1
    fi
fi

FULL_BASE="${LATEST_DB_FILE%.*}"
[[ "$FULL_BASE" =~ \.[0-9]{2}$ ]] && DB_NAME="${FULL_BASE%.*}" || DB_NAME="$FULL_BASE"
LOCAL_VERSION=$(basename "$DB_NAME")
LOCAL_VERSION="${LOCAL_VERSION#*_}"

# 2. Version Check
if wget -q -O GB_Release_Number https://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number; then
    REMOTE_VERSION=$(< GB_Release_Number)
    if [[ "$LOCAL_VERSION" != "$REMOTE_VERSION" ]]; then
        echo "📢 NOTICE: A newer GenBank release ($REMOTE_VERSION) is available. Using local: $LOCAL_VERSION"
    fi
fi

# 3. Setup Files
RAW_INPUT="$1"
OUT_DIR="$2"
mkdir -p "$OUT_DIR"
BASE_NAME=$(basename "${RAW_INPUT%.*}")
MASTER_CSV="${OUT_DIR}/${BASE_NAME}_reblast_summary.csv"
FINAL_FASTA="${OUT_DIR}/${BASE_NAME}_reoriented.fa"

# Clean headers for BLAST compatibility
INPUT_FILE=$(mktemp)
sed 's/[ ;].*//' "$RAW_INPUT" > "$INPUT_FILE"

blast_results=$(mktemp)
seen_list=$(mktemp)

# 4. BLAST Execution
if [[ "$REMOTE_MODE" == "true" ]]; then
  echo "Running Remote BLAST..."
  blastn -query "$INPUT_FILE" -db nt -max_target_seqs 5 -max_hsps 1 \
         -remote -entrez_query "txid10239[Organism]" \
         -outfmt "6 qseqid sacc sstrand bitscore stitle" > "$blast_results"
else
  export BLASTDB=$(dirname "$DB_NAME")
  echo "Running Local BLAST against $DB_NAME..."
  blastn -query "$INPUT_FILE" -db "$DB_NAME" -max_target_seqs 5 -max_hsps 1 \
         -num_threads "$THREADS" \
         -outfmt "6 qseqid sacc sstrand bitscore stitle" > "$blast_results"
fi

# 5. Process Top Hits
echo "Accession,Query_ID,Strand,Bitscore,Full_Title" > "$MASTER_CSV"

while IFS=$'\t' read -r qseqid sacc strand score stitle; do
    if grep -q "^$qseqid$" "$seen_list"; then continue; fi

    if [ -n "$sacc" ]; then
        echo "$sacc,$qseqid,$strand,$score,\"$stitle\"" >> "$MASTER_CSV"
        echo "$qseqid" >> "$seen_list"

        # Reorient to plus strand and save
        if [ "$strand" == "plus" ]; then
            seqkit faidx "$INPUT_FILE" "$qseqid" | seqkit replace -p ".*" -r "${qseqid} | ${stitle}" | seqkit seq -w 0 >> "$FINAL_FASTA"
        else
            seqkit faidx "$INPUT_FILE" "$qseqid" | seqkit replace -p ".*" -r "${qseqid} | ${stitle}" | seqkit seq -t DNA -r -p -w 0 >> "$FINAL_FASTA"
        fi
    fi
done < <(sort -t$'\t' -k1,1 -k4,4rn "$blast_results")

echo "Done! Results in $OUT_DIR"