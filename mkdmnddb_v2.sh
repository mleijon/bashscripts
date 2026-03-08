#!/bin/bash
set -ueo pipefail

# --- CONFIGURATION ---
DB_NAME="nr_cluster_seq"
DEST_DIR="/mnt/micke_ssd/resources"
MAX_VOL=86
THREADS=72

cd "$DEST_DIR"

#echo "Starting download of $DB_NAME volumes..."
#for i in $(seq -f "%02g" 0 $MAX_VOL); do
    # Only download if file doesn't exist or is incomplete
#    wget -c "https://ftp.ncbi.nlm.nih.gov/blast/db/experimental/${DB_NAME}.${i}.tar.gz"
#done

#echo "Extracting volumes in parallel..."
# Using xargs to handle the extraction safely
#ls ${DB_NAME}.*.tar.gz | xargs -n 1 -P 16 tar -xzvf

echo "Building DIAMOND database (Taxonomy-aware)..."
# Piping from blastdbcmd saves ~500GB of disk space
# 512GB RAM allows in-memory sorting
blastdbcmd -db "$DB_NAME" -entry all | \
diamond makedb --in - -d "${DB_NAME}_2026" \
  --taxonmap prot.accession2taxid.FULL.gz \
  --taxonnodes nodes.dmp \
  --taxonnames names.dmp \
  --threads "$THREADS" \

echo "Build complete. Database is at ${DEST_DIR}/${DB_NAME}_2026.dmnd"
