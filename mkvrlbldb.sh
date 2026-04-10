#!/bin/bash
set -ueo pipefail

# --- CONFIGURATION ---
DB_NAME="VRL_270.0"
OUT_DIR="/mnt/micke_ssd/resources"
# NCBI FTP paths
GENBANK_FTP="https://ftp.ncbi.nlm.nih.gov/genbank"
TAXDB_FTP="https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "------------------------------------------------"
echo "Step 1: Downloading Viral ASN.1 files (gbvrl*.aso.gz)"
echo "------------------------------------------------"
# Using --no-clobber (-nc) to avoid re-downloading existing files
wget -nc "${GENBANK_FTP}/gbvrl"*.aso.gz

echo "------------------------------------------------"
echo "Step 2: Downloading Taxonomy Database (taxdb)"
echo "------------------------------------------------"
# taxdb is required for -taxids to work in local BLAST
wget -qO- "$TAXDB_FTP" | tar -xvz

echo "------------------------------------------------"
echo "Step 3: Creating BLAST Database"
echo "------------------------------------------------"
# We quote the wildcard to let makeblastdb handle decompression internally
makeblastdb \
    -in "gbvrl*.aso.gz" \
    -dbtype nucl \
    -input_type asn1_bin \
    -title "Viral GenBank 270.0" \
    -out "$DB_NAME" \
    -parse_seqids

echo "------------------------------------------------"
echo "Success! Database created: $OUT_DIR/$DB_NAME"
