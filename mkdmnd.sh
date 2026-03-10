#!/bin/bash
# Unified script for NCBI Clustered NR download and DIAMOND build
set -ueo pipefail

# --- CONFIGURATION ---
DB_NAME="nr_clst"
DEST_DIR="/mnt/micke_ssd/resources"
NRCLUST_PATH="https://ftp.ncbi.nlm.nih.gov/blast/db/experimental/"
TAXONOMY_PATH="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy"
MAX_VOL=86
THREADS=72
CONDA_ENV="bio-db"

echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

cd "$DEST_DIR"

# 1. Download GenBank version to construct database name
echo "Fetching GenBank version..."
wget -q -O GB_Release_Number https://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number
GENBANK_VERSION=$(< GB_Release_Number)
YEAR=$(date +%Y)
DMND_DB_NAME="${DB_NAME}_${YEAR}_${GENBANK_VERSION}"
echo "Current GenBank version is $GENBANK_VERSION"
echo "Target database name will be ${DMND_DB_NAME}.dmnd"

# 2. Clean up old files
if [ -f "${DMND_DB_NAME}.dmnd" ]; then
    echo "Found existing database. Cleaning up old resources..."
    rm -f "${DMND_DB_NAME}.dmnd" nodes.dmp names.dmp taxdump.tar.gz prot.accession2taxid.FULL.gz
fi

# 3. Download Taxonomy & Metadata
echo "Downloading Taxonomy and version metadata..."

# Helper function for downloads with MD5 verification
download_verify() {
    local base_url=$1
    local file_path=$2
    local base_name=${file_path##*/}

    # Remove trailing slash from base_url if present
    local clean_base_url=${base_url%/}

    echo "Downloading: $base_name"
    wget -c "${clean_base_url}/${file_path}"
    wget -q "${clean_base_url}/${file_path}.md5"
    echo "Verifying $base_name..."
    md5sum -c "$base_name.md5"
}

download_verify "$TAXONOMY_PATH" "accession2taxid/prot.accession2taxid.FULL.gz"
download_verify "$TAXONOMY_PATH" "taxdump.tar.gz"

# 4. Download Clustered NR Volumes
echo "Starting download of $DB_NAME volumes..."
for i in $(seq -f "%02g" 0 $MAX_VOL); do
    download_verify "$NRCLUST_PATH" "${DB_NAME}.${i}.tar.gz"
done

# 5. Extraction
echo "Extracting taxonomy..."
tar -xzvf taxdump.tar.gz nodes.dmp names.dmp

echo "Extracting database volumes in parallel..."
ls ${DB_NAME}.*.tar.gz | xargs -n 1 -P 16 tar -xzvf

# 6. Build DIAMOND DB
echo "Building DIAMOND database..."

blastdbcmd -db "$DB_NAME" -entry all | \
diamond makedb --in - -d "$DMND_DB_NAME" \
  --taxonmap prot.accession2taxid.FULL.gz \
  --taxonnodes nodes.dmp \
  --taxonnames names.dmp \
  --threads "$THREADS"

echo "Build complete. Database is at ${DEST_DIR}/${DMND_DB_NAME}.dmnd"