#!/bin/bash
# Unified script for NCBI Clustered NR download and DIAMOND build
set -ueo pipefail
trap 'echo "❌ Skriptet avbröts på rad $LINENO"' ERR

# --- CONFIGURATION ---
DB_NAME="nr_cluster_seq"
DEST_DIR="/mnt/micke_ssd/resources"
NRCLUST_PATH="https://ftp.ncbi.nlm.nih.gov/blast/db/experimental/"
TAXONOMY_PATH="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy"
THREADS=72
CONDA_ENV="bio-db"

echo "Activating environment..."
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate "$CONDA_ENV"

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

# Helper function for smart downloads with MD5 verification
download_verify() {
    local base_url=$1
    local file_path=$2
    local base_name=${file_path##*/}

    # Remove trailing slash from base_url if present
    local clean_base_url=${base_url%/}

    echo "Fetching latest MD5 for $base_name..."
    # Hämta alltid den absolut senaste MD5-filen och skriv över den gamla lokala
    wget -q -O "${base_name}.md5" "${clean_base_url}/${file_path}.md5"

    if [ -f "$base_name" ]; then
        if md5sum --status -c "${base_name}.md5"; then
            echo "✅ $base_name is up-to-date and verified. Skipping download."
            return 0
        else
            echo "⚠️ MD5 mismatch for $base_name. The file is outdated or incomplete."
        fi
    fi

    echo "⬇️ Downloading/Resuming: $base_name"
    wget -c "${clean_base_url}/${file_path}"
    
    echo "🔍 Verifying newly downloaded $base_name..."
    if md5sum --status -c "${base_name}.md5"; then
        echo "✅ Verification successful."
    else
        echo "❌ ERROR: Checksum failed for $base_name!"
        echo "Deleting corrupted/outdated file so it can be downloaded fresh next run..."
        rm -f "$base_name"
        exit 1
    fi
}
download_verify "$TAXONOMY_PATH" "accession2taxid/prot.accession2taxid.FULL.gz"
download_verify "$TAXONOMY_PATH" "taxdump.tar.gz"

# 4. Download Clustered NR Volumes
echo "Fetching metadata to determine the number of volumes..."
METADATA_URL="${NRCLUST_PATH}${DB_NAME}-prot-metadata.json"

# Ladda ner JSON direkt i minnet, hitta raden med number-of-volumes och rensa bort allt utom siffrorna
NUM_VOLS=$(wget -qO- "$METADATA_URL" | grep '"number-of-volumes"' | tr -cd '0-9')

if [ -z "$NUM_VOLS" ]; then
    echo "❌ ERROR: Could not parse number of volumes from metadata!"
    exit 1
fi

# Eftersom NCBI börjar numrera på 000, blir maxindexet totalen minus 1
MAX_VOL=$((NUM_VOLS - 1))
echo "Metadata reports $NUM_VOLS volumes. Downloading from 000 to $(printf "%03d" $MAX_VOL)..."

echo "Starting download of $DB_NAME volumes..."
# Använd %03g för att alltid få tre siffror (000, 001, ..., 159)
for i in $(seq -f "%03g" 0 $MAX_VOL); do
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
