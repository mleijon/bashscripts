#!/bin/bash
# Script to modify nodes.dmp and rebuild DIAMOND DB without downloading
set -ueo pipefail

# --- CONFIGURATION ---
DB_NAME="nr_cluster_seq"
DEST_DIR="/mnt/micke_ssd/resources"
THREADS=72
CONDA_ENV="megadia"  # Use the active megadia environment

echo "Activating environment: $CONDA_ENV..."
# Robust way to find conda/mamba and activate
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

cd "$DEST_DIR"

# 1. Determine Target Database Name
if [ -f "GB_Release_Number" ]; then
    GENBANK_VERSION=$(< GB_Release_Number)
else
    echo "Fetching GenBank version..."
    wget -q -O GB_Release_Number https://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number
    GENBANK_VERSION=$(< GB_Release_Number)
fi

YEAR=$(date +%Y)
DMND_DB_NAME="${DB_NAME}_${YEAR}_${GENBANK_VERSION}"
echo "Targeting database: ${DMND_DB_NAME}.dmnd"

# 2. Modify nodes.dmp for DIAMOND compatibility
echo "Normalizing ranks in nodes.dmp for Viruses and Domains..."
# Replace 'domain' and 'acellular root' with 'superkingdom'
# Using [[:space:]] to safely handle tabs or spaces in nodes.dmp
# Corrected: Using \t to preserve the tab-delimited format required by Diamond
sed -i 's/\t|\tdomain\t|/\t|\tsuperkingdom\t|/g' nodes.dmp
sed -i 's/\t|\tacellular root\t|/\t|\tsuperkingdom\t|/g' nodes.dmpls

# 3. Build DIAMOND DB
echo "Building DIAMOND database..."
blastdbcmd -db "$DB_NAME" -entry all | \
diamond makedb --in - -d "$DMND_DB_NAME" \
  --taxonmap prot.accession2taxid.FULL.gz \
  --taxonnodes nodes.dmp \
  --taxonnames names.dmp \
  --threads "$THREADS"

echo "Build complete. Database is at ${DEST_DIR}/${DMND_DB_NAME}.dmnd"