#!/bin/bash
# Script to modify nodes.dmp and rebuild DIAMOND DB without downloading
set -ueo pipefail

# --- CONFIGURATION (Matching your mkdmnd.sh) ---
DB_NAME="nr_clst"
DEST_DIR="/mnt/micke_ssd/resources"
THREADS=72
CONDA_ENV="bio-db"

echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

cd "$DEST_DIR"

# 1. Determine Target Database Name (Matching your mkdmnd.sh logic)
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
# We use the literal tab characters inside the sed command
echo "Normalizing ranks in nodes.dmp..."

# Change 'domain' to 'superkingdom'
sed -i 's/|\tdomain\t|/|\tsuperkingdom\t|/g' nodes.dmp

# Change 'acellular root' to 'superkingdom'
sed -i 's/|\tacellular root\t|/|\tsuperkingdom\t|/g' nodes.dmp

# 3. Step 6: Build DIAMOND DB
echo "Building DIAMOND database (Step 6)..."

blastdbcmd -db "$DB_NAME" -entry all | \
diamond makedb --in - -d "$DMND_DB_NAME" \
  --taxonmap prot.accession2taxid.FULL.gz \
  --taxonnodes nodes.dmp \
  --taxonnames names.dmp \
  --threads "$THREADS"

echo "Rebuild complete. Database is at ${DEST_DIR}/${DMND_DB_NAME}.dmnd"