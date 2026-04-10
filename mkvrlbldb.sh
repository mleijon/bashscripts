#!/bin/bash
set -ueo pipefail

# --- CONFIGURATION ---
echo "Fetching GenBank version..."
wget -q -O GB_Release_Number https://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number
GENBANK_VERSION=$(< GB_Release_Number)
DB_NAME="VRL_$GENBANK_VERSION"
OUT_DIR="/mnt/micke_ssd/resources"
# NCBI FTP paths
GENBANK_FTP="https://ftp.ncbi.nlm.nih.gov/ncbi-asn1/"
TAXDB_FTP="https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "------------------------------------------------"
echo "Step 1: Downloading Viral ASN.1 files (gbvrl*.aso.gz)"
echo "------------------------------------------------"
lftp -c "open $GENBANK_FTP; mirror \
     --only-newer \
     --parallel=5 \
     --no-recursion \
     --include='gbvrl.*\.aso\.gz'"

echo "------------------------------------------------"
echo "Step 2: Downloading Taxonomy Database (taxdb)"
echo "------------------------------------------------"
# -N (timestamping) only downloads if the remote file is newer than the local copy
wget -N "$TAXDB_FTP"

# Only extract if the archive is newer than the extracted database files
# (Using taxdb.btd as a marker file)
if [ ! -f "taxdb.btd" ] || [ "taxdb.tar.gz" -nt "taxdb.btd" ]; then
    echo "New taxonomy data found. Extracting..."
    tar -xvzf taxdb.tar.gz taxdb.btd taxdb.bti
else
    echo "Taxonomy files are already up to date."
fi

echo "------------------------------------------------"
echo "Step 3: Checking if database needs to be (re)built"
echo "------------------------------------------------"

FILES=$(echo gbvrl*.aso.gz)
SHOULD_BUILD=false

# Check if the database index file (.nsq) exists
if [ ! -f "$DB_NAME.nsq" ]; then
    echo "Database index not found. Initializing build..."
    SHOULD_BUILD=true
else
    # Check if any .aso.gz file is newer than the existing database index.
    # If find returns a result, it means GenBank has updated files.
    NEW_FILES=$(find . -name "gbvrl*.aso.gz" -newer "$DB_NAME.nsq")
    if [ -n "$NEW_FILES" ]; then
        echo "Updated files detected. Rebuilding database..."
        SHOULD_BUILD=true
    fi
fi

if [ "$SHOULD_BUILD" = true ]; then
    echo "Running makeblastdb..."
    # We pass the expanded $FILES list to ensure all files are processed.
    makeblastdb \
        -in "$FILES" \
        -dbtype nucl \
        -input_type asn1_bin \
        -title "Viral GenBank 270.0" \
        -out "$DB_NAME" \
        -parse_seqids
else
    echo "Database $DB_NAME is already up to date. Skipping build step."
fi

echo "------------------------------------------------"
echo "Success! Database is ready in: $OUT_DIR/$DB_NAME"
