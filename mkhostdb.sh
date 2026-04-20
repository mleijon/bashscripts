#!/bin/bash

# --- Configuration ---
RESOURCE_DIR="/mnt/micke_ssd/resources"
DB_URL="https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv"
DB_FILE="$RESOURCE_DIR/virushostdb.tsv"
VERT_LIST="$RESOURCE_DIR/vertebrate_virus_taxids.txt"
MAMM_LIST="$RESOURCE_DIR/mammal_virus_taxids.txt"

# Ensure the resource directory exists
mkdir -p "$RESOURCE_DIR"

echo "-------------------------------------------------------"
echo "Checking for Virus-Host DB updates..."

# --- Step 1: Conditional Download ---
# -N tells wget to only download if the remote file is newer than the local one
wget -N -P "$RESOURCE_DIR" "$DB_URL"

# Check if the download or check was successful
if [ $? -ne 0 ]; then
    echo "Error: Could not reach the Virus-Host DB server."
    exit 1
fi

echo "Database is up to date. Processing TaxID lists..."

# --- Step 2: Generate TaxID Lists ---
# We use awk with gsub to handle potential Windows line endings (\r)
# Column 10 is the Host Lineage
awk -F'\t' '
    { gsub(/\r/,"") }
    $10 ~ /Vertebrata/ { print $1 > "'"$VERT_LIST"'" }
    $10 ~ /Mammalia/   { print $1 > "'"$MAMM_LIST"'" }
' "$DB_FILE"

# --- Step 3: Sort and Unique (Post-processing) ---
# Awk output is redirected, now we clean them up in place
sort -u -o "$VERT_LIST" "$VERT_LIST"
sort -u -o "$MAMM_LIST" "$MAMM_LIST"

# --- Step 4: Summary ---
V_COUNT=$(wc -l < "$VERT_LIST")
M_COUNT=$(wc -l < "$MAMM_LIST")

echo "Success!"
echo "Vertebrate Virus TaxIDs: $V_COUNT saved to $VERT_LIST"
echo "Mammal Virus TaxIDs:     $M_COUNT saved to $MAMM_LIST"
echo "-------------------------------------------------------"