#!/bin/bash
set -u

# Check if a gene name was provided as an argument
if [ -z "${1:-}" ]; then
    echo "Usage: $0 <gene_name>"
    exit 1
fi

GENE_NAME="$1"

# Query definition
QUERY="${GENE_NAME}[Gene] AND \"Staphylococcus aureus\"[Organism] AND (\"Bos taurus\"[All Fields] OR \"Ovis aries\"[All Fields] OR bovine[All Fields] OR milk[All Fields]) NOT human[All Fields]"

# URL encoding the query
URL_QUERY=$(echo "$QUERY" | sed 's/ /%20/g; s/\[/%5B/g; s/\]/%5D/g; s/"/%22/g; s/(/%28/g; s/)/%29/g; s/+/ %2B/g')

echo "🔍 Searching for $GENE_NAME ID numbers..."
SEARCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=${URL_QUERY}&usehistory=y"
RESPONSE=$(wget -qO- "$SEARCH_URL")

WEBEV=$(echo "$RESPONSE" | grep -oPm1 "(?<=<WebEnv>)[^<]+")
QUERYKEY=$(echo "$RESPONSE" | grep -oPm1 "(?<=<QueryKey>)[^<]+")
TOTAL=$(echo "$RESPONSE" | grep -oPm1 "(?<=<Count>)[^<]+")

if [ "$TOTAL" -eq 0 ]; then
    echo "❌ No sequences found for $GENE_NAME. Skipping."
    exit 0
fi

echo "🧬 Found $TOTAL proteins. Downloading..."

# Fetch ID list
wget -qO- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&WebEnv=${WEBEV}&query_key=${QUERYKEY}&rettype=uilist&retmode=text" > temp_ids.txt

rm -f temp_raw.fasta
count=0

while read -r uid; do
    count=$((count + 1))
    echo -ne "  -> Fetching $count of $TOTAL...\r"
    wget -qO- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$uid&rettype=fasta_cds_na&retmode=text" >> temp_raw.fasta
    sleep 0.5
done < temp_ids.txt

echo -e "\n✅ Download complete."

# Format headers and remove duplicates
echo "🧹 Formatting and removing duplicates for $GENE_NAME..."
awk -v gene="$GENE_NAME" '/^>/ { print ">" gene "_var" ++i } !/^>/ { print }' temp_raw.fasta > temp_clean.fasta
seqkit rmdup -s temp_clean.fasta -o "${GENE_NAME}.fasta"

# Cleanup
rm -f temp_raw.fasta temp_ids.txt temp_clean.fasta

hits=$(grep -c "^>" "${GENE_NAME}.fasta" || true)
echo "🎉 Finished! Saved $hits unique sequences to: ${GENE_NAME}.fasta"
