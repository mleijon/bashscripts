#!/bin/bash
set -u

# --- PARAMETERS ---
# Tolerance for sequence length compared to the median sequence found (%)
TOLERANCE_PCT=70
# Sequence identity threshold for clustering (e.g. 0.90 for 90%)
IDENTITY_THR=0.90
# Minimum number of sequences required in a cluster to keep the representative sequence
MIN_CLUSTER_SIZE=2
# ------------------

# Check if a gene name was provided as an argument
if [ -z "${1:-}" ]; then
    echo "Usage: $0 <gene_name>"
    exit 1
fi

GENE_NAME="$1"

# Query definition - Host criteria removed, broad text search used to catch alternative annotations
QUERY="(${GENE_NAME}[All Fields] OR ${GENE_NAME}[Protein Name]) AND \"Staphylococcus aureus\"[Organism]"

# URL encoding the query
URL_QUERY=$(echo "$QUERY" | sed 's/ /%20/g; s/\[/%5B/g; s/\]/%5D/g; s/"/%22/g; s/(/%28/g; s/)/%29/g; s/+/ %2B/g')

echo "🔍 Searching for $GENE_NAME ID numbers..."
SEARCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=${URL_QUERY}&usehistory=y"
RESPONSE=$(wget -qO- "$SEARCH_URL")

# Use head -n 1 to handle single-line XML responses correctly
WEBEV=$(echo "$RESPONSE" | grep -oP "(?<=<WebEnv>)[^<]+" | head -n 1)
QUERYKEY=$(echo "$RESPONSE" | grep -oP "(?<=<QueryKey>)[^<]+" | head -n 1)
TOTAL=$(echo "$RESPONSE" | grep -oP "(?<=<Count>)[^<]+" | head -n 1)

if [ -z "$TOTAL" ] || [ "$TOTAL" -eq 0 ]; then
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
    sleep 0.2
done < temp_ids.txt

echo -e "\n✅ Download complete."

# ---------------------------------------------------------
# Filter truncated sequences based on TOLERANCE_PCT (90th Percentile)
# ---------------------------------------------------------
echo "📏 Calculating dynamic length threshold based on the 90th percentile..."

# Beräkna längden på sekvensen vid den 90:e percentilen (filtrerar bort gaps med -g)
P90_LEN=$(seqkit fx2tab -g -l -n temp_raw.fasta | cut -f2 | sort -n | awk '
{
    val[NR] = $1
}
END {
    # Hitta index för 90:e percentilen
    idx = int(NR * 0.90)
    if (idx == 0) idx = 1
    print val[idx]
}')

if [ -n "$P90_LEN" ] && [ "$P90_LEN" -gt 0 ]; then
    # Beräkna minsta tillåtna längd baserat på percentilen
    MIN_LEN=$(awk -v p90="$P90_LEN" -v pct="$TOLERANCE_PCT" 'BEGIN { printf "%d", p90 * (pct / 100) }')
    echo "   -> 90th percentile length: $P90_LEN bp. Filtering out sequences < $MIN_LEN bp."
    
    # Kör seqkit med -g för att hantera spaces och -m för att filtrera
    COUNT_BEFORE=$(grep -c "^>" temp_raw.fasta || true)
    seqkit seq -g -m "$MIN_LEN" temp_raw.fasta > temp_filtered.fasta
    COUNT_AFTER=$(grep -c "^>" temp_filtered.fasta || true)
    FILTERED_OUT=$(( COUNT_BEFORE - COUNT_AFTER ))
    
    echo "   -> ✂️ Removed $FILTERED_OUT truncated sequence(s)."
    mv temp_filtered.fasta temp_raw.fasta
else
    echo "⚠️ Could not determine reference length, skipping filtering."
fi
# ---------------------------------------------------------

# ---------------------------------------------------------
# Format headers, remove exact duplicates, and cluster
# ---------------------------------------------------------
echo "🧹 Formatting and removing exact duplicates for $GENE_NAME..."
awk -v gene="$GENE_NAME" '/^>/ { print ">" gene "_raw" ++i } !/^>/ { print }' temp_raw.fasta > temp_clean.fasta
seqkit rmdup -s temp_clean.fasta -o temp_unique.fasta

echo "🧬 Clustering sequences with MMseqs2 (Identity >= ${IDENTITY_THR})..."
mmseqs easy-cluster temp_unique.fasta "${GENE_NAME}_mmseqs" tmp_mmseqs --min-seq-id "$IDENTITY_THR" -c 0.8 --cov-mode 1 > /dev/null 2>&1

# Hitta rätt namn på representantfilen beroende på MMseqs2-version
if [ -f "${GENE_NAME}_mmseqs_rep.fasta" ]; then
    REP_FILE="${GENE_NAME}_mmseqs_rep.fasta"
elif [ -f "${GENE_NAME}_mmseqs_rep_seq.fasta" ]; then
    REP_FILE="${GENE_NAME}_mmseqs_rep_seq.fasta"
else
    echo "❌ Error: MMseqs2 representative file not found!"
    exit 1
fi

# ---------------------------------------------------------
# Filter out small outlier clusters to remove annotation noise
# ---------------------------------------------------------
echo "📊 Filtrerar bort outlier-kluster (krav: minst $MIN_CLUSTER_SIZE sekvenser i klustret)..."

awk '{print $1}' "${GENE_NAME}_mmseqs_cluster.tsv" | sort | uniq -c | awk -v min="$MIN_CLUSTER_SIZE" '$1 >= min {print $2}' > valid_reps.txt

TOTAL_CLUSTERS=$(cut -f1 "${GENE_NAME}_mmseqs_cluster.tsv" | sort -u | wc -l)
VALID_CLUSTERS=$(wc -l < valid_reps.txt)
REMOVED_CLUSTERS=$((TOTAL_CLUSTERS - VALID_CLUSTERS))
echo "   -> ✂️ Kastade $REMOVED_CLUSTERS små kluster."

seqkit grep -f valid_reps.txt "$REP_FILE" > temp_filtered_reps.fasta
awk -v gene="$GENE_NAME" '/^>/ { print ">" gene "_var" ++i } !/^>/ { print }' temp_filtered_reps.fasta > "${GENE_NAME}.fasta"

# Cleanup (inklusive dina faktiska filnamn från outputen)
rm -f temp_raw.fasta temp_ids.txt temp_clean.fasta temp_unique.fasta valid_reps.txt temp_filtered_reps.fasta
rm -rf tmp_mmseqs "${GENE_NAME}_mmseqs_all.fasta" "${GENE_NAME}_mmseqs_cluster.tsv" "${GENE_NAME}_mmseqs_all_seqs.fasta" "${GENE_NAME}_mmseqs_rep.fasta" "${GENE_NAME}_mmseqs_rep_seq.fasta"

hits=$(grep -c "^>" "${GENE_NAME}.fasta" || true)
echo "🎉 Finished! Databasen innehåller $hits validerade varianter sparade till: ${GENE_NAME}.fasta"
