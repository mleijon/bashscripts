#!/bin/bash
set -ueo pipefail

# Default MMseqs2 parameters
MIN_ID=0.99
COV=0.80
MODE=1
THREADS=8

usage() {
    echo "Usage: $0 [-i min-seq-id] [-c coverage] [-m cov-mode] [-t threads]"
    echo "  -i : Min. sequence identity (0.0 to 1.0, default: 0.99)"
    echo "  -c : Min. coverage (0.0 to 1.0, default: 0.80)"
    echo "  -m : Coverage mode (0:bidirectional, 1:target, 2:query, default: 1)"
    echo "  -t : Threads (default: 8)"
    exit 1
}

while getopts "i:c:m:t:h" opt; do
    case ${opt} in
        i ) MIN_ID=$OPTARG ;;
        c ) COV=$OPTARG ;;
        m ) MODE=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        h ) usage ;;
        \? ) usage ;;
    esac
done

# 1. Identify the Diamond folder with GenBank version
DMND_DIR=$(find . -maxdepth 1 -type d -name "diamond_*" | head -n 1)

if [[ -z "$DMND_DIR" ]]; then
    echo "Error: Could not find a folder starting with 'diamond_'."
    exit 1
fi

echo "Using Diamond folder: $DMND_DIR"

# 2. Setup environment
OUT_DIR="./shared_viruses_out"
TMP_DIR="$OUT_DIR/tmp"
mkdir -p "$OUT_DIR"
CONCAT_FASTA="$OUT_DIR/merged_labeled.fasta"
> "$CONCAT_FASTA"

# 3. Label headers and merge
echo "Step 1: Processing samples..."
for file in "$DMND_DIR"/*_S*_viruses.fa; do
    if [[ -f "$file" ]]; then
        # Extract S-value (e.g., S1) from filename
        S_VAL=$(basename "$file" | grep -o "S[0-9]\+")

        # Prepend S-value to headers: >NODE_1 -> >S1_NODE_1
        sed "s/^>/>${S_VAL}_/" "$file" >> "$CONCAT_FASTA"
        echo "   + Added $S_VAL"
    fi
done

# Count how many unique samples were processed
TOTAL_SAMPLES=$(grep "^>" "$CONCAT_FASTA" | cut -d'_' -f1 | sort -u | wc -l)
echo "Total samples found: $TOTAL_SAMPLES"

# 4. Run MMseqs2
echo "Step 2: Clustering sequences..."
mmseqs easy-linclust "$CONCAT_FASTA" "$OUT_DIR/clusters" "$TMP_DIR" \
    --min-seq-id "$MIN_ID" \
    -c "$COV" \
    --cov-mode "$MODE" \
    --threads "$THREADS" \
    --verbosity 0

# 5. Filter for clusters shared across ALL samples and create mapping
echo "Step 3: Identifying contigs present in all $TOTAL_SAMPLES samples..."

# Use awk to:
# 1. Group members by Representative
# 2. Count unique S-prefixes
# 3. Output mapping and final ID list
awk -v total="$TOTAL_SAMPLES" '
    BEGIN { FS="\t"; OFS="," }
    {
        # Get S-prefix from member name (column 2)
        split($2, parts, "_");
        s_val = parts[1];

        # Store clusters
        rep_to_members[$1] = rep_to_members[$1] " " $2;
        cluster_samples[$1][s_val] = 1;
    }
    END {
        print "Cluster_Representative,Sample,Original_Contig_ID" > "'$OUT_DIR'/shared_viruses_map.csv"
        for (rep in cluster_samples) {
            count = 0;
            for (s in cluster_samples[rep]) count++;

            if (count == total) {
                # Save ID for extraction
                print rep > "'$OUT_DIR'/shared_ids.txt";

                # Split the member list we built and write to CSV mapping
                split(rep_to_members[rep], members, " ");
                for (i in members) {
                    m = members[i];
                    if (m != "") {
                        split(m, m_parts, "_");
                        orig_id = m;
                        sub(m_parts[1] "_", "", orig_id); # Remove prefix
                        print rep, m_parts[1], orig_id >> "'$OUT_DIR'/shared_viruses_map.csv"
                    }
                }
            }
        }
    }
' "$OUT_DIR/clusters_cluster.tsv"

# 6. Extract the final sequences
echo "Step 4: Extracting shared sequences to FASTA..."
awk 'FNR==NR{a[$1];next} /^>/{f=0; id=substr($1,2); if(id in a) f=1} f' \
    "$OUT_DIR/shared_ids.txt" "$CONCAT_FASTA" > "$OUT_DIR/final_shared_viruses.fasta"

SHARED_COUNT=$(wc -l < "$OUT_DIR/shared_ids.txt")

echo "------------------------------------------------"
echo "COMPLETED"
echo "Shared Contigs: $SHARED_COUNT"
echo "Mapping File:   $OUT_DIR/shared_viruses_map.csv"
echo "Final Fasta:    $OUT_DIR/final_shared_viruses.fasta"