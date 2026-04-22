#!/bin/bash
# Safety: exit on error, undefined variables, or pipe failures
set -ueo pipefail

# ==============================================================================
# DMND2CLASS.SH (Polished & Optimized)
# Stand-alone sorting of Diamond classification output into Kingdom-specific FASTA
# ==============================================================================

# --- CONFIGURATION ---
# Target directory containing Diamond output (Tabular Outfmt 6)
DMND_DIR="./diamond"

# Check if directory exists
if [[ ! -d "$DMND_DIR" ]]; then
    echo "❌ Error: Directory $DMND_DIR not found."
    exit 1
fi

echo "💎 Starting taxonomic sorting..."

# --- PROCESSING LOOP ---
# Handles Diamond TSV/DAA files (including gzipped)
for f in "$DMND_DIR"/*.tsv* "$DMND_DIR"/*.daa*; do
    [[ -e "$f" ]] || continue

    outname=$(basename "$f")
    outname=${outname%%.tsv*}
    outname=${outname%%.daa*}

    echo "🔍 Sorting: $outname"

    # Define output paths including Archaea
    OUT_VIR="${DMND_DIR}/${outname}_viruses.fa"
    OUT_EUK="${DMND_DIR}/${outname}_eukaryota.fa"
    OUT_BAC="${DMND_DIR}/${outname}_bacteria.fa"
    OUT_ARC="${DMND_DIR}/${outname}_archaea.fa"
    OUT_OTH="${DMND_DIR}/${outname}_other.fa"

    READ_CMD="cat"
    [[ "$f" == *.gz ]] && READ_CMD="zcat"

    $READ_CMD "$f" | awk -v FS='\t' \
        -v v_out="$OUT_VIR" -v e_out="$OUT_EUK" -v b_out="$OUT_BAC" -v a_out="$OUT_ARC" -v o_out="$OUT_OTH" \
        '{
            # Clean spaces in scientific name ($5) and phylum ($8)
            gsub(/ /, "_", $5);
            gsub(/ /, "_", $8);

            # Construct Header: >QueryID:ScientificName:Phylum
            header = ">" $1 ":" $5 ":" $8;
            gsub(/;:/, ";", header);

            seq = $2;

            # Sorting logic using regex match on Lineage ($6)
            if ($6 ~ /Viruses/) {
                print header "\n" seq > v_out;
            } else if ($6 ~ /Eukaryota/) {
                print header "\n" seq > e_out;
            } else if ($6 ~ /Bacteria/) {
                print header "\n" seq > b_out;
            } else if ($6 ~ /Archaea/) {
                print header "\n" seq > a_out;
            } else {
                print header "\n" seq > o_out;
            }
        }'
done

echo "✅ Sorting complete. Results located in: $DMND_DIR"
