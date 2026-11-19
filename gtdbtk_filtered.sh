#!/bin/bash
set -euo pipefail

# --- Configuration ---
COMPLETENESS=50
CONTAMINATION=10
CPUS=72

CHECKM2_REPORT="checkm2_out/quality_report.tsv"
BIN_DIR="out_taxvamb/bins"
MAG_DIR="hq_mq_mags"
SAMPLE_MAP="sample_map.tsv"
OUT_DIR="gtdbtk_out"

export GTDBTK_DATA_PATH="/mnt/micke_ssd/resources/gtdb_r232/release232"

# --- Pre-flight checks ---
if [[ ! -s "$CHECKM2_REPORT" ]]; then
    echo "Error: CheckM2 report not found or empty: ${CHECKM2_REPORT}" >&2
    exit 1
fi

if [[ ! -d "$BIN_DIR" ]]; then
    echo "Error: TaxVAMB bin directory not found: ${BIN_DIR}" >&2
    exit 1
fi

if [[ ! -s "$SAMPLE_MAP" ]]; then
    echo "Error: Sample map not found or empty: ${SAMPLE_MAP}" >&2
    exit 1
fi

# Resolve absolute path once so symlinks stay valid regardless of working directory
ABS_BIN_DIR=$(realpath "$BIN_DIR")

# Recreate target directory to remove stale symlinks from earlier runs
rm -rf "$MAG_DIR"
mkdir -p "$MAG_DIR"

# --- Filter CheckM2 report for HQ/MQ MAGs and create symlinks ---
# MIMAG criteria: completeness >= 50% (col 2) and contamination < 10% (col 3)
awk -F'\t' -v comp="$COMPLETENESS" -v cont="$CONTAMINATION" \
    'NR > 1 && ($2 + 0) >= comp && ($3 + 0) < cont {print $1}' \
    "$CHECKM2_REPORT" | while IFS= read -r bin; do
    src="${ABS_BIN_DIR}/${bin}.fna"
    # ln -s silently creates broken links if the target is missing, so verify explicitly
    if [[ ! -f "$src" ]]; then
        echo "Error: Expected bin file not found: ${src}" >&2
        exit 1
    fi
    ln -s "$src" "${MAG_DIR}/"
done

count=$(find "$MAG_DIR" -maxdepth 1 -type l -name '*.fna' | wc -l)
echo "Linked ${count} HQ/MQ MAGs into ${MAG_DIR}/"

if [[ "$count" -eq 0 ]]; then
    echo "Error: No MAGs met the quality thresholds; aborting GTDB-Tk." >&2
    exit 1
fi

# --- Taxonomic classification with GTDB-Tk ---
gtdbtk classify_wf \
    --genome_dir "$MAG_DIR" \
    --out_dir "$OUT_DIR" \
    --extension fna \
    --cpus "$CPUS"

# --- Annotate GTDB-Tk summary tables with original sample names ---
# Replaces the top-level symlinks in OUT_DIR while leaving the unmodified
# originals under OUT_DIR/classify/ intact.
for summary in "${OUT_DIR}/gtdbtk.ar53.summary.tsv" "${OUT_DIR}/gtdbtk.bac120.summary.tsv"; do
    if [[ ! -f "$summary" ]]; then
        echo "Skipping ${summary} (file not found)."
        continue
    fi

    awk -F'\t' -v OFS='\t' '
        # First pass: load sample_map.tsv (S-ID -> sample name) into lookup array
        NR == FNR {
            sub(/\r$/, "", $2)
            map[$1] = $2
            next
        }
        # Second pass header: prepend "sample" unless already annotated
        FNR == 1 {
            if ($1 == "sample") {
                already_added = 1
                print
            } else {
                print "sample", $0
            }
            next
        }
        already_added {
            print
            next
        }
        # Extract leading sample ID (e.g. S36 from S36C90) and prepend mapped sample name
        {
            sid = ""
            if (match($1, /^S[0-9]+/)) {
                sid = substr($1, RSTART, RLENGTH)
            }
            sample = (sid in map) ? map[sid] : "UNKNOWN"
            print sample, $0
        }
    ' "$SAMPLE_MAP" "$summary" > "${summary}.tmp" && mv "${summary}.tmp" "$summary"

    echo "Added sample column to ${summary}"
done
