#!/bin/bash
set -euo pipefail

# ==============================================================================
# Pipeline: QC/Trimming (fastp) -> Assembly (MEGAHIT) -> Binning (TaxVAMB) -> QC (CheckM2)
# Optimized for 2x150 bp NovaSeq Paired-End Metagenomes
# ==============================================================================

# === Configuration ===
MAX_THREADS=$(grep -c 'processor' /proc/cpuinfo)
THREADS=$(( MAX_THREADS * 2 / 3 ))
FASTP_THREADS=$(( THREADS > 16 ? 16 : THREADS )) # fastp caps at 16 worker threads
MEM_BYTES=$(( 500 * 1024 * 1024 * 1024 ))        # 500 GB in bytes for MEGAHIT
MIN_CONTIG_SIZE=2000

# Directories and Files
FA_DIR="./fa"
TRIM_DIR="./trimmed_fa"
FASTP_REPORT_DIR="./fastp_reports"
ASSEMBLY_BASE="./assemblies"
export VAMB_IN="./contigs.fna.gz"
AEMB_DIR="./aemb"
ABUNDANCE_TSV="./abundance.tsv"
TAX_DIR="./taxonomy"
DIAMOND_OUT="${TAX_DIR}/diamond_out.tsv"
VAMB_TAX="${TAX_DIR}/taxonomy.tsv"
export TAXVAMB_OUT="./out_taxvamb"
export BINS_DIR="./bins_fasta"
CHECKM2_OUT="./checkm2_out"

# Databases
DIAMOND_DB="/mnt/micke_ssd/resources/nr_cluster_seq_2026_274.0.dmnd"
CHECKM2_DB="/mnt/micke_ssd/resources/CheckM2_database/uniref100.KO.1.dmnd"
export TAXONKIT_DB="/mnt/micke_ssd/resources"

echo "Activating Conda environment: vamb5..."
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate vamb5

# Pre-flight check for required binaries in vamb5
for cmd in fastp megahit concatenate.py strobealign merge_aemb.py diamond taxonkit vamb; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: Required command '$cmd' not found in environment 'vamb5'."
        echo "Install missing tools (e.g., 'mamba install -n vamb5 -c bioconda fastp megahit') and retry."
        exit 1
    fi
done

mkdir -p "$TRIM_DIR" "$FASTP_REPORT_DIR" "$ASSEMBLY_BASE" "$AEMB_DIR" "$TAX_DIR"

# === Step 1: Read Trimming (fastp) & Assembly (MEGAHIT) ===
echo "--- Step 1: Running fastp and MEGAHIT for samples in $FA_DIR ---"
CONTIG_FILES=()

for r1 in "$FA_DIR"/*_R1*.fastq.gz; do
    [ -e "$r1" ] || continue
    seqfile=$(basename "$r1")
    sample=${seqfile%%_*}
    ext=${seqfile#*"$sample"}
    r2="${FA_DIR}/${sample}${ext/_R1/_R2}"

    if [ ! -f "$r2" ]; then
        echo "Warning: Missing R2 file for $sample, skipping..."
        continue
    fi

    trim_r1="${TRIM_DIR}/${sample}_R1.trimmed.fastq.gz"
    trim_r2="${TRIM_DIR}/${sample}_R2.trimmed.fastq.gz"

    # 1a. Quality trimming, NovaSeq poly-G removal, and overlap error correction
    if [ ! -f "$trim_r1" ] || [ ! -f "$trim_r2" ]; then
        echo "Trimming and error-correcting reads for $sample with fastp..."
        fastp \
            -i "$r1" \
            -I "$r2" \
            -o "$trim_r1" \
            -O "$trim_r2" \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --cut_right \
            --cut_right_window_size 4 \
            --cut_right_mean_quality 20 \
            --correction \
            --length_required 50 \
            --thread "$FASTP_THREADS" \
            --html "${FASTP_REPORT_DIR}/${sample}_fastp.html" \
            --json "${FASTP_REPORT_DIR}/${sample}_fastp.json"
    else
        echo "Found trimmed reads for $sample, skipping fastp."
    fi

    # 1b. Metagenomic assembly with MEGAHIT
    out_dir="${ASSEMBLY_BASE}/${sample}"
    contig_fasta="${out_dir}/final.contigs.fa"

    if [ ! -f "$contig_fasta" ]; then
        echo "Running MEGAHIT assembly for $sample..."
        # MEGAHIT requires the target directory not to exist prior to launch
        if [ -d "$out_dir" ]; then
            rm -rf "$out_dir"
        fi

        # --k-list capped at 119 for trimmed 2x150bp reads to prevent fragmentation
        # --prune-level 2 (default) balances tip removal with low-depth sensitivity
        megahit \
            -1 "$trim_r1" \
            -2 "$trim_r2" \
            -o "$out_dir" \
            -t "$THREADS" \
            -m "$MEM_BYTES" \
            --min-contig-len 1000 \
            --k-list 21,29,39,59,79,99,119
    else
        echo "Found completed MEGAHIT assembly for $sample, skipping."
    fi

    CONTIG_FILES+=("$contig_fasta")
done

if [ ${#CONTIG_FILES[@]} -eq 0 ]; then
    echo "Error: No contig files were produced or found."
    exit 1
fi

# === Step 2: Concatenation and Filtering for VAMB ===
echo "--- Step 2: Concatenating assemblies (>= ${MIN_CONTIG_SIZE} bp) ---"
if [ -f "$VAMB_IN" ]; then
    rm -f "$VAMB_IN"
fi
concatenate.py "$VAMB_IN" "${CONTIG_FILES[@]}" -m "$MIN_CONTIG_SIZE"

# === Step 3: Mapping Trimmed Reads with strobealign ===
echo "--- Step 3: Mapping trimmed reads with strobealign ---"
for trim_r1 in "$TRIM_DIR"/*_R1.trimmed.fastq.gz; do
    [ -e "$trim_r1" ] || continue
    seqfile=$(basename "$trim_r1")
    sample=${seqfile%%_R1.trimmed.fastq.gz}
    trim_r2="${TRIM_DIR}/${sample}_R2.trimmed.fastq.gz"

    if [ -f "$trim_r2" ]; then
        out_tsv="${AEMB_DIR}/${sample}.tsv"

        if [ ! -f "$out_tsv" ]; then
            echo "Mapping trimmed reads for $sample with strobealign..."
            strobealign -t "$THREADS" --aemb "$VAMB_IN" "$trim_r1" "$trim_r2" > "$out_tsv"
        else
            echo "Found aemb file for $sample, skipping mapping."
        fi
    else
        echo "Warning: Missing trimmed R2 file for $sample in mapping step, skipping..."
    fi
done

# === Step 4: Merge Abundance Files ===
echo "--- Step 4: Merging aemb files into $ABUNDANCE_TSV ---"
if [ -f "$ABUNDANCE_TSV" ]; then
    rm -f "$ABUNDANCE_TSV"
fi
merge_aemb.py "$AEMB_DIR" "$ABUNDANCE_TSV"

# === Step 5: Taxonomic Classification with DIAMOND & TaxonKit ===
echo "--- Step 5: Running taxonomy classification with DIAMOND ---"
if [ ! -f "$DIAMOND_OUT" ]; then
    echo "Mapping contigs against protein reference..."
    diamond blastx \
        -d "$DIAMOND_DB" \
        -q "$VAMB_IN" \
        -o "$DIAMOND_OUT" \
        -f 6 qseqid staxids bitscore \
        --max-target-seqs 1 \
        --max-hsps 1 \
        --evalue 1e-5 \
        --threads "$MAX_THREADS" \
        -b 12.0 \
        -c 1
else
    echo "DIAMOND output already exists, skipping search."
fi

echo "--- Formatting taxonomy for TaxVAMB ---"
if [ ! -f "$VAMB_TAX" ]; then
    TMP_TAX="${TAX_DIR}/tax_tmp.tsv"

    # Extract first taxid if multiple are semicolon-separated, then resolve lineage
    awk -F'\t' 'BEGIN{OFS="\t"} {split($2, a, ";"); $2 = a[1]; if ($2 != "") print $0}' "$DIAMOND_OUT" | \
        taxonkit lineage -i 2 | \
        taxonkit reformat -i 4 -f "{k};{p};{c};{o};{f};{g};{s}" -F | \
        awk -F'\t' '{
            tax=$5;
            gsub(/;+$/, "", tax);
            if(tax != "") {
                print $1 "\t" tax
            }
        }' > "$TMP_TAX"

    echo -e "contigs\tpredictions" > "$VAMB_TAX"
    awk -F'\t' '
        NR==FNR { tax[$1] = $2; next }
        FNR>1 {
            if ($1 in tax) {
                print $1 "\t" tax[$1]
            } else {
                print $1 "\tUnclassified"
            }
        }
    ' "$TMP_TAX" "$ABUNDANCE_TSV" >> "$VAMB_TAX"

    rm -f "$TMP_TAX"
else
    echo "Taxonomy file for VAMB already exists."
fi

# === Step 6: Run TaxVAMB ===
echo "--- Step 6: Starting TaxVAMB Binning ---"
if [ -d "$TAXVAMB_OUT" ]; then
    echo "Removing previous $TAXVAMB_OUT..."
    rm -rf "$TAXVAMB_OUT"
fi

vamb bin taxvamb \
    -m "$MIN_CONTIG_SIZE" \
    -p "$THREADS" \
    --cuda \
    --outdir "$TAXVAMB_OUT" \
    --fasta "$VAMB_IN" \
    --abundance_tsv "$ABUNDANCE_TSV" \
    --taxonomy "$VAMB_TAX"

echo "TaxVAMB binning completed."

# === Step 7: Extract FASTA Bins ===
echo "--- Step 7: Extracting FASTA bins from TaxVAMB clusters ---"
if [ -d "$BINS_DIR" ]; then
    rm -rf "$BINS_DIR"
fi
mkdir -p "$BINS_DIR"

python3 - << 'EOF'
import gzip
import os
import sys

outdir = os.environ["TAXVAMB_OUT"]
vamb_in = os.environ["VAMB_IN"]
bins_dir = os.environ["BINS_DIR"]

cluster_file = os.path.join(outdir, 'vaevae_clusters_split.tsv')
if not os.path.exists(cluster_file):
    cluster_file = os.path.join(outdir, 'clusters.tsv')

if not os.path.exists(cluster_file):
    sys.exit(f"Error: Cluster file not found in {outdir}")

print(f"Reading clusters from {cluster_file}...")
contig_to_bin = {}
bin_names = set()

with open(cluster_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('clustername'):
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            b_id, c_id = parts[0], parts[1]
            contig_to_bin[c_id] = b_id
            bin_names.add(b_id)

print(f"Allocating contigs for {len(bin_names)} bins...")
bin_records = {b_id: [] for b_id in bin_names}

print(f"Streaming sequences from {vamb_in}...")
with gzip.open(vamb_in, 'rt') as f:
    cur_header = None
    cur_seq = []
    
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if cur_header and cur_header in contig_to_bin:
                target_bin = contig_to_bin[cur_header]
                bin_records[target_bin].append(f">{cur_header}\n{''.join(cur_seq)}")
            cur_header = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line)
            
    if cur_header and cur_header in contig_to_bin:
        target_bin = contig_to_bin[cur_header]
        bin_records[target_bin].append(f">{cur_header}\n{''.join(cur_seq)}")

print(f"Writing .fna bin files into {bins_dir}...")
written = 0
for b_id, records in bin_records.items():
    if records:
        out_path = os.path.join(bins_dir, f"{b_id}.fna")
        with open(out_path, 'w') as out:
            out.write('\n'.join(records) + '\n')
        written += 1

print(f"Finished extracting {written} FASTA bins.")
EOF

# === Step 8: Quality Control with CheckM2 ===
echo "--- Step 8: Starting Quality Assessment with CheckM2 ---"

mamba deactivate
mamba activate checkm2

if [ -d "$CHECKM2_OUT" ]; then
    echo "Clearing previous CheckM2 directory..."
    rm -rf "$CHECKM2_OUT"
fi

checkm2 predict \
    --threads "$MAX_THREADS" \
    --input "$BINS_DIR" \
    -x fna \
    --output-directory "$CHECKM2_OUT" \
    --database_path "$CHECKM2_DB"

echo "Pipeline finished successfully!"
echo "CheckM2 results: ${CHECKM2_OUT}/quality_report.tsv"
