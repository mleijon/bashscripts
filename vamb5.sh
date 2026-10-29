#!/bin/bash
set -euo pipefail

# === Configuration ===
THREADS=48
MEM=500
MIN_CONTIG_SIZE=2000

# Directories and files
FA_DIR="./fa"
ASSEMBLY_BASE="./assemblies"
VAMB_IN="./contigs.fna.gz"
AEMB_DIR="./aemb"
ABUNDANCE_TSV="./abundance.tsv"
TAX_DIR="./taxonomy"
DIAMOND_OUT="${TAX_DIR}/diamond_out.tsv"
VAMB_TAX="${TAX_DIR}/taxonomy.tsv"
TAXVAMB_OUT="./out_taxvamb"
BINS_DIR="./bins_fasta"
CHECKM2_OUT="./checkm2_out"

# Databases
DIAMOND_DB="/mnt/micke_ssd/resources/nr_cluster_seq_2026_274.0.dmnd"
# Replace with your actual CheckM2 database path on the SSD
CHECKM2_DB="/mnt/micke_ssd/resources/CheckM2_database/uniref100.KO.1.dmnd" 

echo "Activating Conda environment vamb5..."
source ~/miniforge3/etc/profile.d/mamba.sh 
mamba activate vamb5

mkdir -p "$ASSEMBLY_BASE"
mkdir -p "$AEMB_DIR"
mkdir -p "$TAX_DIR"

# === Step 1: Assembly with MetaSPAdes ===
echo "--- Starting MetaSPAdes for samples in $FA_DIR ---"
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
    
    out_dir="${ASSEMBLY_BASE}/${sample}"
    contig_fasta="${out_dir}/contigs.fasta"
    
    if [ ! -f "$contig_fasta" ]; then
        echo "Running assembly for $sample..."
        spades.py --meta -1 "$r1" -2 "$r2" -o "$out_dir" -t "$THREADS" -m "$MEM"
    else
        echo "Found completed assembly for $sample, skipping SPAdes."
    fi
    
    CONTIG_FILES+=("$contig_fasta")
done

# === Step 2: Concatenation for VAMB ===
echo "--- Concatenating assemblies ---"
if [ -f "$VAMB_IN" ]; then
    rm "$VAMB_IN"
fi
concatenate.py "$VAMB_IN" "${CONTIG_FILES[@]}" -m "$MIN_CONTIG_SIZE"

# === Step 3: Mapping and abundance calculation with strobealign ===
echo "--- Starting mapping with strobealign ---"
for r1 in "$FA_DIR"/*_R1*.fastq.gz; do
    [ -e "$r1" ] || continue
    
    seqfile=$(basename "$r1")
    sample=${seqfile%%_*}
    ext=${seqfile#*"$sample"}
    r2="${FA_DIR}/${sample}${ext/_R1/_R2}"
    
    if [ -f "$r2" ]; then
        out_tsv="${AEMB_DIR}/${sample}.tsv"
        
        if [ ! -f "$out_tsv" ]; then
            echo "Mapping $sample with strobealign..."
            strobealign -t "$THREADS" --aemb "$VAMB_IN" "$r1" "$r2" > "$out_tsv"
        else
            echo "Found aemb file for $sample, skipping mapping."
        fi
    else
        echo "Warning: Missing R2 file for $sample in mapping step, skipping..."
    fi
done

# === Step 4: Merge abundance files ===
echo "--- Merging aemb files into $ABUNDANCE_TSV ---"
if [ -f "$ABUNDANCE_TSV" ]; then
    rm "$ABUNDANCE_TSV"
fi
merge_aemb.py "$AEMB_DIR" "$ABUNDANCE_TSV"

# === Step 5: Taxonomic classification with DIAMOND ===
echo "--- Running taxonomy with DIAMOND ---"
if [ ! -f "$DIAMOND_OUT" ]; then
    echo "Mapping contigs against protein database with DIAMOND..."
    diamond blastx \
        -d "$DIAMOND_DB" \
        -q "$VAMB_IN" \
        -o "$DIAMOND_OUT" \
        -f 6 qseqid staxids bitscore \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --threads "$THREADS"
else
    echo "DIAMOND output already exists, skipping search."
fi

echo "--- Formatting taxonomy for TaxVAMB ---"
if [ ! -f "$VAMB_TAX" ]; then
    TMP_TAX="${TAX_DIR}/tax_tmp.tsv"
    
    cat "$DIAMOND_OUT" | \
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
    
    rm "$TMP_TAX"
else
    echo "Taxonomy file for VAMB already exists."
fi

# === Step 6: Run TaxVAMB ===
echo "--- Starting TaxVAMB ---"

if [ -d "$TAXVAMB_OUT" ]; then
    echo "Removing existing $TAXVAMB_OUT so VAMB can create a new one..."
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

echo "TaxVAMB binning completed!"

# === Step 7: Extract Bins to FASTA files ===
echo "--- Extracting FASTA bins from TaxVAMB output ---"
if [ -d "$BINS_DIR" ]; then
    rm -rf "$BINS_DIR"
fi
mkdir -p "$BINS_DIR"

# Using a pure python snippet to read vaevae_clusters_split.tsv and split contigs.fna.gz into individual bins
python -c "
import gzip
import os

print('Parsing vaevae_clusters_split.tsv...')
cluster_map = {}
cluster_file = os.path.join('$TAXVAMB_OUT', 'vaevae_clusters_split.tsv')
if not os.path.exists(cluster_file):
    print(f'Error: Could not find {cluster_file}')
    exit(1)

with open(cluster_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('clustername'):
            continue # Skip header and empty lines

        c_id, contig_id = line.split('\t')
        if c_id not in cluster_map:
            cluster_map[c_id] = []
        cluster_map[c_id].append(contig_id)

contig_to_c = {contig: c_id for c_id, contigs in cluster_map.items() for contig in contigs}
bins = {c_id: [] for c_id in cluster_map}

print('Extracting sequences from $VAMB_IN...')
with gzip.open('$VAMB_IN', 'rt') as f:
    header = ''
    seq = []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if header and header in contig_to_c:
                bins[contig_to_c[header]].append('>' + header + '\n' + ''.join(seq))
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if header and header in contig_to_c:
        bins[contig_to_c[header]].append('>' + header + '\n' + ''.join(seq))

print('Writing .fna files...')
extracted_count = 0
for c_id, contig_seqs in bins.items():
    if contig_seqs: # Only write non-empty bins
        with open(os.path.join('$BINS_DIR', c_id + '.fna'), 'w') as out:
            out.write('\n'.join(contig_seqs) + '\n')
        extracted_count += 1
print(f'Successfully extracted {extracted_count} valid bins.')
"

# === Step 8: Quality Control with CheckM2 ===
echo "--- Starting Quality Control with CheckM2 ---"

# Deactivate vamb5 and activate the dedicated checkm2 environment
mamba deactivate
mamba activate checkm2

if [ -d "$CHECKM2_OUT" ]; then
    echo "Removing existing CheckM2 output directory..."
    rm -rf "$CHECKM2_OUT"
fi

# Run CheckM2 prediction
checkm2 predict \
    --threads "$THREADS" \
    --input "$BINS_DIR" \
    -x fna \
    --output-directory "$CHECKM2_OUT" \
    --database_path "$CHECKM2_DB"

echo "Pipeline completed successfully! CheckM2 quality report is located in $CHECKM2_OUT/quality_report.tsv"
