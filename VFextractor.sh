#!/bin/bash
# Safety first: exit on error or undefined variables
set -ueo pipefail

# ==============================================================================
# MASTER VIRULENCE FACTOR EXTRACTOR (Modular Folder Edition)
# ==============================================================================

# --- PARAMETERS ---
MIN_COVERAGE=70             # Min % coverage for a gene to be considered "Present"
TANDEM_DIST=2500            # Max distance (bp) for fnbA/B tandem check
REF_DIR="vfextractor_refs"   # Folder containing your reference .fasta files
# ------------------

SUMMARY_FILE="run_summary_$(date +%Y%m%d_%H%M).txt"

# Ensure the reference directory exists
if [ ! -d "$REF_DIR" ]; then
    echo "❌ Error: Directory '$REF_DIR' not found."
    exit 1
fi

# Get list of references from the folder
mapfile -t REFS < <(ls "$REF_DIR"/*.fasta 2>/dev/null)

if [ ${#REFS[@]} -eq 0 ]; then
    echo "❌ Error: No .fasta files found in '$REF_DIR'."
    exit 1
fi

echo "🚀 Starting Modular Extraction (Folder: $REF_DIR)"
echo "🧬 Found ${#REFS[@]} reference genes."

# 1. INITIALIZE SUMMARY
{
    echo "GENOMIC EXTRACTION SUMMARY"
    echo "Date: $(date)"
    echo "Threshold: ${MIN_COVERAGE}% coverage required"
    echo "Reference Folder: $REF_DIR"
    echo "------------------------------------------------"
} > "$SUMMARY_FILE"

# Clean up any existing result files from previous runs
for ref_path in "${REFS[@]}"; do
    gene_name=$(basename "$ref_path" .fasta)
    rm -f "all_${gene_name}_DNA.fasta" "all_${gene_name}_PROT.fasta"
done

# 2. LOOP THROUGH ISOLATE FILES (All .fasta files in current dir)
for isolate in *.fasta; do
    # Skip files that are actually references or result files
    [[ "$isolate" == "$REF_DIR"* ]] && continue
    [[ "$isolate" == "all_"* ]] && continue
    
    sample_id=$(basename "$isolate" .fasta)
    echo "🧪 Processing Sample: $sample_id"
    echo -e "\nSample: $sample_id" >> "$SUMMARY_FILE"

    # Create a working copy for masking
    working_genome="${sample_id}_temp_masked.fasta"
    cp "$isolate" "$working_genome"

    # Trackers for Tandem Check (e.g., fnbA/B)
    fnba_contig=""; fnba_pos=0
    fnbb_contig=""; fnbb_pos=0

    for ref_path in "${REFS[@]}"; do
        gene_name=$(basename "$ref_path" .fasta)

        makeblastdb -in "$working_genome" -dbtype nucl -logfile /dev/null
        
        # BLAST search
        hit=$(blastn -query "$ref_path" -db "$working_genome" -evalue 1e-10 \
              -outfmt "6 sseqid sstart send length qlen" -max_target_seqs 1 -max_hsps 1 2>/dev/null)

        if [ -z "$hit" ]; then
            echo "    - $gene_name: ABSENT" >> "$SUMMARY_FILE"
            rm -f "${working_genome}.n"*
            continue
        fi

        contig=$(echo "$hit" | awk '{print $1}')
        s=$(echo "$hit" | awk '{print $2}')
        e=$(echo "$hit" | awk '{print $3}')
        aln=$(echo "$hit" | awk '{print $4}')
        qlen=$(echo "$hit" | awk '{print $5}')
        cov=$(( 100 * aln / qlen ))

        # Positional tracking for tandem check
        mid=$(( (s + e) / 2 ))
        [[ "$gene_name" == *"fnbA"* && "$cov" -ge "$MIN_COVERAGE" ]] && { fnba_contig=$contig; fnba_pos=$mid; }
        [[ "$gene_name" == *"fnbB"* && "$cov" -ge "$MIN_COVERAGE" ]] && { fnbb_contig=$contig; fnbb_pos=$mid; }

        # 3. THRESHOLD CHECK
        if [ "$cov" -lt "$MIN_COVERAGE" ]; then
            echo "    - $gene_name: PARTIAL ($cov% coverage) at $contig:$s-$e" >> "$SUMMARY_FILE"
            rm -f "${working_genome}.n"*
            continue
        fi

        # 4. EXTRACT & MASK
        # Handle reverse orientation automatically
        if [ "$s" -lt "$e" ]; then
            seq=$(seqkit grep -p "$contig" "$working_genome" | seqkit subseq -r "$s:$e")
            mask_s=$((s-1)); mask_e=$e
        else
            seq=$(seqkit grep -p "$contig" "$working_genome" | seqkit subseq -r "$e:$s" | seqkit seq -r -p)
            mask_s=$((e-1)); mask_e=$s
        fi
        
        # Save results (Flattened DNA & Translated Protein)
        echo "$seq" | seqkit replace -p ".*" -r "$sample_id" | seqkit seq -w 0 >> "all_${gene_name}_DNA.fasta"
        echo "$seq" | seqkit translate -f 1 | seqkit replace -p ".*" -r "$sample_id" | seqkit seq -w 0 >> "all_${gene_name}_PROT.fasta"

        # MASK the found region to avoid redundant hits in the next loop
        echo -e "$contig\t$mask_s\t$mask_e" > temp_mask.bed
        bedtools maskfasta -fi "$working_genome" -bed temp_mask.bed -fo "${working_genome}.new" 2>/dev/null
        mv "${working_genome}.new" "$working_genome"
        
        echo "    - $gene_name: PRESENT ($cov% coverage) at $contig:$s-$e" >> "$SUMMARY_FILE"
        rm -f "${working_genome}.n"* temp_mask.bed
    done

    # 5. TANDEM CHECK (Specifically for fnbA/fnbB proximity)
    if [[ -n "$fnba_contig" && -n "$fnbb_contig" ]]; then
        if [[ "$fnba_contig" == "$fnbb_contig" ]]; then
            diff=$(( fnba_pos - fnbb_pos )); abs_diff=${diff#-};
            if [ "$abs_diff" -le "$TANDEM_DIST" ]; then
                echo "    >>> NOTICE: fnbA and fnbB are TANDEM (Dist: ${abs_diff}bp)" >> "$SUMMARY_FILE"
            else
                echo "    >>> NOTICE: fnbA and fnbB are DISTANT (${abs_diff}bp)" >> "$SUMMARY_FILE"
            fi
        fi
    fi

    rm -f "$working_genome"
done

# 6. FINAL PREVALENCE SUMMARY
{
    echo -e "\n------------------------------------------------"
    echo "FINAL PREVALENCE (Passed ${MIN_COVERAGE}% threshold)"
    echo "------------------------------------------------"
    for ref_path in "${REFS[@]}"; do
        gene_name=$(basename "$ref_path" .fasta)
        count=0
        [ -f "all_${gene_name}_DNA.fasta" ] && count=$(grep -c ">" "all_${gene_name}_DNA.fasta")
        printf "%-25s : %d samples\n" "$gene_name" "$count"
    done
} >> "$SUMMARY_FILE"

echo -e "\n✅ Processing complete. Summary: $SUMMARY_FILE"
