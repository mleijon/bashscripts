#!/bin/bash
# Safety first: exit on error or undefined variables
set -ueo pipefail

# ==============================================================================
# MASTER VIRULENCE FACTOR EXTRACTOR (RBH Validation Edition)
# ==============================================================================

# --- PARAMETERS ---
MIN_COVERAGE=70             # Min % coverage for a gene to be considered "Present"
TANDEM_DIST=2500            # Max distance (bp) for fnbA/B tandem check
REF_DIR="vfextractor_refs"  # Folder containing your reference .fasta files
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

echo "🚀 Starting Reciprocal Best Hit Extraction (Folder: $REF_DIR)"
echo "🧬 Found ${#REFS[@]} reference genes."

# 1. INITIALIZE SUMMARY & COMBINED DATABASE
{
    echo "GENOMIC EXTRACTION SUMMARY (RBH Validated)"
    echo "Date: $(date)"
    echo "Threshold: ${MIN_COVERAGE}% coverage required"
    echo "Reference Folder: $REF_DIR"
    echo "------------------------------------------------"
} > "$SUMMARY_FILE"

# Clean up existing result files and build the master validation database
COMBINED_DB="combined_refs_master.fasta"
rm -f "$COMBINED_DB" "all_"*"_DNA.fasta" "all_"*"_PROT.fasta"

echo "🗄️  Building master reference database for cross-validation..."
for ref_path in "${REFS[@]}"; do
    gene_name=$(basename "$ref_path" .fasta)
    # Prefix every FASTA header with "gene_name::" so we know where it came from
    awk -v g="$gene_name" '/^>/ {print ">" g "::" substr($0,2)} !/^>/ {print}' "$ref_path" >> "$COMBINED_DB"
done
makeblastdb -in "$COMBINED_DB" -dbtype nucl -logfile /dev/null

# 2. LOOP THROUGH ISOLATE FILES
for isolate in *.fasta; do
    # Skip reference or result files
    [[ "$isolate" == "$REF_DIR"* ]] && continue
    [[ "$isolate" == "all_"* ]] && continue
    [[ "$isolate" == "$COMBINED_DB"* ]] && continue
    
    sample_id=$(basename "$isolate" .fasta)
    echo "🧪 Processing Sample: $sample_id"
    echo -e "\nSample: $sample_id" >> "$SUMMARY_FILE"

    working_genome="${sample_id}_temp_masked.fasta"
    cp "$isolate" "$working_genome"

    # Trackers for Tandem Check
    fnba_contig=""; fnba_pos=0
    fnbb_contig=""; fnbb_pos=0

    for ref_path in "${REFS[@]}"; do
        gene_name=$(basename "$ref_path" .fasta)

        makeblastdb -in "$working_genome" -dbtype nucl -logfile /dev/null
        
        # FORWARD BLAST: Find best match in the genome
        hit=$(blastn -query "$ref_path" -db "$working_genome" -evalue 1e-10 \
              -outfmt "6 sseqid sstart send length qlen bitscore qseqid" \
              -max_target_seqs 1 -max_hsps 1 2>/dev/null | sort -k6,6nr | head -n 1)

        if [ -z "$hit" ]; then
            echo "    - $gene_name: ABSENT" >> "$SUMMARY_FILE"
            rm -f "${working_genome}.n"*
            continue
        fi

        # Parse the top forward hit
        contig=$(echo "$hit" | awk '{print $1}')
        s=$(echo "$hit" | awk '{print $2}')
        e=$(echo "$hit" | awk '{print $3}')
        aln=$(echo "$hit" | awk '{print $4}')
        qlen=$(echo "$hit" | awk '{print $5}')
        variant=$(echo "$hit" | awk '{print $7}')
        
        cov=$(( 100 * aln / qlen ))

        if [ "$cov" -lt "$MIN_COVERAGE" ]; then
            echo "    - $gene_name: PARTIAL ($cov% coverage) at $contig:$s-$e [Best match: $variant]" >> "$SUMMARY_FILE"
            rm -f "${working_genome}.n"*
            continue
        fi

        # Extract sequence early for reverse-validation
        if [ "$s" -lt "$e" ]; then
            seq=$(seqkit grep -p "$contig" "$working_genome" 2>/dev/null | seqkit subseq -r "$s:$e" 2>/dev/null)
            mask_s=$((s-1)); mask_e=$e
        else
            seq=$(seqkit grep -p "$contig" "$working_genome" 2>/dev/null | seqkit subseq -r "$e:$s" 2>/dev/null | seqkit seq -r -p)
            mask_s=$((e-1)); mask_e=$s
        fi
        
        # 3. REVERSE BLAST VALIDATION (The Homology Check)
        echo "$seq" > temp_rev_query.fasta
        rev_hit=$(blastn -query temp_rev_query.fasta -db "$COMBINED_DB" -evalue 1e-10 \
                  -outfmt "6 sseqid bitscore" -max_target_seqs 1 2>/dev/null | head -n 1)
        
        if [ -n "$rev_hit" ]; then
            # Extract the true gene identity from our tagged header (e.g., extracting "sdrC" from "sdrC::var1")
            best_ref=$(echo "$rev_hit" | awk '{print $1}' | awk -F'::' '{print $1}')
            
            if [ "$best_ref" != "$gene_name" ]; then
                echo "    - $gene_name: SKIPPED (Reverse-BLAST identified this locus as $best_ref)" >> "$SUMMARY_FILE"
                rm -f "${working_genome}.n"* temp_rev_query.fasta
                continue
            fi
        fi

        # Positional tracking for tandem check
        mid=$(( (s + e) / 2 ))
        [[ "$gene_name" == *"fnbA"* ]] && { fnba_contig=$contig; fnba_pos=$mid; }
        [[ "$gene_name" == *"fnbB"* ]] && { fnbb_contig=$contig; fnbb_pos=$mid; }

        # 4. SAVE & MASK (Only executes if RBH validation passed)
        echo "$seq" | seqkit replace -p ".*" -r "$sample_id" | seqkit seq -w 0 >> "all_${gene_name}_DNA.fasta"
        echo "$seq" | seqkit translate -f 1 | seqkit replace -p ".*" -r "$sample_id" | seqkit seq -w 0 >> "all_${gene_name}_PROT.fasta"

        echo -e "$contig\t$mask_s\t$mask_e" > temp_mask.bed
        bedtools maskfasta -fi "$working_genome" -bed temp_mask.bed -fo "${working_genome}.new" 2>/dev/null
        mv "${working_genome}.new" "$working_genome"
        
        echo "    - $gene_name: PRESENT ($cov% coverage) at $contig:$s-$e [Best match: $variant]" >> "$SUMMARY_FILE"
        rm -f "${working_genome}.n"* temp_mask.bed temp_rev_query.fasta
    done

    # 5. TANDEM CHECK (fnbA/fnbB proximity)
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

# Cleanup Master DB
rm -f "$COMBINED_DB"*

# 6. FINAL PREVALENCE SUMMARY
{
    echo -e "\n------------------------------------------------"
    echo "FINAL PREVALENCE (Passed ${MIN_COVERAGE}% & RBH Check)"
    echo "------------------------------------------------"
    for ref_path in "${REFS[@]}"; do
        gene_name=$(basename "$ref_path" .fasta)
        count=0
        [ -f "all_${gene_name}_DNA.fasta" ] && count=$(grep -c ">" "all_${gene_name}_DNA.fasta")
        printf "%-25s : %d samples\n" "$gene_name" "$count"
    done
} >> "$SUMMARY_FILE"

echo -e "\n✅ Processing complete. Summary: $SUMMARY_FILE"
