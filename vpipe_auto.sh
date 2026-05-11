#!/bin/bash

# ==============================================================================
# V-pipe Automated Pipeline: Alignment, Variant Calling, Annotation, & Dashboard
# Usage: ./vpipe_auto.sh [AF_THRESHOLD] (Default: 0.01)
# ==============================================================================

# 1. Configuration & Parameters
MIN_AF=${1:-0.01}
MIN_DP=100
REFERENCE_ID="IAV_H5N1_HA"
SNPEFF_CONFIG="snpeff.config"
WORKING_DIR=$(pwd)

echo "--- Starting Pipeline Analysis ---"
echo "Threshold: AF > $MIN_AF, DP > $MIN_DP"

# 2. Step 1: Map sequences using Symbolic Links with Clean Filenames
# Requirement: samples/ID/Tag/raw_fasta/[ID]_[Tag]_[R1/R2].fastq.gz
echo "--- Step 1: Mapping and cleaning filenames in raw_fasta folders ---"
mkdir -p samples
for f in fa/*; do
    [ -e "$f" ] || continue
    filename=$(basename "$f")
    
    # Extract components from Illumina format (e.g., 23-MIK134061_S1_L001_R1_001.fastq.gz)
    sample_id=$(echo "$filename" | cut -d'_' -f1)
    s_tag=$(echo "$filename" | cut -d'_' -f2)
    
    # Identify if it is Read 1 or Read 2
    if [[ "$filename" == *"_R1_"* ]]; then
        read_pair="R1"
    else
        read_pair="R2"
    fi
    
    # Construct the cleaned filename: 23-MIK134061_S1_R1.fastq.gz
    new_filename="${sample_id}_${s_tag}_${read_pair}.fastq.gz"
    
    # Create the required subfolder structure
    target_dir="$WORKING_DIR/samples/$sample_id/$s_tag/raw_data"
    mkdir -p "$target_dir"
    
    # Create symbolic link using the CLEANED filename
    ln -sf "$WORKING_DIR/$f" "$target_dir/$new_filename"
    echo "  Mapped: $filename -> $new_filename"
done

# 3. Step 2: Run V-pipe Core Workflow
echo "--- Step 2: Running V-pipe ---"
. ~/miniforge3/bin/activate V-pipe
./vpipe --use-conda --cores all

# 4. Step 2.5: SnpEff Database Build (with Logging)
echo "--- Step 2.5: Rebuilding SnpEff Database ---"
mkdir -p logs
rm -rf "snpeff_db/$REFERENCE_ID/snpEffectPredictor.bin"

# Redirect build output to log
snpEff build -c "$SNPEFF_CONFIG" -gff3 -noCheckCds -noCheckProtein -v "$REFERENCE_ID" > logs/snpeff_build.log 2>&1
echo "  - Build log: logs/snpeff_build.log"

# 5. Step 3: Annotation & Quality Filtering
echo "--- Step 3: Annotation & Filtering ---"
find results -name "snvs.vcf" | while read -r vcf; do
    sample_dir=$(dirname "$vcf")
    sample_name=$(echo "$vcf" | cut -d'/' -f2)
    s_folder=$(echo "$vcf" | cut -d'/' -f3)
    
    annotated_vcf="${vcf%.vcf}.annotated.vcf"
    final_vcf="${vcf%.vcf}.final.vcf"
    
    summary_html="$sample_dir/snpEff_summary.html"

    echo "Processing: $sample_name ($s_folder)"

    # Annotation: Note that -noCheck flags are excluded here
    snpEff -c "$SNPEFF_CONFIG" \
           -s "$summary_html" \
           "$REFERENCE_ID" "$vcf" > "$annotated_vcf"

    # Quality Filtering using bcftools
    bcftools filter -i "AF > $MIN_AF && DP > $MIN_DP" "$annotated_vcf" > "$final_vcf"
    
    count=$(grep -v "^#" "$final_vcf" | wc -l)
    echo "  - Result: $count high-confidence variants."
done

# 6. Step 4: Generating Results Dashboard
echo "--- Step 4: Generating Dashboard ---"
DASHBOARD_FILE="results/dashboard.html"

cat <<EOF > "$DASHBOARD_FILE"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>V-pipe Analysis Dashboard</title>
    <style>
        body { font-family: sans-serif; margin: 0; display: flex; height: 100vh; background: #f8f9fa; }
        #sidebar { width: 350px; background: #2c3e50; color: white; overflow-y: auto; display: flex; flex-direction: column; }
        #sidebar h2 { padding: 20px; margin: 0; background: #1a252f; font-size: 1.1em; text-align: center; }
        .sample-link { 
            padding: 12px 20px; border-bottom: 1px solid #34495e; cursor: pointer;
            display: flex; justify-content: space-between; align-items: center;
            text-decoration: none; color: #ecf0f1; font-size: 0.85em;
        }
        .sample-link:hover { background: #34495e; }
        .sample-link.active { background: #3498db; font-weight: bold; }
        .badge { background: #e67e22; padding: 2px 7px; border-radius: 8px; font-size: 0.8em; font-weight: bold; }
        #viewer { flex-grow: 1; display: flex; flex-direction: column; }
        iframe { width: 100%; height: 100%; border: none; background: white; }
        .msg { padding: 50px; text-align: center; color: #7f8c8d; }
    </style>
</head>
<body>
    <div id="sidebar">
        <h2>🔬 H5N1 Analysis</h2>
        <div id="menu">
EOF

find results -name "snvs.final.vcf" | sort | while read -r vcf_path; do
    s_id=$(echo "$vcf_path" | cut -d'/' -f2)
    s_t=$(echo "$vcf_path" | cut -d'/' -f3)
    rel_html="${s_id}/${s_t}/visualization/snv_calling.html"
    v_count=$(grep -v "^#" "$vcf_path" | wc -l)
    
    echo "            <a class=\"sample-link\" onclick=\"load(this, '$rel_html')\">
                <span>$s_id ($s_t)</span>
                <span class=\"badge\">$v_count var</span>
            </a>" >> "$DASHBOARD_FILE"
done

cat <<EOF >> "$DASHBOARD_FILE"
        </div>
    </div>
    <div id="viewer">
        <div id="welcome" class="msg">
            <h2>Workflow Complete</h2>
            <p>Threshold: AF > $MIN_AF | DP > $MIN_DP</p>
            <p>Select a sample to view the interactive SNV plot.</p>
        </div>
        <iframe id="frame" src="about:blank"></iframe>
    </div>
    <script>
        function load(el, path) {
            document.getElementById('welcome').style.display = 'none';
            document.getElementById('frame').src = path;
            document.querySelectorAll('.sample-link').forEach(l => l.classList.remove('active'));
            el.classList.add('active');
        }
    </script>
</body>
</html>
EOF

echo "--- Pipeline Finished: results/dashboard.html ---"
