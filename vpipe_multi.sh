#!/bin/bash
set -ueo pipefail

# ==============================================================================
# V-pipe Multi-Segment Controller (Dynamic ID & Portable Assets)
# Usage: vpipe_multi.sh [AF_THRESHOLD] (Default: 0.01)
# ==============================================================================

MIN_AF=${1:-0.01}
MASTER_DASHBOARD="dashboard.html"
ASSETS_DIR="dashboard_pages"

# Dynamically determine REFERENCE_ID from the first available segment fasta
FIRST_REF=$(find segments -name "reference.fasta" | head -n 1)
if [ -n "$FIRST_REF" ]; then
    REFERENCE_ID=$(grep ">" "$FIRST_REF" | head -n 1 | cut -d'_' -f1,2 | sed 's/^>//')
else
    REFERENCE_ID="Multi-Segment_Virus"
fi

echo "=== Starting Multi-Segment Run for $REFERENCE_ID ==="

# Prepare portable assets directory
mkdir -p "$ASSETS_DIR"
rm -f "$ASSETS_DIR"/*.html

# 1. Loop through segment folders
for seg_path in segments/*; do
    [ -d "$seg_path" ] || continue
    seg_name=$(basename "$seg_path")

    echo "------------------------------------------------------------------"
    echo " PROCESSING SEGMENT: $seg_name"
    echo "------------------------------------------------------------------"

    mkdir -p references
    rm -f references/*
    cp "$seg_path/reference.fasta" references/reference.fasta
    cp "$seg_path/reference.gff3" references/reference.gff3
    cp "$seg_path/primers.bed" references/primers.bed

    rm -rf results
    vpipe_auto.sh "$MIN_AF"

    # Move logs and reports into results before renaming
    mv snpEff_summary.html snpEff_genes.txt logs results/ 2>/dev/null || true

    echo "Saving results for $seg_name..."
    rm -rf "results_$seg_name"
    mv results "results_$seg_name"
done

# 2. Generate Master Dashboard (Portable & English)
echo "=== Generating Portable Master Dashboard ==="

cat <<EOF > "$MASTER_DASHBOARD"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>$REFERENCE_ID - Genome Dashboard</title>
    <style>
        body { font-family: sans-serif; margin: 0; display: flex; height: 100vh; background: #f8f9fa; }
        #sidebar { width: 360px; background: #2c3e50; color: white; overflow-y: auto; border-right: 1px solid #ddd; }
        #sidebar h2 { padding: 20px; margin: 0; background: #1a252f; font-size: 1.1em; text-align: center; border-bottom: 2px solid #3498db; }
        .sample-group { border-bottom: 1px solid #34495e; }
        .sample-header { 
            padding: 12px 20px; background: #34495e; cursor: pointer; font-weight: bold;
            display: flex; justify-content: space-between; align-items: center;
        }
        .sample-header:hover { background: #4e6a85; }
        .segment-list { display: none; background: #1a252f; list-style: none; margin: 0; padding: 0; }
        .segment-list.open { display: block; }
        .segment-link { 
            padding: 10px 30px; display: flex; justify-content: space-between; 
            text-decoration: none; color: #bdc3c7; font-size: 0.85em; cursor: pointer;
        }
        .segment-link:hover { color: white; background: #2980b9; }
        .segment-link.active { color: white; background: #3498db; font-weight: bold; }
        .badge { background: #e67e22; padding: 2px 7px; border-radius: 8px; font-size: 0.8em; }
        .badge-total { background: #95a5a6; }
        #viewer { flex-grow: 1; display: flex; flex-direction: column; }
        iframe { width: 100%; height: 100%; border: none; background: white; }
        .msg { padding: 50px; text-align: center; color: #7f8c8d; }
    </style>
</head>
<body>
    <div id="sidebar">
        <h2>🔬 $REFERENCE_ID Dashboard</h2>
EOF

# Identify unique sample IDs
sample_ids=$(find results_* -maxdepth 1 -type d -name "23-*" 2>/dev/null | sed 's|.*/||' | sort -u)

for s_id in $sample_ids; do
    s_tag=$(find results_* -maxdepth 2 -type d -path "*/$s_id/S*" -print -quit 2>/dev/null | sed 's|.*/||')
    
    total_var=0
    for vcf in results_*/$s_id/$s_tag/variants/SNVs/snvs.final.vcf; do
        [ -f "$vcf" ] && total_var=$((total_var + $(awk '!/^#/{c++} END{print c+0}' "$vcf")))
    done

    echo "        <div class=\"sample-group\">
            <div class=\"sample-header\" onclick=\"toggle('$s_id')\">
                <span>$s_id ($s_tag)</span>
                <span class=\"badge badge-total\">$total_var var</span>
            </div>
            <div id=\"list-$s_id\" class=\"segment-list\">" >> "$MASTER_DASHBOARD"

    for seg_dir in results_*; do
        [ -d "$seg_dir" ] || continue
        seg_name=${seg_dir#results_}
        
        vcf_file="$seg_dir/$s_id/$s_tag/variants/SNVs/snvs.final.vcf"
        seg_var=$( [ -f "$vcf_file" ] && awk '!/^#/{c++} END{print c+0}' "$vcf_file" || echo 0 )
        
        orig_html="$seg_dir/$s_id/$s_tag/visualization/snv_calling.html"

        if [ -f "$orig_html" ]; then
            # Copy file to assets folder for portability
            portable_file="${s_id}_${seg_name}.html"
            cp "$orig_html" "$ASSETS_DIR/$portable_file"
            
            echo "                <a class=\"segment-link\" onclick=\"load(this, '$ASSETS_DIR/$portable_file')\">
                    <span>Segment: $seg_name</span>
                    <span class=\"badge\">$seg_var var</span>
                </a>" >> "$MASTER_DASHBOARD"
        fi
    done

    echo "            </div>
        </div>" >> "$MASTER_DASHBOARD"
done

cat <<EOF >> "$MASTER_DASHBOARD"
    </div>
    <div id="viewer">
        <div id="welcome" class="msg">
            <h2>Analysis Complete</h2>
            <p>Select a sample from the sidebar and click a segment to view its variants.</p>
        </div>
        <iframe id="frame" src="about:blank"></iframe>
    </div>
    <script>
        function toggle(id) {
            const list = document.getElementById('list-' + id);
            list.classList.toggle('open');
        }
        function load(el, path) {
            document.getElementById('welcome').style.display = 'none';
            document.getElementById('frame').src = path;
            document.querySelectorAll('.segment-link').forEach(l => l.classList.remove('active'));
            el.classList.add('active');
        }
    </script>
</body>
</html>
EOF

echo "=== Process Complete! Open ./$MASTER_DASHBOARD to view results ==="
