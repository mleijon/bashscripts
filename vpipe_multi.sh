#!/bin/bash
set -ueo pipefail

# ==============================================================================
# V-pipe Multi-Segment Controller
# Usage: ./vpipe_multi.sh [AF_THRESHOLD] (Default: 0.01)
# ==============================================================================

MIN_AF=${1:-0.01}
MASTER_DASHBOARD="dashboard.html"

echo "=== Starting Multi-Segment V-pipe Run ==="

# 1. Loopa igenom alla segmentmappar i mappen 'segments'
for seg_path in segments/*; do
    [ -d "$seg_path" ] || continue
    seg_name=$(basename "$seg_path")

    echo ""
    echo "=================================================================="
    echo " PROCESSING SEGMENT: $seg_name"
    echo "=================================================================="
    echo ""

    # 2. Rensa gamla centrala referensfiler och förbered mappen
    mkdir -p references
    rm -f references/*

    # 3. Kopiera in detta segments specifika filer till references/
    # Notera: Förväntar sig namnen reference.fasta, reference.gff3, primers.bed
    cp "$seg_path/reference.fasta" references/reference.fasta
    cp "$seg_path/reference.gff3" references/reference.gff3
    cp "$seg_path/primers.bed" references/primers.bed

    # 4. Se till att V-pipe börjar med en tom/ren results-mapp för detta segment
    rm -rf results

    # 5. Kör det vanliga automatiserade skriptet för detta segment
    ./vpipe_auto.sh "$MIN_AF"

    # 6. Flytta undan resultatet till en segmentspecifik mapp
    echo "Saving results for $seg_name..."
    rm -rf "results_$seg_name"
    mv results "results_$seg_name"
done

# 7. Generera den slutgiltiga två-nivåers master-dashboarden
echo ""
echo "=== Generating Master Multi-Segment Dashboard ==="

cat <<EOF > "$MASTER_DASHBOARD"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>H5N1 Genome Dashboard (All Segments)</title>
    <style>
        body { font-family: sans-serif; margin: 0; display: flex; height: 100vh; background: #f8f9fa; }
        #sidebar { width: 350px; background: #2c3e50; color: white; overflow-y: auto; }
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
        <h2>🔬 H5N1 Multi-Segment Analysis</h2>
EOF

# Hitta alla unika prov-IDn som har körts genom att titta i en av resultatmapparna
# Det här gör att skriptet dynamiskt hittar dina "23-*" prover.
sample_ids=$(find results_* -maxdepth 1 -type d -name "23-*" 2>/dev/null | sed 's|.*/||' | sort -u)

for s_id in $sample_ids; do
    # Hitta S-taggen (t.ex. S1, S2) för detta prov genom att titta i första bästa resultatmapp
    s_tag=$(find results_* -maxdepth 2 -type d -path "*/$s_id/S*" -print -quit 2>/dev/null | sed 's|.*/||')
    
    # Räkna totalt antal varianter för detta prov över ALLA segmentmappar
    total_var=0
    for vcf in results_*/$s_id/$s_tag/variants/SNVs/snvs.final.vcf; do
        [ -f "$vcf" ] && total_var=$((total_var + $(grep -v "^#" "$vcf" | wc -l)))
    done

    echo "        <div class=\"sample-group\">
            <div class=\"sample-header\" onclick=\"toggle('$s_id')\">
                <span>$s_id ($s_tag)</span>
                <span class=\"badge badge-total\">$total_var var</span>
            </div>
            <div id=\"list-$s_id\" class=\"segment-list\">" >> "$MASTER_DASHBOARD"

    # Loopa igenom resultatmapparna för att bygga länkarna till varje enskilt segment
    for seg_dir in results_*; do
        [ -d "$seg_dir" ] || continue
        seg_name=${seg_dir#results_}
        
        vcf_file="$seg_dir/$s_id/$s_tag/variants/SNVs/snvs.final.vcf"
        seg_var=0
        [ -f "$vcf_file" ] && seg_var=$(grep -v "^#" "$vcf_file" | wc -l)
        
        # Sökväg till den specifika HTML-vyn för detta segment
        rel_html="$seg_dir/$s_id/$s_tag/visualization/snv_calling.html"

        echo "                <a class=\"segment-link\" onclick=\"load(this, '$rel_html')\">
                    <span>Segment: $seg_name</span>
                    <span class=\"badge\">$seg_var var</span>
                </a>" >> "$MASTER_DASHBOARD"
    done

    echo "            </div>
        </div>" >> "$MASTER_DASHBOARD"
done

cat <<EOF >> "$MASTER_DASHBOARD"
    </div>
    <div id="viewer">
        <div id="welcome" class="msg">
            <h2>All 8 Segments Processed</h2>
            <p>Select a sample and click on a specific segment to load its custom alignment track.</p>
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

echo "=== Done! Open ./dashboard.html to view your results ==="
