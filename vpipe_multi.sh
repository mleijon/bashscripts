#!/bin/bash
#set -ueo pipefail

# ==============================================================================
# V-pipe Multi-Segment Controller (PATH-aware version)
# Usage: vpipe_multi.sh [AF_THRESHOLD] (Default: 0.01)
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
    cp "$seg_path/reference.fasta" references/reference.fasta
    cp "$seg_path/reference.gff3" references/reference.gff3
    cp "$seg_path/primers.bed" references/primers.bed

    # 4. Se till att V-pipe börjar med en tom/ren results-mapp för detta segment
    rm -rf results

    # 5. Kör det vanliga automatiserade skriptet via PATH (utan ./)
    vpipe_auto.sh "$MIN_AF"

    echo "Organizing logs and reports for $seg_name..."
    mv snpEff_summary.html snpEff_genes.txt logs results/ 2>/dev/null || true

    # 6. Flytta undan resultatet till en segmentspecifik mapp
    echo "Saving results for $seg_name..."
    rm -rf "results_$seg_name"
    mv results "results_$seg_name"
done

# 7. Generera den slutgiltiga två-nivåers master-dashboarden
echo ""
echo "=== Generating Master Multi-Segment Dashboard ==="

# Skapa sidhuvudet
cat <<EOF > "$MASTER_DASHBOARD"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>H5N1 Genome Dashboard (All Segments)</title>
    <style>
        body { font-family: sans-serif; margin: 0; display: flex; height: 100vh; background: #f8f9fa; overflow: hidden; }
        #sidebar { width: 350px; background: #2c3e50; color: white; overflow-y: auto; flex-shrink: 0; }
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
        .badge { background: #e67e22; color: white; padding: 2px 7px; border-radius: 8px; font-size: 0.8em; }
        .badge-total { background: #3498db; }
        #viewer { flex-grow: 1; display: flex; flex-direction: column; background: white; }
        iframe { width: 100%; height: 100%; border: none; }
        .msg { padding: 50px; text-align: center; color: #7f8c8d; }
    </style>
</head>
<body>
    <div id="sidebar">
        <h2>🔬 H5N1 Multi-Segment Analysis</h2>
EOF

# Hitta alla unika Sample IDs (t.ex. 23-MIK...)
sample_ids=$(find results_* -maxdepth 1 -type d -name "23-*" 2>/dev/null | sed 's|.*/||' | sort -u)

for s_id in $sample_ids; do
    # Hitta S-taggen (S1, S2 osv) för detta prov genom att titta i en av resultatmapparna
    s_tag=$(find results_* -maxdepth 2 -type d -path "*/$s_id/S*" -print -quit 2>/dev/null | sed 's|.*/||')
    
    # Beräkna totalt antal varianter över alla 8 segment för detta prov
    total_var=0
    for vcf in results_*/$s_id/$s_tag/variants/SNVs/snvs.final.vcf; do
        if [ -f "$vcf" ]; then
	    count=$(awk '!/^#/{c++} END{print c+0}' "$vcf")
            total_var=$((total_var + count))
        fi
    done

    # Skriv ut prov-headern till HTML
    cat <<EOF >> "$MASTER_DASHBOARD"
        <div class="sample-group">
            <div class="sample-header" onclick="toggle('$s_id')">
                <span>$s_id ($s_tag)</span>
                <span class="badge badge-total">$total_var var</span>
            </div>
            <div id="list-$s_id" class="segment-list">
EOF

    # Loopa igenom varje segmentsmapp och skapa länkar
    for seg_dir in results_*; do
        [ -d "$seg_dir" ] || continue
        seg_name=${seg_dir#results_}
        
        vcf_file="$seg_dir/$s_id/$s_tag/variants/SNVs/snvs.final.vcf"
        seg_var=0
	[ -f "$vcf_file" ] && seg_var=$(awk '!/^#/{c++} END{print c+0}' "$vcf_file")
        
        # Sökvägen är nu direkt relativ till working-mappen
        rel_html="$seg_dir/$s_id/$s_tag/visualization/snv_calling.html"

        if [ -f "$rel_html" ]; then
            echo "                <a class=\"segment-link\" onclick=\"load(this, '$rel_html')\">
                    <span>Segment: $seg_name</span>
                    <span class=\"badge\">$seg_var var</span>
                </a>" >> "$MASTER_DASHBOARD"
        fi
    done

    echo "            </div>
        </div>" >> "$MASTER_DASHBOARD"
done

# Avsluta HTML-filen med JavaScript
cat <<EOF >> "$MASTER_DASHBOARD"
    </div>
    <div id="viewer">
        <div id="welcome" class="msg">
            <h2>Analysen är klar</h2>
            <p>Välj ett prov i menyn till vänster och klicka på ett segment för att se dess alignment och varianter.</p>
        </div>
        <iframe id="frame" src="about:blank"></iframe>
    </div>
    <script>
        function toggle(id) {
            const list = document.getElementById('list-' + id);
            // Stäng andra listor om man vill ha det städat (valfritt)
            // document.querySelectorAll('.segment-list').forEach(el => { if(el.id !== 'list-'+id) el.classList.remove('open'); });
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

echo "=== Klart! Öppna ./$MASTER_DASHBOARD för att se resultaten ==="
