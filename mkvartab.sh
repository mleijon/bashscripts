#!/bin/bash

# Detta skript skapar en tab-separerad matris över alla varianter i varje segment
# Kolumner: Sample name, 1, 2, 3 ... referensens längd
# Data: '-', 'S:REF>ALT' eller 'aaREF>aaALT (REF>ALT)'
# INFO: Ignorerar frameshifts och inframe-indels men loggar dem till terminalen.

echo "=== Genererar Variant-Matriser för alla segment ==="

for seg_dir in results_*; do
    # Hoppa över om det inte är en mapp
    if [ ! -d "$seg_dir" ]; then continue; fi
    
    echo "Bearbetar $seg_dir..."
    
    # Hitta alla SnpEff-annoterade VCF-filer i detta segment
    vcfs=$(find "$seg_dir" -name "snvs.annotated.vcf" 2>/dev/null)
    if [ -z "$vcfs" ]; then
        echo "  -> Inga annoterade VCF-filer hittades. Hoppar över."
        continue
    fi
    
    # 1. Hämta referensens längd dynamiskt
    # Vi plockar den första VCF-filen vi hittar och letar efter length= i headern
    first_vcf=$(echo "$vcfs" | head -n 1)
    ref_len=$(grep -m 1 "^##contig=<" "$first_vcf" | sed -n 's/.*length=\([0-9]*\).*/\1/p')
    
    # Fallback: Om VCF-headern saknar längd, kolla i chrom.size
    if [ -z "$ref_len" ] && [ -f "$seg_dir/chrom.size" ]; then
        ref_len=$(awk '{print $2}' "$seg_dir/chrom.size" | head -n 1)
    fi
    
    if [ -z "$ref_len" ]; then
        echo "  -> FEL: Kunde inte fastställa referensens längd. Hoppar över."
        continue
    fi
    
    echo "  -> Referenslängd hittad: $ref_len baser."
    
    out_file="$seg_dir/variant_matrix.tsv"
    
    # 2. Använd AWK för att bygga och skriva ut matrisen
    # Vi skickar in alla VCF-filer som argument till AWK samtidigt
    awk -v ref_len="$ref_len" '
    BEGIN { 
        FS="\t"; OFS="\t" 
    }
    
    # Körs på allra första raden i varje ny fil vi läser
    FNR==1 {
        # Sökvägen ser ut såhär: results_1_PB2/100-J3PC-GD-OCV/S97/...
        # Vi klyver den vid "/" för att få ut provnamnet i position 2
        split(FILENAME, path_parts, "/")
        sample = path_parts[2]
        samples[sample] = 1
    }
    
    # Hoppa över alla header-rader inuti filerna
    /^#/ { next }
    
    {
        pos = $2
        ref = $4
        alt = $5
        info = $8
        
        # Extrahera Allele Frequency (AF)
        af_val = "N/A"
        if (match(info, /AF=[0-9\.eE\-]+/)) {
            af_val = substr(info, RSTART+3, RLENGTH-3)
        }
        
        # Hitta ANN-fältet från SnpEff
        ann_idx = index(info, "ANN=")
        if (ann_idx > 0) {
            # Klipp ut strängen som börjar efter "ANN=" fram till nästa semikolon
            ann_sub = substr(info, ann_idx + 4)
            semi_idx = index(ann_sub, ";")
            if (semi_idx > 0) ann_sub = substr(ann_sub, 1, semi_idx - 1)
            
            # SnpEff fält är separerade med pipe "|"
            # Fält 2: Typ av mutation (t.ex. synonymous_variant)
            # Fält 11: HGVS.p (t.ex. p.Pro43Ser)
            split(ann_sub, ann_arr, "|")
            annotation = ann_arr[2]
            hgvs_p = ann_arr[11]
            
            val = "" # Återställ värdet för denna rad
            
            # Filtrera och formatera baserat på SnpEff-annotering
            if (annotation ~ /frameshift_/) {
                print "  -> [OBS] Frameshift hittad och ignorerad i matrisen: Prov " sample ", Pos " pos " (" ref ">" alt ")" > "/dev/stderr"
            } else if (annotation ~ /inframe_/) {
                print "  -> [OBS] Inframe-indel hittad och ignorerad i matrisen: Prov " sample ", Pos " pos " (" ref ">" alt ")" > "/dev/stderr"
            } else if (annotation ~ /synonymous_variant/) {
                val = "S:" ref ">" alt
            } else if (annotation ~ /missense_variant/) {
                # Trimma bort "p." från början av p.Pro43Ser
                sub(/^p\./, "", hgvs_p)
                
                # POSIX-säkert sätt att plocka ut aminosyrorna (runt siffrorna)
                match(hgvs_p, /[0-9]+/)
                if (RSTART > 0) {
                    aa_ref = substr(hgvs_p, 1, RSTART-1)
                    aa_alt = substr(hgvs_p, RSTART+RLENGTH)
                    val = aa_ref ">" aa_alt " (" ref ">" alt ")"
                } else {
                    # Om formatet är oväntat, skriv ut hela
                    val = hgvs_p " (" ref ">" alt ")"
                }
            } else {
                # Fallback för andra varianter (t.ex. stop_gained)
                val = annotation " (" ref ">" alt ")"
            }
        } else {
            val = "Unknown (" ref ">" alt ")"
        }
        
        # Om val inte är tom, lägg till den i cellen
        if (val != "") {
            # Lägg till Allele Frequency (AF) i slutet av strängen
            val = val "; " af_val
            
            # Om flera varianter hamnar på exakt samma position, separera med kommatecken
            if (cell[sample, pos] != "") {
                cell[sample, pos] = cell[sample, pos] "," val
            } else {
                cell[sample, pos] = val
            }
        }
    }
    
    END {
        # 1. Bygg och skriv ut rubrikraden (Header)
        printf "Sample name"
        for (p = 1; p <= ref_len; p++) {
            printf "\t%d", p
        }
        printf "\n"
        
        # 2. Sortera proverna alfabetiskt för snyggare utskrift
        n_samples = 0
        for (s in samples) {
            sorted_samples[++n_samples] = s
        }
        
        for (i = 1; i <= n_samples; i++) {
            for (j = i + 1; j <= n_samples; j++) {
                if (sorted_samples[i] > sorted_samples[j]) {
                    tmp = sorted_samples[i]
                    sorted_samples[i] = sorted_samples[j]
                    sorted_samples[j] = tmp
                }
            }
        }
        
        # 3. Skriv ut varje provrad
        for (i = 1; i <= n_samples; i++) {
            s = sorted_samples[i]
            printf "%s", s
            for (p = 1; p <= ref_len; p++) {
                val = cell[s, p]
                if (val == "") val = "-"
                printf "\t%s", val
            }
            printf "\n"
        }
    }
    ' $vcfs > "$out_file"
    
    echo "  -> Sparade tabell: $out_file"
done

echo "=== Alla matriser färdiga! ==="
