#!/bin/bash

# Detta skript skapar en tab-separerad matris över alla varianter i varje segment
# Kolumner: Sample name, 1, 2, 3 ... referensens längd
# Data: '-', 'S:REF>ALT' eller 'aaREF>aaALT (REF>ALT)'
# INFO: Ignorerar frameshifts och inframe-indels men loggar dem till terminalen.
# SKAPAR ÄVEN: En textfil som summerar totala antalet samt de 12 specifika nukleotid-substitutionerna 
# för både synonyma och missense varianter.

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
    summary_file="$seg_dir/substitution_summary.tsv"

    # 2. Använd AWK för att bygga och skriva ut matrisen
    awk -v ref_len="$ref_len" -v summary_file="$summary_file" '
    BEGIN {
        FS="\t"; OFS="\t"
        
        # NYTT: Definiera de 12 möjliga substitutionerna i en array för att få konsekvent kolumnordning
        split("A>C A>G A>T C>A C>G C>T G>A G>C G>T T>A T>C T>G", sub_types, " ")
    }

    # Körs på allra första raden i varje ny fil vi läser
    FNR==1 {
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
        
        # NYTT: Skapa strängen för själva substitutionen (t.ex. "A>G")
        mut_type = ref ">" alt

        # Extrahera Allele Frequency (AF)
        af_val = "N/A"
        if (match(info, /AF=[0-9\.eE\-]+/)) {
            af_val = substr(info, RSTART+3, RLENGTH-3)
        }

        # Hitta ANN-fältet från SnpEff
        ann_idx = index(info, "ANN=")
        if (ann_idx > 0) {
            ann_sub = substr(info, ann_idx + 4)
            semi_idx = index(ann_sub, ";")
            if (semi_idx > 0) ann_sub = substr(ann_sub, 1, semi_idx - 1)

            split(ann_sub, ann_arr, "|")
            annotation = ann_arr[2]
            hgvs_p = ann_arr[11]

            val = "" 

            if (annotation ~ /frameshift_/) {
                print "  -> [OBS] Frameshift hittad och ignorerad i matrisen: Prov " sample ", Pos " pos " (" mut_type ")" > "/dev/stderr"
            } else if (annotation ~ /inframe_/) {
                print "  -> [OBS] Inframe-indel hittad och ignorerad i matrisen: Prov " sample ", Pos " pos " (" mut_type ")" > "/dev/stderr"
            } else if (annotation ~ /synonymous_variant/) {
                val = "S:" mut_type
                
                # NYTT: Räkna upp totalen OCH den specifika substitutionstypen
                syn_count[sample]++ 
                syn_sub_count[sample, mut_type]++ 

            } else if (annotation ~ /missense_variant/) {
                sub(/^p\./, "", hgvs_p)

                match(hgvs_p, /[0-9]+/)
                if (RSTART > 0) {
                    aa_ref = substr(hgvs_p, 1, RSTART-1)
                    aa_alt = substr(hgvs_p, RSTART+RLENGTH)
                    val = aa_ref ">" aa_alt " (" mut_type ")"
                } else {
                    val = hgvs_p " (" mut_type ")"
                }
                
                # NYTT: Räkna upp totalen OCH den specifika substitutionstypen
                mis_count[sample]++ 
                mis_sub_count[sample, mut_type]++ 

            } else {
                val = annotation " (" mut_type ")"
            }
        } else {
            val = "Unknown (" mut_type ")"
        }

        if (val != "") {
            val = val "; " af_val
            if (cell[sample, pos] != "") {
                cell[sample, pos] = cell[sample, pos] "," val
            } else {
                cell[sample, pos] = val
            }
        }
    }

    END {
        # 1. Bygg och skriv ut rubrikraden (Header) för matrisen
        printf "Sample name"
        for (p = 1; p <= ref_len; p++) {
            printf "\t%d", p
        }
        printf "\n"

        # 2. Sortera proverna alfabetiskt
        n_samples = 0
        for (s in samples) {
            sorted_samples[++n_samples] = s
        }
        for (i = 1; i <= n_samples; i++) {
            for (j = i + 1; j <= n_samples; j++) {
                if (sorted_samples[i] > sorted_samples[j]) {
                    tmp = sorted_samples[i]; sorted_samples[i] = sorted_samples[j]; sorted_samples[j] = tmp
                }
            }
        }

        # 3. Skriv ut varje provrad till matrisen
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

        # 4. NYTT: Skapa rubrikraden för summeringsfilen (dynamiskt uppbyggd)
        printf "Sample name\tTotal_Syn" > summary_file
        for (k = 1; k <= 12; k++) printf "\tSyn_%s", sub_types[k] > summary_file
        
        printf "\tTotal_Mis" > summary_file
        for (k = 1; k <= 12; k++) printf "\tMis_%s", sub_types[k] > summary_file
        
        printf "\n" > summary_file

        # 5. NYTT: Skriv ut data-raderna till summeringsfilen
        for (i = 1; i <= n_samples; i++) {
            s = sorted_samples[i]
            
            # Skriv ut Totala Synonyma
            printf "%s\t%d", s, syn_count[s] + 0 > summary_file
            
            # Skriv ut de 12 specifika Synonyma substitutionerna
            for (k = 1; k <= 12; k++) {
                mut = sub_types[k]
                printf "\t%d", syn_sub_count[s, mut] + 0 > summary_file
            }
            
            # Skriv ut Totala Missense
            printf "\t%d", mis_count[s] + 0 > summary_file
            
            # Skriv ut de 12 specifika Missense substitutionerna
            for (k = 1; k <= 12; k++) {
                mut = sub_types[k]
                printf "\t%d", mis_sub_count[s, mut] + 0 > summary_file
            }
            
            # Ny rad för nästa prov
            printf "\n" > summary_file
        }
    }
    ' $vcfs > "$out_file"

    echo "  -> Sparade tabell: $out_file"
    echo "  -> Sparade detaljerad summering: $summary_file"
done

echo "=== Alla matriser och summeringar är färdiga! ==="
