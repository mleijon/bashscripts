#!/bin/bash

# Skapa en mapp för de godkända genomen
mkdir -p hq_mq_mags

# Läs CheckM2-rapporten, hitta rader där kolumn 2 (Completeness) är >= 50, och kopiera motsvarande .fna-fil
awk -F'\t' 'NR>1 && $2 >= 50 {print $1}' checkm2_out/quality_report.tsv | while read bin; do
    cp "out_taxvamb/bins/${bin}.fna" hq_mq_mags/
done

echo "Kopierade $(ls hq_mq_mags | wc -l) MAGs till hq_mq_mags/"
