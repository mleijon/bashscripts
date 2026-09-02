#!/bin/bash
set -ueo pipefail

# ==============================================================================
# Shovill + MLST + RagTag Pipeline för Bakteriella Isolat
# ==============================================================================

# Standardvärden för resurser
THREADS=72
RAM=500
INDIR=""
OUTDIR=""
REF=""

Help() {
    echo "Användning: $0 -i <indata_mapp> -o <utdata_mapp> -r <referens.fasta> [alternativ]"
    echo "Alternativ:"
    echo "  -i    Mapp med indata (FASTQ-filer). Måste innehålla '_R1' och '_R2'."
    echo "  -o    Mapp för utdata (skapas automatiskt om den inte finns)."
    echo "  -r    Referenssekvens (FASTA) för RagTag scaffolding."
    echo "  -t    Antal trådar (standard: $THREADS)"
    echo "  -m    RAM i GB (standard: $RAM)"
    echo "  -h    Visa denna hjälp"
}

while getopts "i:o:r:t:m:h" opt; do
    case ${opt} in
        i) INDIR=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        r) REF=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        m) RAM=$OPTARG ;;
        h) Help; exit 0 ;;
        *) Help; exit 1 ;;
    esac
done

# Validering av indata
if [[ -z "$INDIR" || -z "$OUTDIR" || -z "$REF" ]]; then
    echo "❌ Fel: Saknar obligatoriska parametrar."
    Help
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "❌ Fel: Hittar inte referensfilen '$REF'."
    exit 1
fi

mkdir -p "$OUTDIR"

# Hitta alla R1-filer
shopt -s nullglob
R1_FILES=("$INDIR"/*_R1*.fastq*)
shopt -u nullglob

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "❌ Fel: Inga filer som matchar '*_R1*.fastq*' hittades i $INDIR."
    exit 1
fi

echo "🚀 Startar Shovill + MLST + RagTag Pipeline"
echo "Trådar: $THREADS | RAM: ${RAM}GB | Referens: $REF"

for R1 in "${R1_FILES[@]}"; do
    # Extrahera provnamn och hitta motsvarande R2-fil
    FILENAME=$(basename "$R1")
    SAMPLE="${FILENAME/_R1*/}"
    R2="${R1/_R1/_R2}"

    if [[ ! -f "$R2" ]]; then
        echo "⚠️ Varning: Hittade ingen R2-fil för $R1. Hoppar över $SAMPLE."
        continue
    fi

    echo "---------------------------------------------------"
    echo "🦠 Processar prov: $SAMPLE"
    
    SAMPLE_OUT="$OUTDIR/$SAMPLE"
    SHOVILL_OUT="$SAMPLE_OUT/shovill"
    RAGTAG_OUT="$SAMPLE_OUT/ragtag"

    # Skapa provets mapp i förväg så att loggarna har någonstans att ta vägen
    mkdir -p "$SAMPLE_OUT"

    # Steg 1: Shovill
    if [[ ! -f "$SHOVILL_OUT/contigs.fa" ]]; then
        echo "  🧬 Kör Shovill assembly..."
        shovill --R1 "$R1" --R2 "$R2" --outdir "$SHOVILL_OUT" \
                --cpus "$THREADS" --ram "$RAM" --force > "$SAMPLE_OUT/shovill.log" 2>&1
    else
        echo "  ✅ Shovill-utdata finns redan. Hoppar över assembly."
    fi

    # Steg 1.5: Species Confirmation (MLST)
    if [[ -f "$SHOVILL_OUT/contigs.fa" && ! -f "$SAMPLE_OUT/species_typing.tsv" ]]; then
        echo "  🔍 Verifierar art och ST-typ (MLST)..."
        
        # Kör mlst genom dess isolerade miljö för att undvika beroendekonflikter
        conda run -n mlst_env mlst "$SHOVILL_OUT/contigs.fa" > "$SAMPLE_OUT/species_typing.tsv" 2> /dev/null || true
        
        # Skriv ut resultatet direkt i terminalen för snabb överblick
        if [[ -s "$SAMPLE_OUT/species_typing.tsv" ]]; then
            awk -F'\t' '{print "     Resultat: " $2 " (ST: " $3 ")"}' "$SAMPLE_OUT/species_typing.tsv"
        fi
    fi

    # Steg 2: RagTag
    if [[ -f "$SHOVILL_OUT/contigs.fa" && ! -f "$RAGTAG_OUT/ragtag.scaffold.fasta" ]]; then
        echo "  🧩 Kör RagTag scaffolding..."
        ragtag.py scaffold "$REF" "$SHOVILL_OUT/contigs.fa" -t "$THREADS" -o "$RAGTAG_OUT" > "$SAMPLE_OUT/ragtag.log" 2>&1
        
        echo "  🎯 Klar med $SAMPLE! Färdig scaffold hittas i:"
        echo "     $RAGTAG_OUT/ragtag.scaffold.fasta"
    elif [[ -f "$RAGTAG_OUT/ragtag.scaffold.fasta" ]]; then
        echo "  ✅ RagTag-utdata finns redan."
    else
        echo "  ❌ Fel: Kunde inte starta RagTag eftersom Shovill-contigs saknas."
    fi
done

echo "---------------------------------------------------"
echo "🎉 Pipelinen är färdigkörd!"
