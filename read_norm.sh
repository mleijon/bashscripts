#!/bin/bash
set -ueo pipefail

# ==============================================================================
# Standalone BBTools Normalization for V-pipe
# ==============================================================================

# Standardvärden
INPUT_DIR=""
OUTPUT_DIR=""
NORM_TARGET=5000  # Högre default för V-pipe (kvasispecies)
NORM_MIN=5
NORM_BITS=32

# --- Hjälpfunktion ---
Help() {
    echo "Användning: $0 -i <input_dir> -o <output_dir> [alternativ]"
    echo ""
    echo "Obligatoriska parametrar:"
    echo "  -i [katalog]  Mapp med input-filer (Måste innehålla '_R1' och '_R2' fastq/fastq.gz-filer)"
    echo "  -o [katalog]  Mapp för normaliserade output-filer"
    echo ""
    echo "Frivilliga parametrar:"
    echo "  -t [heltal]   Måltäckning (Target coverage). Default: $NORM_TARGET (anpassat för V-pipe)"
    echo "  -m [heltal]   Minimitäckning (Kastade k-mers). Default: $NORM_MIN"
    echo "  -b [heltal]   Bloom filter bits. Default: $NORM_BITS"
    echo "  -h            Visa denna hjälp"
}

# --- Argumentparsning ---
while getopts "i:o:t:m:b:h" opt; do
    case ${opt} in
        i) INPUT_DIR=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        t) NORM_TARGET=$OPTARG ;;
        m) NORM_MIN=$OPTARG ;;
        b) NORM_BITS=$OPTARG ;;
        h) Help; exit 0 ;;
        *) Help; exit 1 ;;
    esac
done

# --- Validering av input ---
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "❌ Fel: Både input-katalog (-i) och output-katalog (-o) måste anges."
    Help
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "❌ Fel: Input-katalogen '$INPUT_DIR' existerar inte."
    exit 1
fi

# Kontrollera att bbnorm är tillgängligt (t.ex. via din conda-miljö)
if ! command -v bbnorm.sh &> /dev/null; then
    echo "❌ Fel: 'bbnorm.sh' hittades inte. Är rätt Conda-miljö aktiverad?"
    exit 1
fi

# Skapa output-mapp och logg-mapp om de inte existerar
mkdir -p "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$LOG_DIR"

# --- Utför Normalisering ---
echo "📊 Startar normalisering (Måltäckning: ${NORM_TARGET}x, Min: ${NORM_MIN})"

shopt -s nullglob
# Hanterar både okomprimerade och gzippade fastq-filer
INPUT_FILES=("$INPUT_DIR"/*_R1*fastq*)
shopt -u nullglob

if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "⚠️ Inga filer matchande '*_R1*fastq*' hittades i $INPUT_DIR"
    exit 1
fi

for f1 in "${INPUT_FILES[@]}"; do
    [ -e "$f1" ] || continue
    
    # Skapa R2-filnamnet baserat på R1
    f2="${f1/_R1/_R2}"
    
    # Extrahera filnamnet utan sökväg och R1-suffix för snyggare namngivning
    sample_name=$(basename "$f1" | sed 's/_R1.*//')
    out_prefix="$OUTPUT_DIR/${sample_name}"
    
    # Filnamn för output
    out1="${out_prefix}_R1_norm.fastq.gz"
    out2="${out_prefix}_R2_norm.fastq.gz"
    log_file="${LOG_DIR}/${sample_name}_bbnorm.log"

    # Kontrollera om filen redan är processad och matchande R2 finns
    if [[ ! -f "$f2" ]]; then
        echo "⚠️ Varning: Hittade ingen matchande R2-fil för $f1. Hoppar över..."
        continue
    fi

    if [[ ! -f "$out1" ]]; then
        echo "   Processar: $sample_name"
        
        bbnorm.sh \
            in1="$f1" \
            in2="$f2" \
            out1="$out1" \
            out2="$out2" \
            target="$NORM_TARGET" \
            min="$NORM_MIN" \
            bits="$NORM_BITS" \
            unpigz=t \
            &> "$log_file"
            
        # Enkel felhantering om bbnorm kraschar (t.ex. slut på minne)
        if [ $? -ne 0 ]; then
             echo "❌ Fel vid processering av $sample_name. Se logg: $log_file"
        fi
    else
        echo "   Hoppar över $sample_name (redan normaliserad)"
    fi
done

echo "✅ Normalisering klar. Filer sparade i: $OUTPUT_DIR"
