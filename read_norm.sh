#!/bin/bash
set -ueo pipefail

# ==============================================================================
# Standalone BBTools Normalization for V-pipe (Single-End Workaround)
# ==============================================================================

# Standardvärden
INPUT_DIR=""
OUTPUT_DIR=""
NORM_TARGET=5000  # Högre default för V-pipe (kvasispecies)
NORM_MIN=2        # Sänkt till 2 för att bevara lågtäckta regioner (dalar)
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

# Kontrollera att nödvändiga BBTools-verktyg är tillgängliga
if ! command -v bbnorm.sh &> /dev/null; then
    echo "❌ Fel: 'bbnorm.sh' hittades inte."
    exit 1
fi

if ! command -v repair.sh &> /dev/null; then
    echo "❌ Fel: 'repair.sh' hittades inte. Krävs för att synkronisera R1 och R2."
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
    out1="${out_prefix}_R1.fastq.gz"
    out2="${out_prefix}_R2.fastq.gz"
    
    # Separata loggfiler för de två stegen
    log_file1="${LOG_DIR}/${sample_name}_bbnorm_R1.log"
    log_file2="${LOG_DIR}/${sample_name}_repair_R2.log"

    # Kontrollera om filen redan är processad och matchande R2 finns
    if [[ ! -f "$f2" ]]; then
        echo "⚠️ Varning: Hittade ingen matchande R2-fil för $f1. Hoppar över..."
        continue
    fi

    if [[ ! -f "$out1" || ! -f "$out2" ]]; then
        echo "   Processar: $sample_name"
        
        # --- STEG 1: Normalisera enbart R1 (Bryter ankareffekten) ---
        echo "     -> Steg 1/2: Plattar till täckningen för Read 1..."
        bbnorm.sh \
            in="$f1" \
            out="$out1" \
            target="$NORM_TARGET" \
            min="$NORM_MIN" \
            minkmers=0 \
            ecc=f \
	    passes=1 \
	    percentile=54 \
            bits="$NORM_BITS" \
            unpigz=t \
            &> "$log_file1"
            
        if [ $? -ne 0 ]; then
             echo "❌ Fel vid normalisering av $sample_name (Steg 1). Se logg: $log_file1"
             continue
        fi

        # --- STEG 2: Återställ paret för R2 med repair.sh ---
        echo "     -> Steg 2/2: Återskapar par-synkronisering..."
        
        # Vi lägger out1 i en temporär fil under synkroniseringen
        out1_temp="${out_prefix}_R1_temp.fastq.gz"

        repair.sh \
            in="$out1" \
            in2="$f2" \
            out="$out1_temp" \
            out2="$out2" \
            outs=/dev/null \
            &> "$log_file2"

        if [ $? -ne 0 ]; then
             echo "❌ Fel vid synkronisering av $sample_name (Steg 2). Se logg: $log_file2"
        else
             # Skriv över den initiala R1 med den 100% synkroniserade R1-filen
             mv "$out1_temp" "$out1"
        fi
        
    else
        echo "   Hoppar över $sample_name (både R1 och R2 är redan genererade)"
    fi
done

echo "✅ Normalisering och par-synkronisering klar. Filer sparade i: $OUTPUT_DIR"
