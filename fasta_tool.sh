#!/bin/bash
# Safety: exit on error or undefined variables
set -ueo pipefail

Help() {
   echo "################################################################################"
   echo "# fasta_tool.sh: Robust FASTA splitting and extraction                         #"
   echo "#                                                                              #"
   echo "# Syntax: fasta_tool.sh [-s | -e] <input_file> <number> [prefix]               #"
   echo "# options:                                                                     #"
   echo "# -s     Split file into <number> of roughly equal parts (records).            #"
   echo "# -e     Extract first <number> of records from the file.                      #"
   echo "# -h     Print this Help.                                                      #"
   echo "################################################################################"
}

split_file() {
    local input="$1"
    local num_splits="$2"
    local prefix="${3:-split_fa}"

    if [[ "$num_splits" -lt 1 ]]; then echo "Error: Splits must be > 0"; exit 1; fi

    echo "📊 Counting records in $input..."
    local total_recs=$(grep -c '>' "$input")
    local recs_per_file=$(( (total_recs + num_splits - 1) / num_splits ))

    echo "✂️ Splitting $total_recs records into $num_splits files (~$recs_per_file records each)..."

    # Use awk for record-aware splitting (prevents broken records)
    awk -v recs_per_file="$recs_per_file" -v prefix="$prefix" '
        BEGIN {RS=">"; ORS=""}
        NR > 1 {
            if ((NR-2) % recs_per_file == 0) {
                file_idx++;
                out_file = sprintf("%s_%02d.fasta", prefix, file_idx);
            }
            print ">"$0 > out_file
        }
    ' "$input"

    echo "✅ Done! Created $(ls ${prefix}_*.fasta | wc -l) split files."
}

extract_seqs() {
    local input="$1"
    local num_recs="$2"
    local prefix="${3:-extracted}"

    echo "📦 Extracting first $num_recs records from $input..."

    # Use awk to extract records safely
    awk -v limit="$num_recs" -v prefix="$prefix" '
        BEGIN {RS=">"; ORS=""; count=0}
        NR > 1 && count < limit {
            print ">"$0 > (prefix ".fasta");
            count++;
        }
        NR > 1 && count >= limit {
            print ">"$0 > (prefix ".remainder.fasta");
        }
    ' "$input"

    echo "✅ Extraction completed: ${prefix}.fasta"
    echo "📝 Remaining sequences saved to: ${prefix}.remainder.fasta"
}

# --- Argument Parsing ---
do_splits=false
do_extract=false

while getopts "seh" opt; do
   case ${opt} in
      s) do_splits=true ;;
      e) do_extract=true ;;
      h) Help; exit 0 ;;
      *) Help; exit 1 ;;
   esac
done
shift $((OPTIND - 1))

# Validation
if [[ $# -lt 2 ]]; then
    echo "❌ Error: Missing arguments."
    Help; exit 1
fi

input_file="$1"
value="$2"
custom_prefix="${3:-}"

if [[ ! -f "$input_file" ]]; then echo "❌ Error: File $input_file not found."; exit 1; fi
if ! [[ "$value" =~ ^[0-9]+$ ]]; then echo "❌ Error: $value is not a valid number."; exit 1; fi

# Execution
if [[ "$do_splits" == true ]]; then
    split_file "$input_file" "$value" "$custom_prefix"
elif [[ "$do_extract" == true ]]; then
    extract_seqs "$input_file" "$value" "$custom_prefix"
else
    echo "❌ Error: You must specify -s (split) or -e (extract)."
    Help; exit 1
fi