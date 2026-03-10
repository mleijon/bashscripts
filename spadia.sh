#!/bin/bash
# Safety: exit on error, undefined variables, or pipe failures
set -ueo pipefail

# ==============================================================================
# SPADIA PIPELINE V4 (Optimized, Modular & Environment-Aware)
# Trimming -> Assembly (SPAdes/USEARCH) -> Classification (Diamond)
# ==============================================================================

# --- CONDA ENVIRONMENT CHECK ---
# Ensures the 'spadia' environment is active before proceeding
if [[ "${CONDA_DEFAULT_ENV:-}" != "spadia" ]]; then
    echo "❌ Error: The 'spadia' conda environment is not active."
    echo "Please activate it using: conda activate spadia"
    exit 1
fi

# --- DEFAULT CONFIGURATION ---
# Databases and Tool Paths
DIAMOND_DB="/ssd2/classify/nr"
USEARCH_BIN="/ssd2/classify/usearch11.0.667_i86linux64"

# Resources
THREADS=$(grep -c 'processor' /proc/cpuinfo)
# Specific thread cap for Trimmomatic to avoid Java memory issues
TRIM_THREADS=$(( THREADS > 8 ? 8 : THREADS ))
MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
MEM_GB=$((MEM_KB / 1024 / 1024))
SPADES_MEM=$((MEM_GB * 8 / 10))  # Default: 80% of system RAM

# Default Flags
SPADES_FLAG=""       # Standard mode
SEQ_TYPE="c"         # 'c' for contigs, 'r' for reads
RAW_MODE="n"         # Skip trimming if 'y'
REMOVE_TEMP="y"      # Remove intermediate files
DIAMOND_BLOCK=12     # Diamond block size
DIAMOND_CHUNKS=1

# --- HELP FUNCTION ---
Help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -s [c|r]    Sequence type: contigs (c) or reads (r). Default: c"
    echo "  -x [y|n]    Raw reads mode (skip trimming/de-rep). Default: n"
    echo "  -f [flag]   SPAdes mode (meta, plasmid, etc.). Default: standard"
    echo "  -m [0.1-1]  Fraction of RAM for SPAdes. Default: 0.8"
    echo "  -p [int]    Number of threads. Default: $THREADS"
    echo "  -r [y|n]    Remove intermediate files. Default: y"
    echo "  -h          Show this help"
}

# --- ARGUMENT PARSING ---
while getopts "m:p:b:c:f:s:r:x:h" opt; do
    case ${opt} in
        m) SPADES_MEM=$(awk "BEGIN {print int($MEM_GB * $OPTARG)}") ;;
        p) THREADS=$OPTARG ;;
        f) SPADES_FLAG=$OPTARG ;;
        b) DIAMOND_BLOCK=$OPTARG ;;
        c) DIAMOND_CHUNKS=$OPTARG ;;
        s) SEQ_TYPE=$OPTARG ;;
        r) REMOVE_TEMP=$OPTARG ;;
        x) RAW_MODE=$OPTARG ;;
        h) Help; exit 0 ;;
        *) Help; exit 1 ;;
    esac
done

# --- DIRECTORY SETUP ---
TRIM_DIR="./trimmed"
ASM_DIR="./${SPADES_FLAG:-standard}_spades"
MERGED_DIR="./merged"
DMND_DIR="./diamond"

mkdir -p "$TRIM_DIR" "$DMND_DIR"

# ==============================================================================
# 1. TRIMMING (TRIMMOMATIC)
# ==============================================================================
if [[ "$RAW_MODE" == "n" ]]; then
    echo "✂️ Starting Trimmomatic (Threads capped at