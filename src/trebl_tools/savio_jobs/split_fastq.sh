#!/bin/bash

# Usage: ./split_fastq.sh -i <input_fastq> -o <output_dir> -c <chunk1,chunk2,...> [-j <parallel_jobs>]
# Example: ./split_fastq.sh -i sample.fastq.gz -o ./chunks -c 50,100,200 -j 24

usage() {
    echo "Usage: $0 -i <input_fastq> -o <output_dir> -c <chunks> [-j <parallel_jobs>]"
    echo ""
    echo "  -i    Path to input FASTQ file (required)"
    echo "  -o    Path to output directory (required)"
    echo "  -c    Comma-separated list of chunk sizes (required), e.g. 50,100,200"
    echo "  -j    Number of parallel jobs (optional, default: 24)"
    echo ""
    echo "Example:"
    echo "  $0 -i /path/to/sample.fastq.gz -o /path/to/output -c 50,100,200 -j 8"
    exit 1
}

# Defaults
PARALLEL_JOBS=24

# Parse arguments
while getopts "i:o:c:j:h" opt; do
    case $opt in
        i) INPUT_FASTQ="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) IFS=',' read -ra CHUNKS <<< "$OPTARG" ;;
        j) PARALLEL_JOBS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_FASTQ" || -z "$OUTPUT_DIR" || -z "${CHUNKS[*]}" ]]; then
    echo "Error: -i, -o, and -c are required arguments."
    usage
fi

# Validate input file exists
if [[ ! -f "$INPUT_FASTQ" ]]; then
    echo "Error: Input file not found: $INPUT_FASTQ"
    exit 1
fi

# Derive base name from input file (strip path and extensions like .fastq.gz / .fq.gz)
BASENAME=$(basename "$INPUT_FASTQ")
BASENAME="${BASENAME%.fastq.gz}"
BASENAME="${BASENAME%.fq.gz}"
BASENAME="${BASENAME%.fastq}"
BASENAME="${BASENAME%.fq}"

mkdir -p "$OUTPUT_DIR"

# Function to split for one N
split_fastq() {
    local N=$1
    local INPUT_FASTQ=$2
    local OUTPUT_DIR=$3
    local BASENAME=$4

    local ARGS=()
    for i in $(seq 1 $N); do
        ARGS+=("-o" "$OUTPUT_DIR/${BASENAME}_${N}_chunks_part_${i}.fq.gz")
    done

    echo "Splitting into $N chunks..."
    fastqsplitter -i "$INPUT_FASTQ" "${ARGS[@]}"
}

export -f split_fastq

# Run all chunk splits in parallel
parallel -j "$PARALLEL_JOBS" split_fastq {} "$INPUT_FASTQ" "$OUTPUT_DIR" "$BASENAME" ::: "${CHUNKS[@]}"

echo "All splits completed!"