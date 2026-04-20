#!/bin/bash
#SBATCH --job-name=fastp_array
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --cpus-per-task=32
#SBATCH --time=6:00:00

set -euo pipefail
shopt -s nullglob

INPUT_DIR=$1   # First input: input dir
OUTPUT_DIR=$2  # Second output: output dir

# Make output directories if they don't exist
mkdir -p "${OUTPUT_DIR}/logs"

# Grab all .fastq and .fastq.gz files into a sorted Bash array
mapfile -t FILES < <(ls "${INPUT_DIR}"/*.fastq "${INPUT_DIR}"/*.fastq.gz 2>/dev/null | sort)

# Exit if no files found
if [ ${#FILES[@]} -eq 0 ]; then
    echo "No FASTQ files found in ${INPUT_DIR}"
    exit 1
fi

# Ensure array index is valid
INDEX=$((SLURM_ARRAY_TASK_ID-1))
if [ "$INDEX" -ge "${#FILES[@]}" ]; then
    echo "SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} exceeds file count ${#FILES[@]}"
    exit 1
fi

FILE="${FILES[$INDEX]}"

# Get basename without extensions
BASENAME=$(basename "$FILE")
BASENAME=${BASENAME%.gz}
BASENAME=${BASENAME%.fastq}

OUTFILE="${OUTPUT_DIR}/${BASENAME}_fastp.fastq.gz"
LOGFILE="${OUTPUT_DIR}/logs/${BASENAME}_fastp_report.out"
REPORT_JSON="${OUTPUT_DIR}/logs/${BASENAME}_fastp_report.json"
REPORT_HTML="${OUTPUT_DIR}/logs/${BASENAME}_fastp_report.html"

echo "[$(date)] Processing ${FILE} → ${OUTFILE}"

fastp \
  -i "${FILE}" \
  -o "${OUTFILE}" \
  -w ${SLURM_CPUS_PER_TASK} \
  --disable_adapter_trimming \
  --json "${REPORT_JSON}" \
  --html "${REPORT_HTML}" \
  &> "${LOGFILE}"

echo "[$(date)] Finished ${FILE}"