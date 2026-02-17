#!/bin/bash
#SBATCH --job-name=fastp_array
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --cpus-per-task=32
#SBATCH --time=6:00:00

source activate /global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/umi_tools

INPUT_DIR=$1 # First input: input dir
OUTPUT_DIR=$2 # Second output: output dir

# mkdir -p "${OUTPUT_DIR}"
# mkdir -p "${OUTPUT_DIR}/logs"

# Grab all .fastq.gz or .fastq files into a Bash array
FILES=(
  ${INPUT_DIR}/*.fastq
  ${INPUT_DIR}/*.fastq.gz
)

# Get the file for this array task
FILE="${FILES[$SLURM_ARRAY_TASK_ID-1]}"

BASENAME=$(basename "$FILE")
BASENAME=${BASENAME%.gz}
BASENAME=${BASENAME%.fastq}

OUTFILE="${OUTPUT_DIR}/${BASENAME}_fastp.fastq.gz"
LOGFILE="${OUTPUT_DIR}/logs/${BASENAME}_fastp_report.out"
REPORT_JSON="${OUTPUT_DIR}/logs/${BASENAME}_fastp_report.json"
REPORT_HTML="${OUTPUT_DIR}/logs/${BASENAME}_fastp_report.html"

echo "[$(date)] Processing ${FILE} → ${OUTFILE}"

fastp -i "${FILE}" -o "${OUTFILE}" -w ${SLURM_CPUS_PER_TASK} --disable_adapter_trimming --json "${REPORT_JSON}" --html "${REPORT_HTML}" &> "${LOGFILE}"

# echo "fastp -i \"${FILE}\" -o \"${OUTFILE}\" -w ${SLURM_CPUS_PER_TASK} --disable_adapter_trimming --dont_overwrite"

echo "[$(date)] Finished ${FILE}"
