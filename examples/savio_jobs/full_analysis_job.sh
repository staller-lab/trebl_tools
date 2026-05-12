#!/bin/bash
#SBATCH --job-name=trebl_full_analysis
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=logs/full_analysis_%j.out
#SBATCH --error=logs/full_analysis_%j.err

# TREBL Full Analysis Job
#
# Before submitting:
#   1. Edit run_full_analysis.py — update DATA_DIR, barcode sequences,
#      file paths, and any other settings in the CONFIGURATION section.
#   2. Create required directories:
#        mkdir -p output/full_analysis db logs
#   3. Submit this script:
#        sbatch examples/savio_jobs/full_analysis_job.sh
#
# To run with the bundled example data, submit from the repo root without
# any changes — run_full_analysis.py already points to examples/data/.
#
# Note: This job requires more time and resources than quick_start_job.sh
#       due to error correction and complex UMI deduplication.

# Load required modules and activate conda environment
module load python
source activate trebl_tools_env  # Or: conda activate trebl_tools_env

# Set number of threads for multiprocessing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "=========================================="
echo "TREBL Full Analysis"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Running on node: $(hostname)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "=========================================="

# Run the analysis script (edit that file to configure your paths/barcodes)
python "$(dirname "$0")/run_full_analysis.py"

echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="


echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="
