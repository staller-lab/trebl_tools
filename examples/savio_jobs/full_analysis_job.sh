#!/bin/bash
#SBATCH --job-name=trebl_full_analysis
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1:00:00
#SBATCH --output=logs/full_analysis.out
# TREBL Full Analysis Job
#
# Before submitting:
#   1. Edit run_full_analysis.py — update DATA_DIR, barcode sequences,
#      file paths, and any other settings in the CONFIGURATION section.
#   2. Create required directories:
#        mkdir -p ../output/full_analysis db logs
#   3. Change time as needed, then submit this script:
#        sbatch examples/savio_jobs/full_analysis_job.sh
#
# To run with the bundled example data, submit from the repo root without
# any changes — run_full_analysis.py already points to examples/data/.
#
# Note: This job requires more time and resources than quick_start_job.sh
#       due to error correction and complex UMI deduplication.

# Set number of threads for multiprocessing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the analysis script (edit that file to configure your paths/barcodes)
#/global/scratch/projects/fc_mvslab/conda/trebl_tools/bin/python "run_full_analysis.py"
/global/scratch/projects/fc_mvslab/OpenProjects/Sanjana/conda/trebl_env/bin/python "run_full_analysis.py"

echo "Job completed at: $(date)"
