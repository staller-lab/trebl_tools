#!/bin/bash
#SBATCH --job-name=trebl_quick_start
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/quick_start.out
# TREBL Quick Start Analysis Job
#
# Before submitting:
#   1. Edit run_quick_start.py — update DATA_DIR, barcode sequences,
#      file paths, and any other settings in the CONFIGURATION section.
#   2. Create required directories:
#        mkdir -p output/quick_start db logs
#   3. Change time as needed, then submit this script:
#        sbatch examples/savio_jobs/quick_start_job.sh
#
# To run with the bundled example data, submit from the repo root without
# any changes — run_quick_start.py already points to examples/data/.

# Set number of threads for multiprocessing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the analysis script (edit that file to configure your paths/barcodes)
/global/scratch/projects/fc_mvslab/conda/trebl_tools/bin/python "run_quick_start.py"

echo "Job completed at: $(date)"