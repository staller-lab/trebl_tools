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
# This script runs a comprehensive TREBL analysis with:
# - Error correction enabled (improves accuracy)
# - Both simple and directional/complex UMI deduplication
#
# Usage:
#   sbatch full_analysis_job.sh
#
# Before running:
# 1. Update paths in the Python script section below
# 2. Ensure output and logs directories exist:
#    mkdir -p output/full_analysis logs
#
# Note: This job requires more time and resources than quick_start_job.sh
#       due to error correction and complex UMI deduplication.

# Load required modules and activate conda environment
module load python
source activate trebl_env  # Or: conda activate trebl_env

# Set number of threads for multiprocessing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "=========================================="
echo "TREBL Full Analysis"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Running on node: $(hostname)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "=========================================="

# Run the analysis
python << 'EOF'
import sys
import os
import glob

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for cluster
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
from tqdm import tqdm

from trebl_tools import (
    initial_map,
    map_refiner,
    complexity,
    finder,
    preprocess,
    error_correct,
    plotting,
    umi_deduplicate,
    pipelines
)

# ==========================================
# CONFIGURATION - UPDATE THESE PATHS
# ==========================================
DESIGN_FILE = "/path/to/your/design_file.txt"
STEP1_SEQ_FILE = "/path/to/your/step1_sequencing_file.fastq"
AD_SEQ_FILES_PATTERN = "/path/to/AD_assembled/*"
RT_SEQ_FILES_PATTERN = "/path/to/RT_assembled/*"
OUTPUT_DIR = "output/full_analysis"

# ==========================================
# Initialize Pipeline
# ==========================================
print("\n[1/5] Initializing pipeline...")
pipeline = pipelines.TreblPipeline(
    db_path="full_analysis.db",
    design_file_path=DESIGN_FILE,
    error_correction=True,  # Full analysis: enable error correction
    output_path=OUTPUT_DIR
)

# ==========================================
# Define Barcodes
# ==========================================
print("\n[2/5] Defining barcodes...")
AD = finder.Barcode(
    name="AD",
    preceder="GGCTAGC",
    post="TGACTAG",
    length=120
)

AD_BC = finder.Barcode(
    name="AD_BC",
    preceder="CGCGCC",
    post="GGGCCC",
    length=11
)

RT_BC = finder.Barcode(
    name="RT_BC",
    preceder="CTCGAG",
    post="GGCCGC",
    length=14
)

bc_objects = [AD, AD_BC, RT_BC]

# Define UMIs
AD_UMI = finder.Barcode(
    name="UMI",
    preceder="TGATTT",
    post="",
    length=12
)

RT_UMI = finder.Barcode(
    name="UMI",
    preceder="TGTCAC",
    post="",
    length=12
)

# ==========================================
# Step 1: TREBL Mapping with Error Correction
# ==========================================
print("\n[3/5] Running Step 1 mapping with error correction...")

# Plot reads distribution
print("  - Plotting reads distribution...")
pipeline.step1_reads_distribution(
    seq_file=STEP1_SEQ_FILE,
    bc_objects=bc_objects,
    reverse_complement=True
)

# Run Step 1 mapping with error correction
print("  - Running Step 1 mapping with error correction...")
step1_map = pipeline.run_step_1(
    seq_file=STEP1_SEQ_FILE,
    bc_objects=bc_objects,
    column_pairs=[("RT_BC", "AD")],
    reads_threshold=10,
    reverse_complement=False
)
print(f"  - Step 1 complete: {len(step1_map)} entries")

# ==========================================
# TREBL Experiment Analysis
# ==========================================
print("\n[4/5] Running TREBL experiment analysis with both UMI deduplication methods...")

# Collect sequencing files
trebl_AD_seq_files = glob.glob(AD_SEQ_FILES_PATTERN)
trebl_RT_seq_files = glob.glob(RT_SEQ_FILES_PATTERN)

print(f"  - Found {len(trebl_AD_seq_files)} AD files")
print(f"  - Found {len(trebl_RT_seq_files)} RT files")

AD_bc_objects = [AD, AD_BC]
RT_bc_objects = [RT_BC]

# Plot reads distribution
print("  - Plotting TREBL experiment reads distribution...")
pipeline.trebl_experiment_reads_distribution(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True
)

# Run TREBL experiment with both simple and directional/complex UMI deduplication
print("  - Running TREBL experiment with both simple and directional UMI deduplication...")
print("    (This may take significant time for large datasets)")
trebl_results = pipeline.trebl_experiment_analysis(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    step1_map_csv_path=f"{OUTPUT_DIR}/step1_AD_AD_BC_RT_BC_error_corrected_designed.csv",
    AD_umi_object=AD_UMI,
    RT_umi_object=RT_UMI,
    umi_deduplication='both'  # Full analysis: both simple and directional deduplication
)

AD_results = trebl_results["AD_results"]
RT_results = trebl_results["RT_results"]

print(f"  - AD results: {len(AD_results)} entries")
print(f"  - RT results: {len(RT_results)} entries")

# ==========================================
# Summary Statistics
# ==========================================
print("\n[5/5] Generating summary...")
print("\nAnalysis Summary:")
print("-" * 50)
print(f"Error correction: Enabled")
print(f"UMI deduplication: Both simple and directional")
print(f"Step 1 map entries: {len(step1_map)}")
print(f"AD results entries: {len(AD_results)}")
print(f"RT results entries: {len(RT_results)}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Database file: full_analysis.db")
print("-" * 50)

print("\nAnalysis complete!")

EOF

echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="
