#!/bin/bash
#SBATCH --job-name=trebl_quick_start
#SBATCH --account=fc_mvslab
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/quick_start_%j.out
#SBATCH --error=logs/quick_start_%j.err

# TREBL Quick Start Analysis Job
# This script runs a TREBL analysis with:
# - No error correction (faster)
# - Simple UMI deduplication only
#
# Usage:
#   sbatch quick_start_job.sh
#
# Before running:
# 1. Update paths in the Python script section below
# 2. Ensure output and logs directories exist:
#    mkdir -p output/quick_start logs

# Load required modules and activate conda environment
module load python
source activate trebl_env  # Or: conda activate trebl_env

# Set number of threads for multiprocessing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "=========================================="
echo "TREBL Quick Start Analysis"
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
STEP2_AD_SEQ_FILE = "/path/to/your/step2_AD_file.fastq"
STEP2_RT_SEQ_FILE = "/path/to/your/step2_RT_file.fastq"
AD_SEQ_FILES_PATTERN = "/path/to/AD_assembled/*"
RT_SEQ_FILES_PATTERN = "/path/to/RT_assembled/*"
OUTPUT_DIR = "output/quick_start"

# ==========================================
# Initialize Pipeline
# ==========================================
print("\n[1/6] Initializing pipeline...")
pipeline = pipelines.TreblPipeline(
    db_path="quick_start.db",
    design_file_path=DESIGN_FILE,
    error_correction=False,  # Quick start: no error correction
    output_path=OUTPUT_DIR
)

# ==========================================
# Define Barcodes
# ==========================================
print("\n[2/6] Defining barcodes...")
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
# Step 1: TREBL Mapping
# ==========================================
print("\n[3/6] Running Step 1 mapping...")

# Plot reads distribution
print("  - Plotting reads distribution...")
pipeline.step1_reads_distribution(
    seq_file=STEP1_SEQ_FILE,
    bc_objects=bc_objects,
    reverse_complement=True
)

# Run Step 1 mapping
print("  - Running Step 1 mapping...")
step1_map = pipeline.run_step_1(
    seq_file=STEP1_SEQ_FILE,
    bc_objects=bc_objects,
    column_pairs=[("RT_BC", "AD")],
    reads_threshold=10,
    reverse_complement=False
)
print(f"  - Step 1 complete: {len(step1_map)} entries")

# ==========================================
# Step 2: TREBL Step 2 Mapping
# ==========================================
print("\n[4/6] Running Step 2 mapping...")

AD_bc_objects = [AD, AD_BC]
RT_bc_objects = [RT_BC]

# Plot Step 2 reads distribution
print("  - Plotting Step 2 reads distribution...")
pipeline.step2_reads_distribution(
    AD_seq_file=STEP2_AD_SEQ_FILE,
    AD_bc_objects=AD_bc_objects,
    RT_seq_file=STEP2_RT_SEQ_FILE,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True
)

# Run Step 2 mapping
print("  - Running Step 2 mapping...")
step2 = pipeline.run_step_2(
    AD_seq_file=STEP2_AD_SEQ_FILE,
    AD_bc_objects=AD_bc_objects,
    RT_seq_file=STEP2_RT_SEQ_FILE,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    reads_threshold_AD=10,
    reads_threshold_RT=10,
    step1_map_csv_path=f"{OUTPUT_DIR}/step1_AD_AD_BC_RT_BC_designed.csv"
)

AD_step2 = step2["AD_step2"]
RT_step2 = step2["RT_step2"]
step1_overlap = step2["step1_overlap"]

print(f"  - AD Step 2: {len(AD_step2)} entries")
print(f"  - RT Step 2: {len(RT_step2)} entries")

# ==========================================
# TREBL Experiment Analysis
# ==========================================
print("\n[5/6] Running TREBL experiment analysis...")

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

# Run TREBL experiment with simple UMI deduplication
print("  - Running TREBL experiment with simple UMI deduplication...")
trebl_results = pipeline.trebl_experiment_analysis(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    step1_map_csv_path=f"{OUTPUT_DIR}/step1_AD_AD_BC_RT_BC_designed.csv",
    AD_umi_object=AD_UMI,
    RT_umi_object=RT_UMI,
    umi_deduplication='simple'  # Quick start: simple deduplication only
)

AD_results = trebl_results["AD_results"]
RT_results = trebl_results["RT_results"]

print(f"  - AD results: {len(AD_results)} entries")
print(f"  - RT results: {len(RT_results)} entries")

# ==========================================
# Complete
# ==========================================
print("\n[6/6] Analysis complete!")
print(f"  - Results saved to: {OUTPUT_DIR}")
print(f"  - Database file: quick_start.db")

EOF

echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="
