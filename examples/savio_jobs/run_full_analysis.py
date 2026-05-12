"""
TREBL Full Analysis Script
===========================
Edit the CONFIGURATION section below before submitting the Savio job.

What to update before running:
  1. DATA_DIR         - folder containing your FASTQ files and design file
  2. Barcode objects  - update preceder/post/length to match your experiment
  3. UMI objects      - update preceder/length to match your experiment
  4. File path variables (DESIGN_FILE, STEP1_SEQ_FILE, etc.)
  5. OUTPUT_DIR / db_path as needed

Run via:
  sbatch full_analysis_job.sh
Or directly (outside Savio):
  python run_full_analysis.py
"""

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
    pipelines,
)

# ==========================================
# CONFIGURATION — EDIT THIS SECTION
# ==========================================

# Path to the folder containing your FASTQ files and design file.
# The default points to the example data bundled with the repo.
# Change to an absolute path for your own data, e.g.:
#   DATA_DIR = "/global/scratch/projects/fc_mvslab/data/my_experiment"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "../data")

# Input files — update filenames to match your sequencing data
DESIGN_FILE = f"{DATA_DIR}/design_file.txt"
STEP1_SEQ_FILE = f"{DATA_DIR}/step1_ChopTFs_sample.fastq"
STEP2_AD_SEQ_FILE = f"{DATA_DIR}/step2_ChopTFs_AD_sample.fastq"
STEP2_RT_SEQ_FILE = f"{DATA_DIR}/step2_ChopTFs_RT_sample.fastq"

# TREBL experiment sequencing files (one entry per time point / replicate)
# File names must include _tXX and _rX so activity scores can be parsed, e.g.:
#   myexperiment_AD_t10_r1.fastq  →  timepoint 10, replicate 1
AD_SEQ_FILES = [
    f"{DATA_DIR}/trebl_experiment_ChopTFs_AD_t10_r2.fastq",
    f"{DATA_DIR}/trebl_experiment_ChopTFs_AD_t60_r2.fastq",
]
RT_SEQ_FILES = [
    f"{DATA_DIR}/trebl_experiment_ChopTFs_RT_t10_r2.fastq",
    f"{DATA_DIR}/trebl_experiment_ChopTFs_RT_t60_r2.fastq",
]

# Where to save results and the DuckDB working database
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../output/full_analysis")
DB_PATH = os.path.join(SCRIPT_DIR, "../db/full_analysis.db")

# ------------------------------------------
# Barcode definitions — update for your experiment
# Each Barcode has:
#   name     : identifier used in output columns
#   preceder : sequence immediately BEFORE the barcode in the read
#   post     : sequence immediately AFTER the barcode in the read
#   length   : expected barcode/AD length in bases
# ------------------------------------------
AD = finder.Barcode(
    name="AD",
    preceder="GGCTAGC",   # ← update if your flanking sequences differ
    post="TGACTAG",
    length=120,
)

AD_BC = finder.Barcode(
    name="AD_BC",
    preceder="CGCGCC",    # ← update if your flanking sequences differ
    post="GGGCCC",
    length=11,
)

RT_BC = finder.Barcode(
    name="RT_BC",
    preceder="CTCGAG",    # ← update if your flanking sequences differ
    post="GGCCGC",
    length=14,
)

# UMI definitions — update preceder/length for your experiment
AD_UMI = finder.Barcode(
    name="UMI",
    preceder="TGATTT",    # ← update if your UMI flanking sequence differs
    post="",
    length=12,
)

RT_UMI = finder.Barcode(
    name="UMI",
    preceder="TGTCAC",    # ← update if your UMI flanking sequence differs
    post="",
    length=12,
)

# Minimum reads to keep a barcode (adjust based on your reads distribution)
READS_THRESHOLD = 1
READS_THRESHOLD_AD = 1
READS_THRESHOLD_RT = 1

# ==========================================
# (No edits needed below this line)
# ==========================================

# Derived barcode groups
bc_objects = [AD, AD_BC, RT_BC]
AD_bc_objects = [AD, AD_BC]
RT_bc_objects = [RT_BC]

# ==========================================
# Initialize Pipeline
# ==========================================
print("\n[1/6] Initializing pipeline...")
pipeline = pipelines.TreblPipeline(
    db_path=DB_PATH,
    design_file_path=DESIGN_FILE,
    error_correction=True,  # Full analysis: enable error correction
    output_path=OUTPUT_DIR,
    # test_n_reads=100000  # Uncomment to test with first 100k reads
)

# ==========================================
# Step 1: TREBL Mapping with Error Correction
# ==========================================
print("\n[2/6] Running Step 1 mapping with error correction...")

print("  - Plotting reads distribution...")
pipeline.step1_reads_distribution(
    seq_file=STEP1_SEQ_FILE,
    bc_objects=bc_objects,
    reverse_complement=True,
)

print("  - Running Step 1 mapping with error correction...")
step1_map = pipeline.run_step_1(
    seq_file=STEP1_SEQ_FILE,
    bc_objects=bc_objects,
    column_pairs=[("RT_BC", "AD")],
    reads_threshold=READS_THRESHOLD,
    reverse_complement=False,
)
print(f"  - Step 1 complete: {len(step1_map)} entries")

# ==========================================
# Step 2: TREBL Step 2 Mapping with Error Correction
# ==========================================
print("\n[3/6] Running Step 2 mapping with error correction...")

print("  - Plotting Step 2 reads distribution...")
pipeline.step2_reads_distribution(
    AD_seq_file=STEP2_AD_SEQ_FILE,
    AD_bc_objects=AD_bc_objects,
    RT_seq_file=STEP2_RT_SEQ_FILE,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
)

print("  - Running Step 2 mapping with error correction...")
step2 = pipeline.run_step_2(
    AD_seq_file=STEP2_AD_SEQ_FILE,
    AD_bc_objects=AD_bc_objects,
    RT_seq_file=STEP2_RT_SEQ_FILE,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    reads_threshold_AD=READS_THRESHOLD_AD,
    reads_threshold_RT=READS_THRESHOLD_RT,
    step1_map_csv_path=f"{OUTPUT_DIR}/step1.csv",
)

AD_step2 = step2["AD_step2"]
RT_step2 = step2["RT_step2"]
step1_overlap = step2["step1_overlap"]

print(f"  - AD Step 2: {len(AD_step2)} entries")
print(f"  - RT Step 2: {len(RT_step2)} entries")

# ==========================================
# TREBL Experiment Analysis
# ==========================================
print("\n[4/6] Running TREBL experiment analysis with both UMI deduplication methods...")

print(f"  - Using {len(AD_SEQ_FILES)} AD files and {len(RT_SEQ_FILES)} RT files")

print("  - Running TREBL experiment with both simple and directional UMI deduplication...")
print("    (This may take significant time for large datasets)")
trebl_results = pipeline.trebl_experiment_analysis(
    AD_seq_files=AD_SEQ_FILES,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=RT_SEQ_FILES,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    step1_map_csv_path=f"{OUTPUT_DIR}/step1.csv",
    AD_umi_object=AD_UMI,
    RT_umi_object=RT_UMI,
    umi_deduplication="both",  # Full analysis: both simple and directional deduplication
)

print(f"AD results and RT results ready.")

print("\n[5/6] Calculating activity scores...")
ad_activity_per_barcode_df = pipeline.calculate_activity_scores(
    step1_path=f"{OUTPUT_DIR}/step1.csv",
    AD_bc_objects=AD_bc_objects,
    RT_bc_objects=RT_bc_objects,
    time_regex=r"_t(\d+)",
    rep_regex=r"_r(\d+)",
)
print(f"  - Activity rows (per barcode): {len(ad_activity_per_barcode_df)}")

print("\n[6/6] Activity score tables written by calculate_activity_scores().")

print("\nAnalysis complete!")
