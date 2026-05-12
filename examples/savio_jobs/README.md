# Savio Job Scripts

This directory contains SLURM batch scripts and matching Python analysis scripts for running complete TREBL analysis workflows on the Savio cluster.

## File Structure

Each workflow is split into two files — **edit the `.py` file**, then submit the `.sh` file.

| File | Purpose |
|------|---------|
| `run_quick_start.py` | **Edit this** — all paths, barcodes, and parameters for the quick-start workflow |
| `quick_start_job.sh` | SLURM wrapper — submits `run_quick_start.py` as a cluster job |
| `run_full_analysis.py` | **Edit this** — all paths, barcodes, and parameters for the full analysis workflow |
| `full_analysis_job.sh` | SLURM wrapper — submits `run_full_analysis.py` as a cluster job |

Each workflow covers all three analysis steps:
1. **Step 1**: Initial TREBL mapping to establish barcode relationships
2. **Step 2**: Process separated AD and RT libraries
3. **TREBL Experiment**: Full experiment analysis with UMI deduplication and activity scores

---

## When to Use a Savio Job

Use a Savio batch job instead of a Jupyter notebook when:

- **Any sequencing file is large (>10M reads)** — reads distribution plots and mapping steps can take 30+ minutes per file interactively.
- **You are running error correction** — the full analysis workflow is significantly slower; a job guarantees it won't be interrupted.
- **You have many TREBL experiment files** — processing 10+ FASTQ files sequentially benefits from dedicated compute time.
- **You need reproducibility** — job logs capture everything that was run and any errors.

For small datasets or quick exploration, the interactive Jupyter notebooks (`../notebooks/`) are easier to use.

---

## Available Workflows

### 1. Quick Start (`run_quick_start.py` + `quick_start_job.sh`)

Runs a quick TREBL analysis workflow optimized for speed.

**Configuration:**
- **CPUs:** 8 cores
- **Time limit:** 4 hours
- **Features:** No error correction, simple UMI deduplication only

**Best for:**
- Initial data exploration
- Testing your barcode/path configuration
- Time-sensitive analysis

---

### 2. Full Analysis (`run_full_analysis.py` + `full_analysis_job.sh`)

Runs a comprehensive TREBL analysis with maximum accuracy.

**Configuration:**
- **CPUs:** 16 cores
- **Time limit:** 12 hours
- **Features:** Error correction enabled, both simple and directional UMI deduplication

**Best for:**
- Final, publication-quality results
- Variable-quality data
- Comprehensive UMI deduplication comparison

---

## How to Run

### Step 1: Edit the Python analysis file

Open `run_quick_start.py` (or `run_full_analysis.py`) and update the **CONFIGURATION** section at the top:

```python
# Path to your data folder
DATA_DIR = "/global/scratch/projects/fc_mvslab/data/my_experiment"

# Input files
DESIGN_FILE    = f"{DATA_DIR}/design_file.txt"
STEP1_SEQ_FILE = f"{DATA_DIR}/step1.fastq"
# ... etc.

# Barcode flanking sequences — update for your experiment
AD = finder.Barcode(
    name="AD",
    preceder="GGCTAGC",  # ← sequence just BEFORE the AD in your reads
    post="TGACTAG",      # ← sequence just AFTER the AD
    length=120,
)
# ... and the other barcodes / UMIs
```

### Step 2: Create required directories

```bash
mkdir -p output/quick_start db logs
# or for full analysis:
mkdir -p output/full_analysis db logs
```

### Step 3: Submit the job

```bash
sbatch examples/savio_jobs/quick_start_job.sh
# or
sbatch examples/savio_jobs/full_analysis_job.sh
```

To run with the bundled example data without any edits, submit from the repo root — the default `DATA_DIR` already points to `examples/data/`.

### Step 4: Monitor the job

```bash
# Check queue status
squeue -u $USER

# Watch live output
tail -f logs/quick_start_*.out

# Check for errors
tail -f logs/quick_start_*.err
```

### Step 5: Review results

Once the job completes, check:
- CSV files and PNG plots in your `OUTPUT_DIR`
- Log files for any warnings or errors

---

## Key Differences Between Workflows

| Feature | Quick Start | Full Analysis |
|---------|-------------|---------------|
| `error_correction` | `False` | `True` |
| `umi_deduplication` | `'simple'` | `'both'` |
| Processing time | ~2–4 hours | ~6–12 hours |
| CPU cores | 8 | 16 |

---

## Customization Guide

### Adjusting SLURM Resource Allocation

Edit the `#SBATCH` header lines in the `.sh` file:

```bash
#SBATCH --cpus-per-task=32   # More CPUs for very large datasets (>50M reads)
#SBATCH --time=24:00:00      # Extend time limit if needed
```

### Testing with a Subset of Reads

In the `.py` file, uncomment the `test_n_reads` line:

```python
pipeline = pipelines.TreblPipeline(
    ...
    test_n_reads=100000  # Test with first 100k reads
)
```

### Adjusting Reads Thresholds

Edit `READS_THRESHOLD`, `READS_THRESHOLD_AD`, and `READS_THRESHOLD_RT` in the `.py` file based on the reads distribution plots from Step 1 / Step 2. Start with `1` and increase if you want to filter low-count barcodes.

### Using a Different Conda Environment

Edit the `.sh` file:
```bash
source activate YOUR_ENV_NAME
```

---

## Troubleshooting

### Job Fails Immediately

```bash
cat logs/quick_start_*.err
```

Common causes:
- Conda environment not found → check environment name in the `.sh` file
- File not found → verify all paths in the `.py` CONFIGURATION section
- Permissions issue → `ls -l /path/to/your/files`

### Job Times Out

- Increase `--time` in the `.sh` SBATCH header
- Test with `test_n_reads=100000` first to estimate run time
- Use more CPUs: `--cpus-per-task=32`

### Out of Memory

- Increase `--cpus-per-task` (more CPUs = more memory allocated on Savio)
- Test with `test_n_reads` first

### Wrong Output Files / File Not Found

- Check `OUTPUT_DIR` and `DATA_DIR` in the `.py` file
- Confirm `mkdir -p output/... db logs` was run before submitting

### Import Errors

```bash
conda activate trebl_tools_env
pip show trebl_tools
```

---

## Job Monitoring Commands

```bash
squeue -u $USER          # Check queue position / status
scontrol show job JOBID  # Detailed job info
scancel JOBID            # Cancel a job
seff JOBID               # Resource usage after job completes
```

---

## Additional Resources

- **Main Examples README:** `../README.md`
- **Notebook Examples:** `../notebooks/`
- **Read the Docs:** User Guide → Running on Savio
- **Savio Documentation:** https://docs-research-it.berkeley.edu/services/high-performance-computing/
- **SLURM Documentation:** https://slurm.schedmd.com/sbatch.html
