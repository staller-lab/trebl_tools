# Running on Savio (Cluster Jobs)

For large datasets, running TREBL analysis interactively in a Jupyter notebook is impractical. This page explains when to use a Savio batch job and how to run the provided scripts.

## When to Use a Savio Job

**Use a Savio batch job instead of a notebook when any of the following apply:**

- **Your sequencing files are large (>10M reads each)** — the reads distribution plot and mapping steps can each take 30+ minutes per file when run interactively.
- **You are running the full analysis with error correction** — error correction substantially increases processing time and should not be run in a session that can time out or be interrupted.
- **You have many TREBL experiment files** — processing 10+ FASTQ files sequentially benefits from dedicated cluster compute time.
- **You need a reproducible record** — job logs capture the complete run and any errors.

For small datasets or initial exploration, use the interactive notebooks:
[`quick_start_example.ipynb`](https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/quick_start_example.ipynb) |
[`full_analysis_example.ipynb`](https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/full_analysis_example.ipynb)

---

## Job Files

The job scripts are in `examples/savio_jobs/`. Each workflow is split into two files:

| File | What it is | What you do |
|------|-----------|-------------|
| `run_quick_start.py` | Python analysis script | **Edit this** — update paths and barcodes |
| `quick_start_job.sh` | SLURM wrapper | Submit this with `sbatch` |
| `run_full_analysis.py` | Python analysis script | **Edit this** — update paths and barcodes |
| `full_analysis_job.sh` | SLURM wrapper | Submit this with `sbatch` |

---

## Step-by-Step Guide

### 1. Choose a workflow

| | Quick Start | Full Analysis |
|--|-------------|---------------|
| **Error correction** | No | Yes |
| **UMI deduplication** | Simple + directional | Simple + directional |
| **Time limit** | 4 hours | 12 hours |
| **Best for** | Exploration / testing | Publication results |

### 2. Edit the Python analysis file

Open `run_quick_start.py` or `run_full_analysis.py` and update the **CONFIGURATION** section at the top of the file. This is the only section you need to edit.

**Paths:**

```python
# Set this to your data folder
DATA_DIR = "/global/scratch/projects/fc_mvslab/data/my_experiment"

DESIGN_FILE    = f"{DATA_DIR}/design_file.txt"
STEP1_SEQ_FILE = f"{DATA_DIR}/step1.fastq"
STEP2_AD_SEQ_FILE = f"{DATA_DIR}/step2_AD.fastq"
STEP2_RT_SEQ_FILE = f"{DATA_DIR}/step2_RT.fastq"

# List all TREBL experiment FASTQ files
# Names must include _tXX and _rX (e.g. sample_t10_r1.fastq)
AD_SEQ_FILES = [
    f"{DATA_DIR}/experiment_AD_t10_r1.fastq",
    f"{DATA_DIR}/experiment_AD_t60_r1.fastq",
]
RT_SEQ_FILES = [
    f"{DATA_DIR}/experiment_RT_t10_r1.fastq",
    f"{DATA_DIR}/experiment_RT_t60_r1.fastq",
]
```

**Barcodes** — update `preceder`, `post`, and `length` to match your library design:

```python
AD = finder.Barcode(
    name="AD",
    preceder="GGCTAGC",   # ← sequence immediately before the AD in your reads
    post="TGACTAG",       # ← sequence immediately after the AD
    length=120,
)

AD_BC = finder.Barcode(
    name="AD_BC",
    preceder="CGCGCC",
    post="GGGCCC",
    length=11,
)

RT_BC = finder.Barcode(
    name="RT_BC",
    preceder="CTCGAG",
    post="GGCCGC",
    length=14,
)

# UMI flanking sequences
AD_UMI = finder.Barcode(name="UMI", preceder="TGATTT", post="", length=12)
RT_UMI = finder.Barcode(name="UMI", preceder="TGTCAC", post="", length=12)
```

Refer to [Step 1: TREBL Mapping](step1.md) for a full explanation of barcode extraction parameters.

**Reads thresholds:**

```python
READS_THRESHOLD    = 1   # Step 1 — increase to filter low-read barcodes
READS_THRESHOLD_AD = 1   # Step 2 AD
READS_THRESHOLD_RT = 1   # Step 2 RT
```

### 3. Create required directories

```bash
mkdir -p output/quick_start db logs
# or for full analysis:
mkdir -p output/full_analysis db logs
```

### 4. Submit the job

```bash
sbatch examples/savio_jobs/quick_start_job.sh
# or
sbatch examples/savio_jobs/full_analysis_job.sh
```

### 5. Monitor the job

```bash
# Check queue status
squeue -u $USER

# Watch live output
tail -f logs/quick_start_*.out

# Check for errors
tail -f logs/quick_start_*.err
```

### 6. Review results

Once complete, outputs are in your `OUTPUT_DIR`:
- `step1.csv` — Step 1 barcode map (needed for Step 2 and TREBL experiment)
- `step2_AD.csv` / `step2_RT.csv` — Step 2 mappings
- `AD_activity_scores_per_barcode.csv` — per-barcode activity scores (`log10(directional/simple)`)
- PNG plots for each step

---

## Testing Before Running on Full Data

Add `test_n_reads` to process only the first N reads — useful to verify paths and barcodes before committing to a full run:

```python
pipeline = pipelines.TreblPipeline(
    ...
    test_n_reads=100000  # Process first 100k reads only
)
```

---

## Adjusting SLURM Resources

Edit the `#SBATCH` lines in the `.sh` file:

```bash
#SBATCH --cpus-per-task=32   # More CPUs for very large datasets
#SBATCH --time=24:00:00      # Extend time if needed
```

On Savio, memory is automatically allocated per CPU, so requesting more CPUs also gives more memory.

---

## Running with the Example Data

To test the scripts with the bundled example data, submit from the repo root without editing the `.py` file — `DATA_DIR` already points to `examples/data/`:

```bash
# From the repo root
mkdir -p output/quick_start db logs
sbatch examples/savio_jobs/quick_start_job.sh
```

---

## Troubleshooting

**Job fails immediately:**
```bash
cat logs/quick_start_*.err
```
- Check that all paths in the `.py` CONFIGURATION section exist
- Verify the conda environment name in the `.sh` file

**Job times out:** Increase `--time` in the `.sh` file, or use `test_n_reads` to estimate run time first.

**Out of memory:** Increase `--cpus-per-task` in the `.sh` file.

**Import errors:** Verify `trebl_tools` is installed in the active environment:
```bash
conda activate trebl_tools_env
pip show trebl_tools
```
