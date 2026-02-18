# Savio Job Scripts

This directory contains SLURM batch scripts for running complete TREBL analysis workflows on the Savio cluster.

Each script covers all three analysis steps:
1. **Step 1**: Initial TREBL mapping to establish barcode relationships
2. **Step 2**: Process separated AD and RT libraries
3. **TREBL Experiment**: Full experiment analysis with UMI deduplication

## Available Scripts

### 1. `quick_start_job.sh`

Runs a quick TREBL analysis workflow optimized for speed.

**Configuration:**
- **CPUs:** 8 cores
- **Time limit:** 4 hours
- **Memory:** Auto-allocated (~16-32 GB)
- **Features:** No error correction, simple UMI deduplication only

**Best for:**
- Initial data exploration
- Testing pipelines
- Time-sensitive analysis

---

### 2. `full_analysis_job.sh`

Runs a comprehensive TREBL analysis with maximum accuracy.

**Configuration:**
- **CPUs:** 16 cores
- **Time limit:** 12 hours
- **Memory:** Auto-allocated (~32-64 GB)
- **Features:** Error correction enabled, both simple and directional/complex UMI deduplication

**Best for:**
- Final, publication-quality analysis
- Variable-quality data
- Comprehensive UMI deduplication comparison

---

## How to Use These Scripts

### Step 1: Prepare Your Environment

1. **Create necessary directories:**
   ```bash
   mkdir -p output/quick_start logs
   # or
   mkdir -p output/full_analysis logs
   ```

2. **Ensure conda environment is set up:**
   ```bash
   conda activate trebl_tools_env
   # Verify installation
   pip show trebl_tools
   ```

### Step 2: Customize the Script

1. **Copy the script to your working directory:**
   ```bash
   cp examples/savio_jobs/quick_start_job.sh my_analysis.sh
   # or
   cp examples/savio_jobs/full_analysis_job.sh my_analysis.sh
   ```

2. **Edit the configuration section** in the Python code within the script:

   ```python
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
   ```

3. **Optional: Customize barcode sequences** if your experiment differs from the examples:
   ```python
   AD = finder.Barcode(
       name="AD",
       preceder="YOUR_PRECEDER",
       post="YOUR_POST",
       length=YOUR_LENGTH
   )
   ```

4. **Optional: Adjust SLURM parameters** if needed:
   ```bash
   #SBATCH --time=4:00:00      # Increase if your data is very large
   #SBATCH --cpus-per-task=8   # More CPUs = more memory
   ```

### Step 3: Submit the Job

```bash
sbatch my_analysis.sh
```

### Step 4: Monitor the Job

**Check job status:**
```bash
squeue -u $USER
```

**View job output in real-time:**
```bash
tail -f logs/quick_start_*.out
# or
tail -f logs/full_analysis_*.out
```

**View any errors:**
```bash
tail -f logs/quick_start_*.err
# or
tail -f logs/full_analysis_*.err
```

### Step 5: Review Results

Once the job completes, check:
- Output CSV files in your output directory
- PNG visualizations of loss tables and distributions
- Log files for any warnings or errors

---

## Script Details

### SLURM Header Parameters

Both scripts use these common SLURM directives:

```bash
#SBATCH --account=fc_mvslab           # Billing account
#SBATCH --partition=savio3            # Partition (node type)
#SBATCH --nodes=1                     # Single node
#SBATCH --output=logs/NAME_%j.out     # Standard output log
#SBATCH --error=logs/NAME_%j.err      # Standard error log
```

**Difference between scripts:**

| Parameter | Quick Start | Full Analysis |
|-----------|-------------|---------------|
| `--cpus-per-task` | 8 | 16 |
| `--time` | 4:00:00 | 12:00:00 |

### Python Script Structure

Both scripts follow this workflow:

1. **Initialize pipeline** with appropriate settings
2. **Define barcodes** for AD, AD_BC, RT_BC, and UMIs
3. **Run Step 1:**
   - Plot reads distribution
   - Run mapping
4. **Run Step 2:**
   - Plot Step 2 reads distribution
   - Run Step 2 mapping
5. **Run TREBL experiment:**
   - Plot TREBL experiment reads distribution
   - Run analysis with UMI deduplication
6. **Output summary statistics**

### Key Differences Between Scripts

| Feature | Quick Start | Full Analysis |
|---------|-------------|---------------|
| `error_correction` | `False` | `True` |
| `umi_deduplication` | `'simple'` | `'both'` |
| Processing time | ~2-4 hours | ~6-12 hours |
| CPU cores | 8 | 16 |

---

## Customization Guide

### Adjusting Resource Allocation

#### For Larger Datasets

If you have very large datasets (>50M reads per file):

```bash
#SBATCH --cpus-per-task=32     # More CPUs for parallel processing
#SBATCH --time=24:00:00        # Extend time limit
```

#### For Smaller Datasets

If you have small datasets (<5M reads per file):

```bash
#SBATCH --cpus-per-task=4      # Fewer CPUs needed
#SBATCH --time=2:00:00         # Reduce time limit
```

### Testing with Subset of Data

Add `test_n_reads` parameter to process only first N reads:

```python
pipeline = pipelines.TreblPipeline(
    db_path="test.db",
    design_file_path=DESIGN_FILE,
    error_correction=False,
    output_path=OUTPUT_DIR,
    test_n_reads=100000  # Test with 100k reads
)
```

### Changing Output Locations

Update the `OUTPUT_DIR` variable and ensure the directory exists:

```python
OUTPUT_DIR = "/path/to/your/custom/output"
```

```bash
mkdir -p /path/to/your/custom/output
```

### Using Different Conda Environment

If your environment has a different name:

```bash
# Change this line in the script
source activate trebl_tools_env
# to
source activate YOUR_ENV_NAME
```

---

## Troubleshooting

### Job Fails Immediately

**Check the error log:**
```bash
cat logs/quick_start_*.err
# or
cat logs/full_analysis_*.err
```

**Common causes:**
- Conda environment not found → Check environment name
- Module not available → Ensure Python module is loaded
- Permissions issue → Check file/directory permissions

### Job Times Out

**Symptoms:** Job ends with TIMEOUT error

**Solutions:**
- Increase `--time` in SBATCH header
- Reduce dataset size for testing
- Use more CPUs to speed up processing

### Out of Memory

**Symptoms:** Job ends with OOM (Out of Memory) error

**Solutions:**
- Increase `--cpus-per-task` (more CPUs = more memory)
- Use `test_n_reads` to process subset first
- Split large files into smaller chunks

### File Not Found Errors

**Check:**
- All paths in CONFIGURATION section are correct
- Files have read permissions: `ls -l /path/to/file`
- Wildcards expand correctly: `ls /path/to/pattern/*`

### Import Errors

**Check:**
- Conda environment activated correctly
- trebl_tools installed: `conda activate trebl_tools_env && pip show trebl_tools`
- All dependencies installed

### Wrong Results Directory

**Issue:** Can't find output files

**Check:**
- `OUTPUT_DIR` variable in script
- Current working directory when job was submitted
- Use absolute paths for clarity

---

## Best Practices

### 1. Always Test First

Before running on full dataset:
```python
test_n_reads=100000  # Test with subset
```

### 2. Use Descriptive Job Names

Modify the job name to track different analyses:
```bash
#SBATCH --job-name=trebl_experiment1_quick
```

### 3. Organize Output

Create separate output directories for different experiments:
```bash
mkdir -p output/experiment1/quick_start
mkdir -p output/experiment1/full_analysis
```

### 4. Keep Logs

Don't delete log files until you've verified results:
```bash
logs/
├── quick_start_12345.out
├── quick_start_12345.err
├── full_analysis_67890.out
└── full_analysis_67890.err
```

### 5. Document Your Analysis

Add comments to your customized script:
```bash
# Analysis for Experiment 1 - ChopTF data from June 2025
# Using standard barcodes with 12bp UMIs
```

---

## Job Monitoring Tips

### Check Queue Position

```bash
squeue -u $USER
```

### Check Job Details

```bash
scontrol show job JOBID
```

### Cancel a Job

```bash
scancel JOBID
```

### View Resource Usage

After job completes:
```bash
seff JOBID
```

---

## Example Workflow

Here's a complete example of running a quick start analysis:

```bash
# 1. Set up directories
mkdir -p output/my_experiment/quick_start logs

# 2. Copy and customize script
cp examples/savio_jobs/quick_start_job.sh my_experiment_quick.sh
nano my_experiment_quick.sh  # Edit paths

# 3. Submit job
sbatch my_experiment_quick.sh

# 4. Monitor progress
tail -f logs/quick_start_*.out

# 5. Check results
ls output/my_experiment/quick_start/
```

---

## Additional Resources

- **Main Examples README:** `../README.md`
- **Notebook Examples:** `../notebooks/`
- **Savio Documentation:** https://docs-research-it.berkeley.edu/services/high-performance-computing/
- **SLURM Documentation:** https://slurm.schedmd.com/sbatch.html
