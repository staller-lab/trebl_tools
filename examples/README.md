# TREBL Tools Examples

This directory contains example notebooks and Savio job scripts for running TREBL analysis workflows. Two complete workflows are provided to accommodate different analysis needs.

## Directory Structure

```
examples/
├── notebooks/          # Jupyter notebooks for interactive analysis
├── savio_jobs/         # SLURM job scripts for cluster submission
├── data/              # Example data and documentation
└── README.md          # This file
```

## Workflows

### 1. Quick Start Workflow

**Best for:** Initial data exploration, testing, or when processing time is a priority

**Features:**
- No error correction (faster processing)
- Simple UMI deduplication only
- Typical runtime: 2-4 hours for standard datasets

**Files:**
- Notebook: `notebooks/quick_start_example.ipynb`
- Savio job: `savio_jobs/quick_start_job.sh`

### 2. Full Analysis Workflow

**Best for:** Final, publication-quality analysis when accuracy is the priority

**Features:**
- Error correction enabled (salvages more reads)
- Both simple and directional/complex UMI deduplication
- Typical runtime: 6-12 hours for standard datasets

**Files:**
- Notebook: `notebooks/full_analysis_example.ipynb`
- Savio job: `savio_jobs/full_analysis_job.sh`

## Key Differences

| Feature | Quick Start | Full Analysis |
|---------|-------------|---------------|
| Error Correction | ❌ No | ✅ Yes |
| UMI Deduplication | Simple only | Simple + Directional/Complex |
| Processing Time | Faster (~2-4 hrs) | Slower (~6-12 hrs) |
| Accuracy | Good | Best |
| Use Case | Initial exploration | Final analysis |

## Getting Started

### Option 1: Interactive Analysis (Jupyter Notebook)

1. Start a Jupyter session on Savio OOD
2. Navigate to `examples/notebooks/`
3. Open either `quick_start_example.ipynb` or `full_analysis_example.ipynb`
4. Update the file paths in the notebook to point to your data
5. Select the `trebl_env` kernel
6. Run the cells sequentially

**Note:** For large files (>10M reads), plotting steps can be time-consuming. Consider using the Savio job submission method instead.

### Option 2: Batch Job Submission (Recommended for Large Datasets)

1. Copy and customize the appropriate job script:
   ```bash
   cp examples/savio_jobs/quick_start_job.sh my_analysis_job.sh
   # or
   cp examples/savio_jobs/full_analysis_job.sh my_analysis_job.sh
   ```

2. Edit `my_analysis_job.sh` to update:
   - File paths (DESIGN_FILE, STEP1_SEQ_FILE, etc.)
   - Barcode sequences if different from examples
   - Resource allocation if needed (time, CPUs)

3. Create output and log directories:
   ```bash
   mkdir -p output/quick_start logs
   # or
   mkdir -p output/full_analysis logs
   ```

4. Submit the job:
   ```bash
   sbatch my_analysis_job.sh
   ```

5. Monitor progress:
   ```bash
   squeue -u $USER
   tail -f logs/quick_start_*.out
   # or
   tail -f logs/full_analysis_*.out
   ```

## Resource Requirements

### Quick Start
- **CPUs:** 8 cores
- **Memory:** ~16-32 GB (auto-allocated by Savio)
- **Time:** 4 hours (adjust if needed)
- **Partition:** savio3

### Full Analysis
- **CPUs:** 16 cores (more beneficial for error correction)
- **Memory:** ~32-64 GB (auto-allocated by Savio)
- **Time:** 12 hours (adjust if needed)
- **Partition:** savio3

## Customization

All example files use placeholder sequences and paths. Before running, you must:

1. **Update file paths** to point to your actual sequencing data
2. **Update barcode sequences** if your experiment uses different flanking sequences
3. **Adjust resource allocation** based on your dataset size
4. **Update the design file path** to your AD design sequences

## Output Files

Both workflows generate:

### Step 1 Outputs
- `step1_*.csv` - Mapping of barcodes from Step 1
- `step1_*_loss_table.png` - Visualization of filtering steps
- `step1_reads_distribution.png` - Histogram of reads per barcode

### TREBL Experiment Outputs
- `AD_trebl_experiment_results.csv` - AD library results with UMI counts
- `RT_trebl_experiment_results.csv` - RT library results with UMI counts
- `trebl_experiment_reads_distribution_*.png` - Histograms for all files
- `trebl_experiment_loss_*.png` - Loss summaries for AD and RT

### Database
- `*.db` - DuckDB database file (can be deleted after analysis)

## Troubleshooting

### Job Fails Due to Time Limit
- Increase the `--time` parameter in the SBATCH header
- For very large datasets, consider requesting up to 24 hours

### Out of Memory
- Error correction and complex UMI deduplication are memory-intensive
- Request more CPUs (memory is allocated per CPU on Savio)
- Consider using `test_n_reads` parameter for initial testing

### Module/Package Not Found
- Ensure you've activated the correct conda environment
- Check that all dependencies are installed: `pip list | grep trebl_tools`
- See main repository README for installation instructions

### Large File Warnings
- For files with >10M reads, expect longer processing times
- The plotting steps (`*_reads_distribution`) can take 10-30 minutes for large files
- Consider submitting as a job rather than running interactively

## Further Reading

- **Full Documentation:** See `docs/user_guide/` in the main repository
- **Installation Guide:** `docs/user_guide/installation.md`
- **Advanced Usage:** `docs/user_guide/advanced_usage.md`
- **API Reference:** See Read the Docs (link in main README)

## Questions or Issues?

If you encounter problems:
1. Check the log files in `logs/` directory
2. Review the main documentation in `docs/user_guide/`
3. Open an issue on the GitHub repository
