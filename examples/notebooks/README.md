# Example Notebooks

This directory contains Jupyter notebooks demonstrating TREBL analysis workflows.

## Available Notebooks

### 1. `quick_start_example.ipynb`

A streamlined workflow for quick analysis and initial data exploration.

**Key Features:**
- No error correction (faster processing)
- Simple UMI deduplication only
- Ideal for testing or time-sensitive analysis

**When to Use:**
- Initial data exploration
- Testing new datasets
- When speed is more important than maximum accuracy
- When you have high-quality, deeply-sequenced data

**Expected Runtime:** ~2-4 hours for typical datasets

---

### 2. `full_analysis_example.ipynb`

A comprehensive workflow for publication-quality analysis.

**Key Features:**
- Error correction enabled (salvages reads with sequence errors)
- Both simple and directional/complex UMI deduplication
- Maximum accuracy for final results

**When to Use:**
- Final analysis for publication
- When data quality is variable
- When maximum accuracy is required
- When you need both simple and complex UMI deduplication for comparison

**Expected Runtime:** ~6-12 hours for typical datasets

---

## How to Use These Notebooks

### Prerequisites

1. **Environment Setup:**
   - Ensure you have installed `trebl_tools` and dependencies
   - See main repository README for installation instructions
   - Activate your conda environment: `conda activate trebl_env`

2. **Jupyter Setup on Savio:**
   - Start a Jupyter session via Savio Open OnDemand (OOD)
   - Install the kernel if not already done:
     ```bash
     python -m ipykernel install --user --name trebl_env
     ```
   - Select the `trebl_env` kernel when opening a notebook

### Running a Notebook

1. **Open the notebook** in Jupyter
2. **Update paths** in the configuration cells:
   - Design file path
   - Sequencing file paths
   - Output directory paths
3. **Customize barcodes** if your experiment uses different sequences
4. **Run cells sequentially** from top to bottom
5. **Review outputs** as they are generated

### Important Notes

#### For Large Files

If your sequencing files are large (>10 million reads), the plotting steps can take significant time:
- `step1_reads_distribution()` - Can take 10-30 minutes
- `trebl_experiment_reads_distribution()` - Can take 30-60 minutes for multiple files

**Recommendation:** For very large datasets, consider using the Savio job scripts in `../savio_jobs/` instead, which are optimized for cluster execution.

#### Comments About Savio Jobs

Throughout the notebooks, you'll see comments like:
```python
# NOTE: For large files (>10M reads), consider submitting this as a Savio job
```

These indicate steps that can be particularly time-consuming. The corresponding job scripts in `../savio_jobs/` handle these steps in a non-interactive environment, which is more efficient for large-scale processing.

### Comparing the Two Workflows

Both notebooks follow the same overall structure:
1. Setup and imports
2. Initialize pipeline
3. Define barcodes
4. Run Step 1 mapping
5. Run TREBL experiment analysis

**The key differences are:**

| Aspect | Quick Start | Full Analysis |
|--------|-------------|---------------|
| `error_correction` parameter | `False` | `True` |
| `umi_deduplication` parameter | `'simple'` | `'both'` |
| Loss table steps | Standard filtering | Includes 'error_corrected' step |
| Output columns | Simple UMI counts | Both simple and complex UMI counts |
| Processing time | Faster | Slower |

### Output Locations

By default, outputs are saved to:
- **Quick Start:** `output/quick_start/`
- **Full Analysis:** `output/full_analysis/`

These directories will contain:
- CSV files with mapping and count results
- PNG files with visualizations
- Loss table summaries

### Customization Tips

1. **Test with Subset:**
   Uncomment the `test_n_reads` parameter to process only the first N reads:
   ```python
   pipeline = pipelines.TreblPipeline(
       ...
       test_n_reads=100000  # Test with first 100k reads
   )
   ```

2. **Adjust Thresholds:**
   Use the reads distribution plots to pick appropriate `reads_threshold` values:
   ```python
   reads_threshold=10  # Adjust based on your distribution
   ```

3. **Modify Barcodes:**
   Update the `preceder`, `post`, and `length` parameters to match your experimental design:
   ```python
   AD = finder.Barcode(
       name="AD",
       preceder="YOUR_PRECEDER_SEQ",
       post="YOUR_POST_SEQ",
       length=YOUR_LENGTH
   )
   ```

## Troubleshooting

### Kernel Not Found
- Make sure you've installed the IPython kernel for your conda environment
- Run: `python -m ipykernel install --user --name trebl_env`
- Refresh your Jupyter session

### Import Errors
- Verify trebl_tools is installed: `pip show trebl_tools`
- Make sure you selected the correct kernel (`trebl_env`)
- Check that all dependencies are installed

### Long Running Cells
- For large files, some cells may take 10+ minutes
- Look for progress bars (tqdm) to monitor progress
- If too slow, consider using the Savio job scripts instead

### File Not Found Errors
- Double-check all path variables point to existing files
- Use absolute paths when possible
- Verify files have read permissions

## Next Steps

After running these notebooks:

1. **Review the outputs** in your output directory
2. **Examine loss tables** to understand data filtering
3. **Compare simple vs complex UMI counts** (full analysis only)
4. **Calculate activity scores** using the exported CSV files
5. **For production runs**, adapt the code into Savio job scripts for reproducibility

## Additional Resources

- **Main Examples README:** `../README.md`
- **User Guide:** `../../docs/user_guide/`
- **Savio Job Scripts:** `../savio_jobs/`
- **API Documentation:** See repository Read the Docs
