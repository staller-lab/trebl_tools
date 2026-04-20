# TREBL Experiment with UMIs

This step processes a full TREBL-seq experiment, mapping AD and RT barcodes across all sequencing files and generating counts. Here, we assume you will be using UMIs to account for PCR artifacts. See "TREBL Experiment without UMIs" in the [Advanced Usage](advanced_usage.md) page to do analysis without UMIs.

**Note:** Must have successfully run step 1 before.

## Choosing Your UMI Deduplication Mode

Pass `umi_deduplication` to `trebl_experiment_analysis` to control deduplication:

| Value | What it does | When to use |
|-------|-------------|-------------|
| `"simple"` | Count unique UMIs per barcode | Quick Start — faster, less disk usage |
| `"both"` | Run simple **and** directional (UMI-tools) deduplication | Full Analysis — publication-quality results |

> **Tip:** See the [Analysis Setup](analysis_setup.md) page for Quick Start vs Full Analysis presets, or jump straight to the example notebooks:
> [`quick_start_example.ipynb`](https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/quick_start_example.ipynb) |
> [`full_analysis_example.ipynb`](https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/full_analysis_example.ipynb)

## Copy-Paste Ready Block

### Quick Start (simple deduplication)

```python
# TREBL Experiment: Quick Start — simple UMI deduplication only

# Define barcodes (same as previous steps)
AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="TGACTAG", length=120)
AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="GGGCCC", length=11)
RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="GGCCGC", length=14)

AD_bc_objects = [AD, AD_BC]
RT_bc_objects = [RPTR_BC]

AD_UMI = finder.Barcode(name="UMI", preceder="TGATTT", post="", length=12)
RT_UMI = finder.Barcode(name="UMI", preceder="TGTCAC", post="", length=12)

trebl_AD_seq_files = glob.glob("/path/to/AD_Assembled/*")
trebl_RT_seq_files = glob.glob("/path/to/RPTR_Assembled/*")

pipeline.trebl_experiment_reads_distribution(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True
)

pipeline.trebl_experiment_analysis(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    step1_map_csv_path="output/step1.csv",   # Path to Step 1 CSV
    AD_umi_object=AD_UMI,
    RT_umi_object=RT_UMI,
    umi_deduplication="simple",              # Quick Start: simple only
)
```

### Full Analysis (both simple and directional deduplication)

```python
# TREBL Experiment: Full Analysis — simple + directional UMI deduplication

# (Same barcode and file definitions as Quick Start above)

pipeline.trebl_experiment_analysis(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    step1_map_csv_path="output/step1.csv",   # Path to Step 1 CSV
    AD_umi_object=AD_UMI,
    RT_umi_object=RT_UMI,
    umi_deduplication="both",               # Full Analysis: simple + directional
)
```

## Code Blocks Explained

### 1. Define Barcodes

```python
# Define AD and reporter barcodes (same as previous steps)
AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="TGACTAG", length=120)
AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="GGGCCC", length=11)
RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="GGCCGC", length=14)

# Group barcodes for convenience
AD_bc_objects = [AD, AD_BC]  # AD barcodes
RT_bc_objects = [RPTR_BC]    # Reporter barcodes

# Define UMIs for deduplication
AD_UMI = finder.Barcode(name="UMI", preceder="TGATTT", post="", length=12)
RT_UMI = finder.Barcode(name="UMI", preceder="TGTCAC", post="", length=12)
```

- AD, AD BC, and RT BC are defined the same way as Step 1 and 2.
- Now, we must also define UMI objects. Their name must be "UMI."
- If you want to take the last N bases, use `length=N`.
- Or, to use flanking sequences, refer to "Define barcodes" section of step 1 as a reminder for how extraction works.

### 2. Sequencing Files

```python
trebl_AD_seq_files = glob.glob("/path/to/AD_Assembled/*")
trebl_RT_seq_files = glob.glob("/path/to/RPTR_Assembled/*")
```

- Here we use `glob` to generate a list of fastq or fastq.gz files from a folder.
- Allows multiple sequencing runs to be analyzed together.

### 3. Reads Distribution

```python
pipeline.trebl_experiment_reads_distribution(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=[AD, AD_BC],
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=[RPTR_BC],
    reverse_complement=True
)
```

- Produces reads histograms for all AD and RT files.
- Quick QC: check read counts.

### 4. TREBL Experiment Analysis

```python
pipeline.trebl_experiment_analysis(
    AD_seq_files=trebl_AD_seq_files,            # AD sequencing files
    AD_bc_objects=AD_bc_objects,                # AD barcodes
    RT_seq_files=trebl_RT_seq_files,            # RT sequencing files
    RT_bc_objects=RT_bc_objects,                # RT barcodes
    reverse_complement=True,                    # Reverse complement?
    step1_map_csv_path="output/step1.csv",      # Path to Step 1 CSV file
    AD_umi_object=AD_UMI,                       # AD UMI deduplication
    RT_umi_object=RT_UMI,                       # RT UMI deduplication
    umi_deduplication="both",                   # "simple" or "both"
)
```

**`umi_deduplication` options:**
- `"simple"`: Count unique UMIs per barcode. Faster; used in the Quick Start notebook.
- `"both"`: Run simple deduplication **and** UMI-tools directional deduplication. Provides both
  counts for comparison; used in the Full Analysis notebook.

**Other outputs:**
- Refines mappings for AD and RT reads, filters low-read counts and keeps barcodes with correct lengths.
- Automatically saves CSVs and figures in `output_path`.
- Saves loss summary figure in output of total reads, reads after barcode quality check, and reads after comparing to step 1.

### 5. Tips

- Keep `step1_map_csv_path` pointing to the correct Step 1 CSV → ensures TREBL experiment maps correctly to Step 1.
- Figures help verify mappings, read distributions, and detect barcode collisions.
- Large experiments may take 30+ minutes, consider submitting as a Savio job.

## Final Steps

### Calculating Activity Scores

After completing the TREBL experiment analysis, calculate final activity scores by integrating AD and RT data through the Step 1 mapping. This produces gene-level activity metrics across timepoints and replicates.

```python
# Calculate activity scores
activity_scores = pipeline.calculate_activity_scores(
    step1_path="output/step1.csv",          # Path to Step 1 mapping CSV
    AD_bc_objects=AD_bc_objects,            # AD barcode objects
    RT_bc_objects=RT_bc_objects,            # RT barcode objects
    time_regex=r"_t(\d+)",                  # Regex to extract timepoint
    rep_regex=r"_r(\d+)"                    # Regex to extract replicate
)
```

#### Understanding the Parameters

- **`step1_path`**: Path to the Step 1 CSV file that maps AD barcodes to genes
- **`AD_bc_objects`** and **`RT_bc_objects`**: Same barcode objects used in previous steps
- **`time_regex`** and **`rep_regex`**: Regular expressions to extract metadata from sample names

#### Regex Convention for Time and Replicate Extraction

The function uses regular expressions to parse metadata from your sample filenames:

**Time Regex** (`r"_t(\d+)"`):
- Captures the numeric timepoint value from sample names
- Pattern matches: `_t` followed by one or more digits
- Example: `ChopTFs_AD_t10_r2.fastq` → extracts timepoint = `10`
- The parentheses `(\d+)` create a capture group for the numeric value

**Replicate Regex** (`r"_r(\d+)"`):
- Captures the numeric replicate identifier from sample names  
- Pattern matches: `_r` followed by one or more digits
- Example: `ChopTFs_AD_t10_r2.fastq` → extracts replicate = `2`
- The parentheses `(\d+)` create a capture group for the numeric value

**Important**: Your sample names must follow this naming convention:
- Include `_tXX` where XX is the timepoint (e.g., `_t0`, `_t10`, `_t60`)
- Include `_rX` where X is the replicate number (e.g., `_r1`, `_r2`, `_r3`)
- Example valid names: `sample_t24_r1.fastq`, `experiment_AD_t0_r2.fastq`, `data_t120_r3.fastq`

If your naming convention differs, adjust the regex patterns accordingly. For example:
- For `sample_time10_rep2.fastq` use: `time_regex=r"_time(\d+)"`, `rep_regex=r"_rep(\d+)"`
- For `sample_T10_R2.fastq` use: `time_regex=r"_T(\d+)"`, `rep_regex=r"_R(\d+)"`

#### Output

The function returns a multi-indexed DataFrame with:
- **Index**: (AD sequence, replicate) pairs
- **Columns**: Hierarchical structure with (timepoint, metric)

**Metrics calculated per timepoint**:
- `bc_activity_avg`: Mean per-barcode activity (RT_UMI / AD_UMI)
- `bc_activity_std`: Standard deviation of per-barcode activity
- `pooled_activity`: AD-level activity (sum of RT_UMI / sum of AD_UMI)

**Files saved** (if `output_path` configured):
- `bc_activities.csv`: Raw per-barcode activity scores
- `AD_activities.csv`: Consolidated activity scores table

The function integrates AD and RT experiment results to calculate activity ratios, providing both barcode-level (averaged) and gene-level (pooled) activity assessments.

### Cleanup

Once you have finished your analysis and saved all results, you can delete the DuckDB database file to free up disk space:

```python
import os
os.remove("test.db")  # Replace with your db_path
```
