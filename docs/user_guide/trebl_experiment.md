# TREBL Experiment with UMIs

This step processes a full TREBL-seq experiment, mapping AD and RT barcodes across all sequencing files and generating counts. Here, we assume you will be using UMIs to account for PCR artifacts. See "TREBL Experiment without UMIs" in the [Advanced Usage](advanced_usage.md) page to do analysis without UMIs.

**Note:** Must have successfully run step 1 before.

Always use `umi_deduplication="both"` — this runs both simple and directional (UMI-tools) deduplication and gives you both sets of counts. The only choice is whether to enable **error correction** when you initialise the pipeline (see [Analysis Setup](analysis_setup.md)).

> **Tip:** Jump straight to the example notebooks:
> [`quick_start_example.ipynb`](https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/quick_start_example.ipynb) (no error correction) |
> [`full_analysis_example.ipynb`](https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/full_analysis_example.ipynb) (with error correction)

## Copy-Paste Ready Block

```python
# TREBL Experiment: Reads Distribution and Analysis
# Works for both Quick Start (error_correction=False) and
# Full Analysis (error_correction=True) — set that flag in TreblPipeline()

# Define barcodes (same as previous steps)
AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="TGACTAG", length=120)
AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="GGGCCC", length=11)
RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="GGCCGC", length=14)

AD_bc_objects = [AD, AD_BC]
RT_bc_objects = [RPTR_BC]

AD_UMI = finder.Barcode(name="UMI", preceder="TGATTT", post="", length=12)
RT_UMI = finder.Barcode(name="UMI", preceder="TGTCAC", post="", length=12)

trebl_AD_seq_files = glob.glob("/path/to/AD_Assembled/*")   # ← update with your path
trebl_RT_seq_files = glob.glob("/path/to/RPTR_Assembled/*")  # ← update with your path

# Check reads distributions to pick threshold values
pipeline.trebl_experiment_reads_distribution(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True
)

# Run full TREBL experiment — always use umi_deduplication="both"
pipeline.trebl_experiment_analysis(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    step1_map_csv_path="output/step1.csv",   # Path to Step 1 CSV
    AD_umi_object=AD_UMI,
    RT_umi_object=RT_UMI,
    umi_deduplication="both",               # Always use both
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
trebl_AD_seq_files = glob.glob("/path/to/AD_Assembled/*")   # ← update with your path
trebl_RT_seq_files = glob.glob("/path/to/RPTR_Assembled/*")  # ← update with your path
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
    umi_deduplication="both",                   # Always use "both"
)
```

**`umi_deduplication="both"`** runs two deduplication passes and saves both sets of counts:
- **Simple:** Count unique UMIs per barcode (fast baseline)
- **Directional:** UMI-tools directional deduplication (more accurate, corrects for PCR amplification bias)

Always use `"both"` — it takes a little longer but gives you the full picture.

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

After completing TREBL experiment analysis with `umi_deduplication="both"`, calculate activity scores from AD simple vs directional counts.

```python
# Calculate activity scores from AD simple + directional count tables
ad_activity_per_barcode_df = pipeline.calculate_activity_scores(
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

Activity score is computed as:

```python
np.log10(count_directional / count_simple)
```

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

The function returns one DataFrame:
- **`ad_activity_per_barcode_df`**: Per-barcode, per-sample merged AD simple/directional counts with activity scores

**Files saved** (if `output_path` configured):
- `AD_activity_scores_per_barcode.csv`

### Cleanup

Once you have finished your analysis and saved all results, you can delete the DuckDB database file to free up disk space:

```python
import os
os.remove("test.db")  # Replace with your db_path
```
