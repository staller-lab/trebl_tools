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
RPTR_BC = finder.Barcode(name="RT_BC", preceder="CTCGAG", post="GGCCGC", length=14)

AD_bc_objects = [AD, AD_BC]
RT_bc_objects = [RPTR_BC]

AD_UMI = finder.Barcode(name="UMI", preceder="TGATTT", post="", length=12)
RT_UMI = finder.Barcode(name="UMI", preceder="TGTCAC", post="", length=12)

trebl_AD_seq_files = glob.glob("/path/to/AD_Assembled/*")   # ← update with your path
trebl_RT_seq_files = glob.glob("/path/to/RPTR_Assembled/*")  # ← update with your path

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
RPTR_BC = finder.Barcode(name="RT_BC", preceder="CTCGAG", post="GGCCGC", length=14)

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

### 3. TREBL Experiment Analysis

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

### 4. Tips

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

## Output

The TREBL experiment analysis saves several per-sample and experiment-level files automatically inside `output_path`.

### Per-sample outputs

For each sequencing file, a folder is created:

```text
output_path/
└── trebl_experiment_<suffix>/
    └── trebl_experiment_<suffix><sample_name>/
```

Inside each sample folder, the following files are saved:

#### UMI count tables

- `*_simple_umi_counts.tsv`
    - Counts each unique UMI separately.
    - Fast baseline deduplication.

- `*_directional_umi_counts.tsv`
    - Generated by UMI-tools directional deduplication.
    - Collapses related UMIs using sequence similarity and read counts.
    - Recommended for downstream analysis.

Both are generated automatically when:

```python
umi_deduplication="both"
```

#### Reads-per-UMI diagnostics

- `*_reads_per_umi.tsv`
    - Number of reads observed per UMI.
    - Useful for assessing PCR duplication and sequencing saturation.

#### QC figures

- Mapping/refinement loss plots
- Error-correction diagnostic plots (if `error_correction=True`)

These help visualize:
- total reads
- reads passing barcode quality filters
- reads retained after error correction
- reads retained after thresholding

---

### Experiment-level outputs

Saved directly in `output_path/`.

#### Activity score table

After running:

```python
pipeline.calculate_activity_scores(...)
```

the pipeline saves:

```text
trebl_experiment_activity_scores_per_barcode.csv
```

This file contains:
- AD and RT barcode identities
- replicate and timepoint information
- simple deduplicated counts
- directional deduplicated counts
- calculated activity scores

Key columns include:

```text
AD_count_simple
RT_count_simple
activity_simple
AD_count_directional
RT_count_directional
activity_directional
rep
time
```

Where "AD_count_simple" and "RT_count_simple" are the UMI counts from simple UMI deduplication, and "activity_simple" is the corresponding log10(RT/AD) activity score. "AD_count_directional" and "RT_count_directional" are the UMI counts from directional UMI deduplication, and "activity_directional" is the corresponding log10(RT/AD) activity score. 


#### TREBL experiment loss summary

If `step1_map_csv_path` is provided, the pipeline also saves experiment-wide QC figures in:

```text
output_path/figures/
```

including plots showing:
- total reads
- reads passing barcode quality checks
- reads matching the Step 1 map

across all sequencing files.

#### DuckDB intermediate tables

The pipeline also stores many intermediate refinement tables inside the DuckDB database file (`db_path`), including tables for:

```text
_initial
_quality
_error_corrected
_grouped
_thresholded
_designed
```

These are useful for debugging or inspecting intermediate filtering steps.

### Cleanup

Once you have finished your analysis and saved all results, you can delete the DuckDB database file to free up disk space:

```python
import os
os.remove("test.db")  # Replace with your db_path
```
