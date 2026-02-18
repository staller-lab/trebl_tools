# TREBL Experiment with UMIs

This step processes a full TREBL-seq experiment, mapping AD and RT barcodes across all sequencing files and generating counts. Here, we assume you will be using UMIs to account for PCR artifacts. See "TREBL Experiment without UMIs" in the [Advanced Usage](advanced_usage.md) page to do analysis without UMIs.

**Note:** Must have successfully run step 1 before.

## Copy-Paste Ready Block

```python
# TREBL Experiment: Reads Distribution and Analysis

# Define barcodes (same as previous steps)
AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="TGACTAG", length=120)
AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="GGGCCC", length=11)
RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="GGCCGC", length=14)

# Separate objects by AD and RT
AD_bc_objects = [AD, AD_BC]  # AD and AD barcodes, only use AD BC if human
RT_bc_objects = [RT_BC]      # Reporter barcodes

# Separate AD and RT objects to search for
AD_UMI = finder.Barcode(name="UMI", preceder="TGATTT", post="", length=12)
RT_UMI = finder.Barcode(name="UMI", preceder="TGTCAC", post="", length=12)

# Use glob to collect all .fastq files for AD and RT
trebl_AD_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/data/sequencing/2025-10-02_ChopTF_TREBL_spike/results/AD_Assembled/*")
trebl_RT_seq_files = glob.glob("/global/scratch/projects/fc_mvslab/data/sequencing/2025-10-02_ChopTF_TREBL_spike/results/RPTR_Assembled/*")

# Run TREBL experiment reads distribution
pipeline.trebl_experiment_reads_distribution(
    AD_seq_files=trebl_AD_seq_files,       # List of AD sequencing files
    AD_bc_objects=[AD, AD_BC],             # AD barcodes
    RT_seq_files=trebl_RT_seq_files,       # List of RT sequencing files
    RT_bc_objects=[RPTR_BC],               # RT barcodes
    reverse_complement=True                 # Reverse complement reads before searching
)
# Generates histograms for all AD and RT files
# Helps visualize distribution across the full experiment

# Run full TREBL experiment analysis
trebl_results = pipeline.trebl_experiment_analysis(
    AD_seq_files=trebl_AD_seq_files,       # AD sequencing files
    AD_bc_objects=[AD, AD_BC],             # AD barcodes
    RT_seq_files=trebl_RT_seq_files,       # RT sequencing files
    RT_bc_objects=[RPTR_BC],               # RT barcodes
    reverse_complement=True,               # Reverse complement reads
    step1_map_name="step1_AD_AD_BC_RPTR_BC_designed", # Step 1 mapping to link reads
    AD_umi_object=AD_UMI,                  # Optional: AD UMI deduplication
    RT_umi_object=RT_UMI,                  # Optional: RT UMI deduplication
    step_name_prefix="trebl_experiment_"   # Prefix for output files
)
# Produces final TREBL experiment output:
# - CSVs of AD and RT counts
# - Visualizations of distributions and mappings
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

# New: Define UMIs for deduplication
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
trebl_results = pipeline.trebl_experiment_analysis(
    AD_seq_files=trebl_AD_seq_files,            # AD sequencing files
    AD_bc_objects=AD_bc_objects,                # AD barcodes
    RT_seq_files=trebl_RT_seq_files,            # RT sequencing files
    RT_bc_objects=RT_bc_objects,                # RT barcodes
    reverse_complement=True,                    # Reverse complement?
    step1_map_name="step1_AD_AD_BC_RPTR_BC_designed",  # Step 1 mapping
    AD_umi_object=AD_UMI,                       # AD UMI deduplication
    RT_umi_object=RT_UMI,                       # RT UMI deduplication
    step_name_prefix="trebl_experiment_"        # Prefix for output files
)
```

- Performs both simple and complex UMI deduplication.
  - **Simple:** Count unique UMIs
  - **Complex:** Use umi tools directional deduplication
- Refines mappings for AD and RT reads, filters low-read counts and keeps barcodes with correct lengths.
- Automatically saves CSVs and figures in `output_path`.
- Saves loss summary figure in output of total reads, reads after barcode quality check, and reads after comparing to step 1.
- Saves trebl_results DataFrame with aggregated counts to output path

### 5. Tips

- Keep `step1_map_name` correct → ensures TREBL experiment maps correctly to Step 1.
- Figures help verify mappings, read distributions, and detect barcode collisions.
- Large experiments may take 30+ minutes, consider submitting as a Savio job.

## Final Steps

[To-do] Get final activities

Finally, you can delete the duck db file.
