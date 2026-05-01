# Advanced Usage

## Automatic Threshold Barcode Error Correction

Rather than manually choosing a threshold for error correction, UMI Tools supports automatic threshold detection based on the reads distribution. 

**Important:** This ONLY works if the data is very deeply sequenced. This means that the reads threshold you would pick manually would be over 100 and the average reads per sequence likely needs to be over 1000.

### To Use This Feature:

1. When initializing the pipeline, set `error_correction` to `True`.
2. In each step (step1, step2, trebl_experiment) set all reads thresholds to `1`.

## TREBL Experiment without UMIs

To analyze TREBL experiment data without UMIs, leave the `AD_umi_object` and `RT_umi_object` as `None` (the default). You must then specify `reads_threshold_AD` and `reads_threshold_RT` to filter low-read-count barcodes.

### Copy-Paste Ready Block

```python
# Pipeline setup — no error correction needed for no-UMI workflows
pipeline = pipelines.TreblPipeline(
    db_path="no_umi.db",
    design_file_path="design_file.txt",
    error_correction=False,
    output_path="output/no_umi"
)

# Define barcodes (same as Steps 1 and 2; no UMI objects needed)
AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="TGACTAG", length=120)
AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="GGGCCC", length=11)
RPTR_BC = finder.Barcode(name="RPTR_BC", preceder="CTCGAG", post="GGCCGC", length=14)

AD_bc_objects = [AD, AD_BC]
RT_bc_objects = [RPTR_BC]

trebl_AD_seq_files = glob.glob("/path/to/AD_Assembled/*")   # ← update with your path
trebl_RT_seq_files = glob.glob("/path/to/RPTR_Assembled/*")  # ← update with your path

# (Optional) Check reads distributions to pick threshold values
pipeline.trebl_experiment_reads_distribution(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True
)

# Run TREBL experiment — no UMI objects, use reads thresholds instead
pipeline.trebl_experiment_analysis(
    AD_seq_files=trebl_AD_seq_files,
    AD_bc_objects=AD_bc_objects,
    RT_seq_files=trebl_RT_seq_files,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    step1_map_csv_path="output/step1.csv",  # Path to Step 1 CSV
    AD_umi_object=None,                     # No UMI deduplication
    RT_umi_object=None,                     # No UMI deduplication
    reads_threshold_AD=10,                  # Set based on your reads distribution
    reads_threshold_RT=10,                  # Set based on your reads distribution
)
```

Use the reads distribution histograms to pick appropriate `reads_threshold_AD` and `reads_threshold_RT` values before running the analysis.

## Repeating an Analysis Step with Multiple Datasets

If you have multiple datasets, such as a spike-in vs a full sequencing run for step 1, you should initialize a separate pipeline, with a new db path and output path, for each new dataset. This ensures that all intermediate tables remain distinct.

For example, you could initialize a `spike_in_pipeline` where you use spike-in data for step1, step2, and the trebl experiment.

```python
## For spike in analysis
spike_in_pipeline = pipelines.TreblPipeline(
    # Path to save duck db file
    db_path = "spike_in.db",  
    # Path to design file
    design_file_path = "design_file.txt",
    # Whether to error correct sequences
    error_correction = False,   
    # Where to save output 
    output_path = "output/spike_in"
)
# Optional: test_n_reads = N, to try N reads first

# Run spike_in_pipeline.run_step_1() etc.
```

Then, once you have full sequencing data, create a new `full_pipeline` object with a new duckdb path before analyzing.

```python
## For full sequencing data analysis
full_pipeline = pipelines.TreblPipeline(
    # Path to save duck db file
    db_path = "full.db",  
    # Path to design file
    design_file_path = "design_file.txt",
    # Whether to error correct sequences
    error_correction = False,   
    # Where to save output 
    output_path = "output/full"
)
# Optional: test_n_reads = N, to try N reads first

# Run full_pipeline.run_step_1() etc.
```
