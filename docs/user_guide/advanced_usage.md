# Advanced Usage

## Automatic Threshold Barcode Error Correction

Rather than manually choosing a threshold for error correction, UMI Tools supports automatic threshold detection based on the reads distribution. 

**Important:** This ONLY works if the data is very deeply sequenced. This means that the reads threshold you would pick manually would be over 100 and the average reads per sequence likely needs to be over 1000.

### To Use This Feature:

1. When initializing the pipeline, set `error_correction` to `True`.
2. In each step (step1, step2, trebl_experiment) set all reads thresholds to `1`.

## TREBL Experiment without UMIs

To analyze TREBL experiment data without UMIs, leave the AD and UMI objects as `None`. However, must now specify `AD_reads_threshold` and `RT_reads_threshold` to be used across all files.

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
