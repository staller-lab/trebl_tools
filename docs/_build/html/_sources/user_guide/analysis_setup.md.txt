# Analysis Setup

## Recommended Presets

Choose one of the two presets below based on your goal. You can follow along with the
matching example notebook on GitHub:
[`quick_start_example.ipynb`](https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/quick_start_example.ipynb) |
[`full_analysis_example.ipynb`](https://github.com/staller-lab/trebl_tools/blob/main/examples/notebooks/full_analysis_example.ipynb)

**Both presets always use `umi_deduplication="both"`** (simple + directional deduplication). The only thing you choose is whether to enable error correction.

### What are the UMI deduplication options?
**UMI deduplication** determines how UMIs are collapsed into counts.  

- **Simple deduplication** counts every unique UMI independently.  
- **Directional deduplication** follows the approach implemented in [UMI-tools](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html), using both UMI sequence similarity and read counts to collapse likely sequencing errors into a single UMI.  

Directional deduplication is generally more accurate and is recommended for most analyses. It adds only a modest amount of runtime.

### What are the error correction options?

**Error correction** attempts to rescue low-count reads that are likely caused by sequencing errors. The method assumes that high-count sequences (above a selected read-count threshold) represent real sequences, while low-count sequences below the threshold may be erroneous variants. Similar low-count sequences are then corrected to their higher-count counterparts.

Error correction is especially useful for libraries containing many highly similar sequences and can increase the number of usable reads. However, it works best with deeply sequenced datasets. In practice, a clear threshold often requires >100 reads per sequence, which may not always be achievable. You can either specify the read-count threshold manually or allow the tool to determine it automatically, as described in `advanced_usage.md`.

### Quick Start (no error correction, faster, recommended)

This is the standard and recommended pipeline. 

```python
pipeline = pipelines.TreblPipeline(
    db_path="quick_start.db",
    design_file_path="design_file.txt",
    error_correction=False,          # Skip error correction — faster
    output_path="output/quick_start"
)
# Optional: test_n_reads=100000 to test with a subset first
```

Later, when calling `trebl_experiment_analysis`, always use:

```python
umi_deduplication="both"             # Simple + directional (complex) deduplication
```

### Full Analysis (with error correction)

Use for error correction to attempt to salvage more reads. Usually not preferred. 

```python
pipeline = pipelines.TreblPipeline(
    db_path="full_analysis.db",
    design_file_path="design_file.txt",
    error_correction=True,           # Enable UMI-tools error correction
    output_path="output/full_analysis"
)
# Optional: test_n_reads=100000 to test with a subset first
```

Later, when calling `trebl_experiment_analysis`, always use:

```python
umi_deduplication="both"             # Simple + directional (complex) deduplication
```

**Important:** These variables remain fixed for downstream analysis. If you change any of these parameters, delete the DuckDB file and rerun the pipeline.

## Initialize the Pipeline

```python
pipeline = pipelines.TreblPipeline(
    # Path to save duck db file
    db_path = "test.db",  
    # Path to design file
    design_file_path = "design_file.txt",
    # Whether to error correct sequences
    error_correction = False,   
    # Where to save output 
    output_path = "../../output/NKX2-2/pipeline_test/no_err_corr"
)
# Optional: test_n_reads = N, to try N reads first
```

## Parameters

- **db_path**: Location to save the DuckDB database. Can be anywhere you have write permission. Okay to delete after you finish running analysis.
- **design_file_path**: Each line contains a DNA sequence of a designed AD. Examples: `/global/scratch/projects/fc_mvslab/OpenProjects/Marissa/DesignFiles/`.
- **error_correction**: Whether to apply error correction at each step. The error correction works by considering sequences to the right of a specified threshold in a reads distribution to be 'real.' Then, sequences to the left of the threshold which are similar enough will be corrected to a real sequence.
- **output_path**: Folder where CSVs and plots are saved.
- **test_n_reads**: Optional. Use a subset of reads for testing before full analysis.

## Tips

- Test with a small number of reads first using `test_n_reads`, and when switching to the full dataset, either use a separate `db_path` or delete the test `db_path` first—this prevents conflicts since the DuckDB database stores intermediate results and can't safely be reused with different parameters or more reads.
- Error correction will take longer, but it may salvage more reads. I would usually recommend skipping error correction.
- Make sure all paths (`db_path`, `design_file_path`, `output_path`) exist and are writable.
