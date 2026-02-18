# Analysis Setup

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

**Important:** These variables remain fixed for downstream analysis. If you change any of these parameters, delete the DuckDB file and rerun the pipeline.

## Parameters

- **db_path**: Location to save the DuckDB database. Can be anywhere you have write permission. Okay to delete after you finish running analysis.
- **design_file_path**: Each line contains a DNA sequence of a designed AD. Examples: `/global/scratch/projects/fc_mvslab/OpenProjects/Marissa/DesignFiles/`.
- **error_correction**: Whether to apply error correction at each step. The error correction works by considering sequences to the right of a specified threshold in a reads distribution to be 'real.' Then, sequences to the left of the threshold which are similar enough will be corrected to a real sequence.
- **output_path**: Folder where CSVs and plots are saved.
- **test_n_reads**: Optional. Use a subset of reads for testing before full analysis.

## Tips

- Test with a small number of reads first using `test_n_reads`, and when switching to the full dataset, either use a separate `db_path` or delete the test `db_path` first—this prevents conflicts since the DuckDB database stores intermediate results and can't safely be reused with different parameters or more reads.
- Error correction will take longer, but it may salvage more reads. If you are time constrained, I would recommend skipping error correction.
- Make sure all paths (`db_path`, `design_file_path`, `output_path`) exist and are writable.
