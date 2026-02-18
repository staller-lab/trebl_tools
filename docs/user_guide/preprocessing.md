# Preprocessing Data

Raw sequencing data should be filtered so that only high quality reads are kept.

## For Paired End Reads

For paired end reads, we recommend running [PEAR](https://cme.h-its.org/exelixis/web/software/pear/doc.html). Most of the data will be paired end reads, except for UMI reporter reads.

## For Single End Reads

For single end reads, we recommend running [fastp](https://github.com/OpenGene/fastp).

### Parameters

- **input_dir**: folder with `.fastq` and/or `.fastq.gz` files to process
- **output_dir**: where to save filtered files

```python
preprocess.run_fastp(input_dir, output_dir)
```

### Get Summary

For a summary of how many reads fastp removed:

```python
preprocess.fastp_summary_df(output_dir)
```
