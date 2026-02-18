# Step 2: TREBL Step 2 Mapping

The example below is again Marissa's Step 2 data. In step 2, the barcodes are now separated so are sequenced in different files.

**Note:** Must have first run step 1 successfully.

## Copy-Paste Ready Block

```python
# Step 2: TREBL Step 2 Mapping (Full Block)

# Define barcodes (same as Step 1)
AD = finder.Barcode(name="AD", preceder="GGCTAGC", post="TGACTAG", length=120)
AD_BC = finder.Barcode(name="AD_BC", preceder="CGCGCC", post="GGGCCC", length=11)
RT_BC = finder.Barcode(name="RT_BC", preceder="CTCGAG", post="GGCCGC", length=14)

# Separate objects by AD and RT
AD_bc_objects = [AD, AD_BC]  # AD and AD barcodes, only use AD BC if human
RT_bc_objects = [RT_BC]      # Reporter barcodes

# Sequencing files for AD and RT
step2_AD_seq_file = "/global/scratch/projects/fc_mvslab/data/sequencing/czb_new_sept2025/MAZ10/ChopTF/results/AD_Assembled/ChopTFstep2_2_AD_concat.fastq"
step2_RT_seq_file = "/global/scratch/projects/fc_mvslab/data/sequencing/czb_new_sept2025/MAZ10/ChopTF/results/RPTR_Assembled/ChopTFstep2_2_RPTR_concat.fastq"

# Plot reads distribution before thresholding
pipeline.step2_reads_distribution(
    AD_seq_file=step2_AD_seq_file,       # AD sequencing file
    AD_bc_objects=AD_bc_objects,         # AD barcodes to search
    RT_seq_file=step2_RT_seq_file,       # RT sequencing file
    RT_bc_objects=RT_bc_objects,         # RT barcodes to search
    reverse_complement=True              # Search reverse complement
)
# Produces histograms for AD and RT reads
# Helps pick reads_threshold_AD and reads_threshold_RT

# Run Step 2 mapping
step2 = pipeline.run_step_2(
    AD_seq_file=step2_AD_seq_file,       # AD sequencing file
    AD_bc_objects=AD_bc_objects,         # AD barcodes
    RT_seq_file=step2_RT_seq_file,       # RT sequencing file
    RT_bc_objects=RT_bc_objects,         # RT barcodes
    reverse_complement=True,             # Search reverse complement
    reads_threshold_AD=10,               # Minimum reads for AD
    reads_threshold_RT=10,               # Minimum reads for RT
    step1_map_name="step1_AD_AD_BC_RPTR_BC_designed"  # Previous step 1 map in DuckDB
)

# Extract outputs
AD_step2 = step2["AD_step2"]                # Step 2 AD map
RT_step2 = step2["RT_step2"]                # Step 2 RT map
step1_step2_overlap = step2["step1_overlap"] # Overlap with Step 1 map
```

## Code Blocks Explained

### 1. Define Barcodes You Will Search For in AD and RT File

```python
# Same barcodes as step 1 
AD = finder.Barcode(
    name="AD", 
    preceder="GGCTAGC", 
    post="TGACTAG", 
    length=120
)
AD_BC = finder.Barcode(
    name="AD_BC", 
    preceder="CGCGCC", 
    post="GGGCCC", 
    length=11
)
RT_BC = finder.Barcode(
    name="RT_BC", 
    preceder="CTCGAG", 
    post="GGCCGC", 
    length=14
)
```

Refer to "Define barcodes" section of step 1 page for explanation of how barcode extraction works.

For human analysis, unlike the yeast example above, you will only have AD BC and RT BC, in your step 2 data. This is because the AD and AD BC will no longer be adjacent.

**Tips:**
- Almost always the same as step 1. Should only be different if you need to use shorter/longer flanking sequences.

### 2. Split Objects by AD and RT

```python
AD_bc_objects = [AD, AD_BC]
RT_bc_objects = [RT_BC]
```

In step 2, AD and RT barcodes are in separate files so must be analyzed separately.

### 3. Specify Sequencing Files

```python
step2_AD_seq_file = "path/to/AD_reads.fastq"
step2_RT_seq_file = "path/to/RT_reads.fastq"
```

- Separate files for AD and RT sequencing reads.
- Can be a single file or list of files; `.fastq` or `.fastq.gz`.

### 4. Plot Reads Distribution

```python
pipeline.step2_reads_distribution(
    AD_seq_file=step2_AD_seq_file,
    AD_bc_objects=AD_bc_objects,
    RT_seq_file=step2_RT_seq_file,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True
)
```

- Generates AD and RT reads histograms. Will appear inline in jupyter notebook. The figure will be saved to the output folder.
- Helps pick `reads_threshold_AD` and `reads_threshold_RT` to get high quality data.
- Visual check for sequencing quality and read coverage.

### 5. Run Step 2 Mapping

```python
step2 = pipeline.run_step_2(
    AD_seq_file=step2_AD_seq_file,
    AD_bc_objects=AD_bc_objects,
    RT_seq_file=step2_RT_seq_file,
    RT_bc_objects=RT_bc_objects,
    reverse_complement=True,
    reads_threshold_AD=10,
    reads_threshold_RT=10,
    step1_map_name="step1_AD_AD_BC_RPTR_BC_designed"
)
```

- Applies reads thresholds to ensure low-read noise is removed.
- Checks how Step 2 reads overlap with Step 1 mapping.

### 6. Outputs

```python
AD_step2 = step2["AD_step2"]
RT_step2 = step2["RT_step2"]
step1_step2_overlap = step2["step1_overlap"]
```

- **AD_step2** → DataFrame of extracted sequences from AD file after Step 2 processing. One row per unique combination with reads count.
- **RT_step2** → DataFrame of RT sequences after Step 2 processing. One row per unique combination with reads count.
- **step1_step2_overlap** → DataFrame showing overlap between Step 1 map and Step 2 reads.
- CSVs and PNGs are saved automatically.
