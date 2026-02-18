# Step 1: TREBL Mapping

The example below is Marissa's Chop TFs data on Savio. Feel free to follow along, this step 1 code will likely need around 10-20 minutes to run.

## Copy-Paste Ready Block

```python
# Step 1: TREBL Mapping (Full Block)

# Define barcodes to search for in reads
AD = finder.Barcode(
    name="AD",             # Name identifier for activation domain
    preceder="GGCTAGC",    # Sequence just before the AD in reads
    post="TGACTAG",        # Sequence just after the AD in reads
    length=120             # Expected length of the AD
)

AD_BC = finder.Barcode(
    name="AD_BC",          # Name for activation domain barcode
    preceder="CGCGCC",     # Sequence before barcode
    post="GGGCCC",         # Sequence after barcode
    length=11              # Barcode length
)

RT_BC = finder.Barcode(
    name="RT_BC",          # Reporter barcode
    preceder="CTCGAG",     # Sequence before reporter barcode
    post="GGCCGC",         # Sequence after reporter barcode
    length=14              # Reporter barcode length
)

# Combine all barcodes into a single list for easier input
bc_objects = [AD, AD_BC, RT_BC]

# Specify sequencing file(s)
step1_seq_file = "/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2024/A10_A11/results/A11_S2.fastq.gz.assembled.fastq"
# Can be a single file (string) or multiple files (list of strings)
# Supported formats: .fastq or .fastq.gz

# Plot reads distribution
pipeline.step1_reads_distribution(
    seq_file=step1_seq_file,   # FASTQ file to analyze
    bc_objects=bc_objects,     # List of barcodes to search for
    reverse_complement=True    # Search both forward and reverse complement
)
# Produces histogram of reads per barcode
# Helps pick appropriate reads_threshold for filtering
# In Jupyter: plot shows inline; in scripts/jobs: plot is saved as PNG

# Run Step 1 mapping
step1_map = pipeline.run_step_1(
    seq_file=step1_seq_file,         # FASTQ file input
    bc_objects=bc_objects,           # Barcodes to map
    column_pairs=[("RT_BC", "AD")],  # Check for collisions between reporter barcode and AD
    reads_threshold=10,              # Minimum number of reads to keep a barcode
    reverse_complement=False         # Do not reverse complement reads for mapping
)
# Returns a DataFrame of Step 1 mapping
# Saves:
# - CSV of step1_map (needed for Step 2)
# - PNG of loss_table visualization
# - Optional CSV of loss_table
```

## Code Blocks Explained

### 1. Define Barcodes

```python
# Sequences to search for in reads
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

bc_objects = [AD, AD_BC, RT_BC]  # Combine barcodes for input
```

Each barcode object defines what TREBL should extract from reads. `preceder` and `post` sequences help locate the barcode/AD in the read.

**Extraction rules:**
- If both `preceder` and `post` are provided, I extract the sequence between them.
- If only one of `preceder` or `post` is provided, I extract a sequence starting from the provided flanking sequence and extending `length` bases.
- If both `preceder` and `post` are empty, TREBL extracts the last `length` bases of the read.

**Tips:**
- Recommend providing both `preceder` and `post` for most stringent search. Can use only `pre` or `post` for a more flexible search, will use more reads.
- Always provide the `length`.

### 2. Specify Sequencing Files

```python
step1_seq_file = "/global/scratch/projects/fc_mvslab/data/sequencing/CZB_Feb2024/A10_A11/results/A11_S2.fastq.gz.assembled.fastq"
```

- Can be a single file or multiple files in a list.
- Can be in `.fastq` or `.fastq.gz` formats.

### 3. Plot Reads Distribution

```python
pipeline.step1_reads_distribution(
    seq_file=step1_seq_file,
    bc_objects=bc_objects,
    reverse_complement=True
)
```

- Generates histogram of reads per barcode. Will appear inline in jupyter notebook. The figure will be saved to the output folder.
- Can use this distribution to decide `reads_threshold` for filtering.
- **Tip:** Unusual distributions may indicate sequencing issues. If there is no clear exponential drop off with a peak after, the data may be undersequenced.

### 4. Run Step 1 Mapping

```python
step1_map = pipeline.run_step_1(
    seq_file=step1_seq_file,
    bc_objects=bc_objects,
    column_pairs=[("RT_BC", "AD")],
    reads_threshold=10,
    reverse_complement=False
)
```

**Explanation of `column_pairs`:**

`column_pairs` specifies which barcodes to check for collisions in the reads. We check that >90% of reads from the first column correspond to a single entry in the second column. This ensures that each barcode can be reliably mapped to its corresponding activation domain (AD).

**Column pair advanced usage:**

Advanced usage: Can check multiple pairs of columns, or even concatenate multiple columns. Check documentation for `run_step_1` function for more information.

**Examples:**

| Example Type | AD → RT Mapping | Comment |
|-------------|-----------------|---------|
| **Keep** - Clear 1:1 | AD 1 → RT 1 (100 reads) | RT maps clearly to a single AD |
| **Remove!** - Ambiguous | AD 1 → RT 1 (100 reads)<br>AD 2 → RT 1 (100 reads) | RT 1 cannot be assigned uniquely → collision |
| **Keep** - Mostly Clear | AD 1 → RT 1 (100 reads)<br>AD 2 → RT 1 (10 reads) | >90% of reads from one AD → acceptable mapping |
| **Keep** - Multiple RTs | AD 1 → RT 1 (100 reads)<br>AD 1 → RT 2 (100 reads) | One AD mapping to multiple RTs → okay (replicates) |

**Outputs:**

- Produces Step 1 map DataFrame.
- CSV and loss table are saved automatically.
- `column_pairs` ensures that unexpected barcode collisions are detected.

### 5. Interpreting Outputs

Loss table rows represent steps of filtering:
- **grouped** → unique reads counted
- **thresholded** → reads above threshold
- **barcode_exists** → only reads with non-empty barcodes
- **unique_target** → checks for barcode collisions
- **quality** → correct length reads
- **designed** → only ADs present in design file

The darker colors/number on bottom represent unique rows while the number to the right and above represents the total number of reads represented by those unique reads.

- Check whether too many reads are being filtered out.
- Keep CSV for Step 2 mapping.
