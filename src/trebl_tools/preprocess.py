import dask.dataframe as dd
import pandas as pd
import subprocess
from Bio.Seq import Seq
from dask.diagnostics import ProgressBar
import os
from dask import delayed, compute
import pathlib
from itertools import islice
from tqdm import tqdm
import time
from functools import wraps
import gzip
from pathlib import Path
import re

def time_it(func):
    """Decorator to print how long a method takes to run.
    
    Args:
        func (callable): The function to be decorated.
        
    Returns:
        callable: The wrapped function that prints execution time.
        
    Note:
        Time is displayed in seconds if less than 60 seconds, otherwise in minutes.
    """    
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        start = time.time()
        result = func(self, *args, **kwargs)
        elapsed = time.time() - start
        if elapsed < 60:
            print(f"Done in {elapsed:.2f} seconds.\n")
        else:
            print(f"Done in {elapsed/60:.2f} minutes.\n")
        return result
    return wrapper


def load_and_shorten_files(seq_files, reverse_complement):
    """Load sequence files (FASTQ, gzipped FASTQ, or TXT) as Dask DataFrame.
    
    For FASTQ files, keeps only every 4th line starting from line 1 (the sequence line).
    
    Args:
        seq_files (list of str or str): List of file paths or single file path.
        reverse_complement (bool): Whether to apply reverse complement to sequences.
    
    Returns:
        dask.DataFrame: Loaded sequences as a Dask DataFrame with column 'sequence'.
        
    Note:
        Supported file extensions: .fastq, .fq, .gz for FASTQ files.
        Other extensions are treated as plain text files with one sequence per line.
    """
    # Normalize input
    if isinstance(seq_files, str):
        seq_files = [seq_files]

    dfs = []
    for f in seq_files:
        ext = pathlib.Path(f).suffix.lower()
        
        if ext in [".fastq", ".fq", ".gz"]:
            # Read all lines
            df = dd.read_csv(f, header=None, names=["raw"], dtype=str)

            # Keep only rows 1, 5, 9, ... → sequence lines
            df = df.map_partitions(lambda d: d.iloc[1::4])
            df = df.rename(columns={"raw": "sequence"})
        else:
            # Assume plain TXT file with one sequence per line
            df = dd.read_csv(f, header=None, names=["sequence"], dtype=str)

        dfs.append(df)

    # Concatenate all inputs
    all_seq_df = dd.concat(dfs)

    # Optionally apply reverse complement
    if reverse_complement:
        all_seq_df["sequence"] = all_seq_df["sequence"].map_partitions(
            reverse_complement_series,
            meta=("sequence", str)
        )

    return all_seq_df

def reverse_complement_series(s: pd.Series) -> pd.Series:
    """Return the reverse complement of a pandas Series of DNA sequences.

    Args:
        s (pd.Series): Series of DNA sequences.

    Returns:
        pd.Series: Series of reverse complemented sequences.

    Example:
        >>> s = pd.Series(["ATGC", "GATTACA"])
        >>> reverse_complement_series(s)
        0       GCAT
        1    TGTAATC
        dtype: object
        
    Note:
        If an error occurs during reverse complement conversion, 
        an empty string is returned for that sequence.
    """
    def safe_rc(seq):
        try:
            return str(Seq(str(seq)).reverse_complement())
        except Exception:
            print(seq)
            return ""
    return s.map(safe_rc)

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence.
    
    Args:
        seq (str): DNA sequence string.
        
    Returns:
        str: Reverse complement of the input sequence.
        
    Example:
        >>> reverse_complement("ATGC")
        'GCAT'
    """
    return str(Seq(str(seq)).reverse_complement())

    
def shorten_seq_file(infile, outfile, chunk_size=4*100000):
    """Extract sequences from FASTQ files to a text file (fastest method).
    
    Extracts the 2nd line of each FASTQ record (sequence line) and writes
    one sequence per line to the output file. Works with both plain and 
    gzipped FASTQ files.
    
    Args:
        infile (str): Input FASTQ or FASTQ.GZ file path.
        outfile (str): Output text file path (1 sequence per line).
        chunk_size (int, optional): Number of lines to read at once. 
            Must be multiple of 4. Defaults to 400000.
            
    Note:
        Progress is displayed using tqdm progress bar.
        
    Example:
        >>> shorten_seq_file("input.fastq.gz", "sequences.txt")
        Shortening input.fastq.gz to sequences.txt...
    """
    print(f"Shortening {infile} to {outfile}...\n")
    opener = gzip.open if infile.endswith(".gz") else open

    with opener(infile, "rt") as fin, open(outfile, "w") as fout:
        pbar = tqdm(desc="Processing chunks", unit="chunk")
        while True:
            lines = list(islice(fin, chunk_size))
            if not lines:
                break
            fout.writelines(lines[1::4])  # Write every 2nd line of 4-line record
            pbar.update(1)
        pbar.close()

                
def save_parquet(df, path):
    """Save a DataFrame to a Parquet file with progress bar.

    Args:
        df (dask.DataFrame): DataFrame to save.
        path (str): Path where the Parquet file will be saved.

    Returns:
        None

    Warning:
        If the path already exists, data will be appended to it.

    Example:
        >>> save_parquet(my_dataframe, "output.parquet")
        Warning: The path 'output.parquet' already exists and will be added to.
    """
    if os.path.exists(path):
        print(f"Warning: The path '{path}' already exists and will be added to.\n")

    with ProgressBar():
        df.to_parquet(
            path,
            engine='pyarrow',
            write_index=False
        )

def run_fastp(input_dir, output_dir, script_path="savio_jobs/fastp.sh"):
    """Run fastp on all .fastq.gz files using SLURM array job.
    
    Submits a SLURM array job to process all .fastq.gz files in the input
    directory using an existing fastp.sh script. Outputs processed files
    as .fastq files.
    
    Args:
        input_dir (str): Path to folder containing .fastq or .fastq.gz files.
        output_dir (str): Path to folder where output should be written.
        script_path (str, optional): Path to the existing fastp.sh script.
            Defaults to "savio_jobs/fastp.sh".
    
    Raises:
        FileNotFoundError: If the script_path does not exist.
        ValueError: If no .fastq.gz files are found in input_dir.
        
    Note:
        Creates a logs subdirectory in output_dir for SLURM job logs.
        
    Example:
        >>> run_fastp("/path/to/input", "/path/to/output")
        Submitted SLURM array job for 10 files using /abs/path/fastp.sh.
    """
    
    # Resolve absolute paths
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    script_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), script_path)
    )
    
    # Check that script exists
    if not os.path.isfile(script_path):
        raise FileNotFoundError(f"Script not found: {script_path}")
    
    # Grab all .fastq.gz files to determine array size
    files = sorted([f for f in os.listdir(input_dir) if f.endswith(".fastq.gz") or f.endswith(".fastq")])

    if not files:
        raise ValueError(f"No .fastq or .fastq.gz files found in {input_dir}")
    
    if not files:
        raise ValueError(f"No .fastq.gz files found in {input_dir}")
    
    array_range = f"1-{len(files)}"
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    
    # Submit the job with SLURM array option
    submit_cmd = [
        "sbatch",
        f"--array={array_range}",
        f"--output={log_dir}/fastp_%A_%a.out",
        script_path,
        input_dir,
        output_dir
    ]
    
    subprocess.run(submit_cmd)
    print(f"Submitted SLURM array job for {len(files)} files using {script_path}.")

def write_fastp_summary(output_dir, filename="fastp_summary.csv"):
    """Parse fastp logs and write summary table to logs directory."""

    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    log_files = [
        f for f in os.listdir(log_dir)
        if f.endswith(".out") and not f.lower().startswith("fastp_")
    ]

    samples = []
    reads_passed = []
    reads_filtered = []

    for f in tqdm(log_files):
        path = os.path.join(log_dir, f)

        with open(path) as fh:
            text = fh.read()

        sample = f.replace("_fastp_report.out", "").replace(".out","")

        match_passed = re.search(r"reads passed filter:\s+([\d,]+)", text)
        match_failed_low = re.search(r"reads failed due to low quality:\s+([\d,]+)", text)
        match_failed_N = re.search(r"reads failed due to too many N:\s+([\d,]+)", text)
        match_failed_short = re.search(r"reads failed due to too short:\s+([\d,]+)", text)

        passed = int(match_passed.group(1).replace(",", "")) if match_passed else 0

        failed = 0
        for m in [match_failed_low, match_failed_N, match_failed_short]:
            if m:
                failed += int(m.group(1).replace(",", ""))

        samples.append(sample)
        reads_passed.append(passed)
        reads_filtered.append(failed)

    df = pd.DataFrame({
        "sample": samples,
        "reads_passed": reads_passed,
        "reads_filtered": reads_filtered
    })

    df["total_reads"] = df["reads_passed"] + df["reads_filtered"]
    df["filtered_percent"] = df["reads_filtered"] / df["total_reads"] * 100
    df = df.sort_values("total_reads", ascending=False).reset_index(drop=True)

    output_path = os.path.join(log_dir, filename)
    df.to_csv(output_path, index=False)

    print(f"fastp summary written to: {output_path}")

def split_fastq(
    input_fastq: str,
    output_dir: str,
    chunks: list[int],
    jobs: int = 24,
):
    """Splits a FASTQ file into multiple chunk sets using fastqsplitter.

    Derives the output file base name from the input filename and runs
    fastqsplitter in parallel for each requested chunk count via a bash script.

    Args:
        input_fastq: Path to the input FASTQ file (.fastq.gz / .fq.gz).
        output_dir: Path to the output directory. Created if it doesn't exist.
        chunks: List of chunk counts to split into, e.g. [50, 100, 200].
        jobs: Number of parallel jobs to run. Defaults to 24.

    Raises:
        FileNotFoundError: If the input FASTQ file does not exist.
        FileNotFoundError: If the split_fastq.sh script cannot be found.
        subprocess.CalledProcessError: If the bash script exits with a non-zero status.
    """
    input_path = Path(input_fastq).resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"Input FASTQ file not found: {input_path}")

    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    # Strip known FASTQ extensions to derive a clean base name
    basename = input_path.name
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if basename.endswith(ext):
            basename = basename[: -len(ext)]
            break

    chunks_str = ",".join(str(c) for c in chunks)

    script_path = (Path(__file__).parent / "savio_jobs" / "split_fastq.sh").resolve()
    if not script_path.is_file():
        raise FileNotFoundError(f"split_fastq.sh not found at: {script_path}")

    cmd = [
        "bash", str(script_path),
        "-i", str(input_path),
        "-o", str(output_path),
        "-c", chunks_str,
        "-j", str(jobs),
    ]

    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print("Done!")