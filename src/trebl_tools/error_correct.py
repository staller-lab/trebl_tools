import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from trebl_tools.preprocess import time_it
from pathlib import Path
import pandas as pd
import matplotlib
import os
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
from Bio.Seq import Seq
import shutil

def rc_and_swap(preceder: str, post: str, length: int) -> (str, str):
    """Generate reverse complement patterns for barcode extraction.
    
    Transforms preceder and post patterns for reverse complement mode by
    reverse complementing and swapping their positions. This is used when
    sequences need to be processed in reverse orientation.
    
    Args:
        preceder (str): Original DNA sequence pattern that precedes the barcode.
        post (str): Original DNA sequence pattern that follows the barcode.
        length (int): Length parameter (currently unused but maintained for compatibility).
        
    Returns:
        tuple[str, str]: A tuple containing:
            - preceder_rc (str): Reverse complement of the original post sequence
            - post_rc (str): Reverse complement of the original preceder sequence
            
    Example:
        >>> preceder, post = "ATCG", "GCTA"
        >>> rc_preceder, rc_post = rc_and_swap(preceder, post, 20)
        >>> print(f"RC preceder: {rc_preceder}, RC post: {rc_post}")
        RC preceder: TAGC, RC post: CGAT
        
    Note:
        The swap occurs because in reverse complement orientation, the roles
        of preceder and post patterns are reversed.
    """
    preceder_rc = str( Seq(post).reverse_complement() )
    post_rc     = str( Seq(preceder).reverse_complement() )
    return preceder_rc, post_rc
    
@time_it
def concatenate_fastqs(input_fastq_list, output_dir):
    """Concatenate multiple FASTQ files into a single temporary file.
    
    Combines multiple FASTQ files into one temporary file for batch processing
    with UMI-tools. Handles both single file inputs and lists of files.
    
    Args:
        input_fastq_list (str or list[str]): Path to single FASTQ file or list
            of FASTQ file paths to concatenate.
        output_dir (str or Path): Directory where temporary concatenated file
            will be created.
            
    Returns:
        str: Path to the concatenated FASTQ file. If input was a single file,
            returns the original file path unchanged.
            
    Note:
        Creates a temporary file named "tmp_combined.fastq" in the output directory.
        
    Example:
        >>> files = ["sample1.fastq", "sample2.fastq", "sample3.fastq"]
        >>> combined_path = concatenate_fastqs(files, "temp_dir/")
        >>> print(f"Combined file: {combined_path}")
        Combined file: temp_dir/tmp_combined.fastq
        Done in 5.23 seconds.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_fastq = output_dir / "tmp_combined.fastq"
    
    # If only one file, just return it
    if isinstance(input_fastq_list, str):
        return input_fastq_list
    
    input_fastq_list = [str(Path(f).resolve()) for f in input_fastq_list]
    with open(tmp_fastq, "wb") as outfile:
        for f in input_fastq_list:
            with open(f, "rb") as infile:
                outfile.write(infile.read())
    return str(tmp_fastq)
    
def run_whitelist_on_concat_domains(fastq_path, output_dir, set_cell_number=None, prefix=None):
    """Execute UMI-tools whitelist on concatenated barcode sequences.

    Runs UMI-tools whitelist algorithm on a FASTQ file containing extracted
    and concatenated barcode sequences to identify canonical sequences and
    group similar sequences for error correction.

    Args:
        fastq_path (str or Path): Path to FASTQ file containing extracted barcodes.
        output_dir (str or Path): Directory to save whitelist outputs including
            log files, whitelist text files, and diagnostic plots.
        set_cell_number (int, optional): Force a specific number of expected
            cell barcodes instead of using automatic knee detection. Defaults to None.
        prefix (str, optional): Prefix for output file naming. If None, uses
            the FASTQ filename stem. Defaults to None.

    Returns:
        dict: Dictionary with paths to generated files:
            - "log": Path to the whitelist log file
            - "whitelist": Path to the whitelist text output file  
            - "plot_prefix": Prefix for generated diagnostic plots

    Raises:
        subprocess.CalledProcessError: If UMI-tools whitelist command fails.

    Example:
        >>> result = run_whitelist_on_concat_domains(
        ...     "barcodes.fastq", 
        ...     "output/", 
        ...     set_cell_number=1000,
        ...     prefix="sample1"
        ... )
        Running umi_tools whitelist on barcodes.fastq ...
        Whitelist complete.
        - Log: output/sample1_whitelist.log
        - Output: output/sample1_whitelist.txt
        - Plots: output/sample1_plots_*.png

    Note:
        Uses density-based knee detection method by default. The barcode pattern
        "(?P<umi_1>.{1})(?P<cell_1>.*)" extracts the entire sequence as cell barcode
        with a single dummy UMI base. Requires UMI-tools conda environment.
    """
    fastq_path = Path(fastq_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True)

    if prefix is None:
        prefix = fastq_path.stem  # e.g. "HAR_AD_filtered_barcode"

    log_path = output_dir / f"{prefix}_whitelist.log"
    whitelist_path = output_dir / f"{prefix}_whitelist.txt"
    plot_prefix = output_dir / f"{prefix}_plots"

    umi_tools_exe = shutil.which("umi_tools")
    if umi_tools_exe is None:
        raise RuntimeError("umi_tools not found on PATH. Activate the appropriate conda environment.")

    cmd = [
        umi_tools_exe,
        "whitelist",
        "--knee-method=density",
        "--stdin", str(fastq_path),
        "--bc-pattern", "(?P<umi_1>.{1})(?P<cell_1>.*)",
        "--extract-method=regex",
        "--method=reads",
        "--log", str(log_path),
        "--stdout", str(whitelist_path),
        "--plot-prefix", str(plot_prefix)
    ]

    if set_cell_number is not None:
        print("Using custom cell number.")
        cmd.extend(["--set-cell-number", str(set_cell_number)])

    print(f"Running umi_tools whitelist on {fastq_path.name} ...")
    subprocess.run(cmd, check=True)
    
    print(f"Whitelist complete.\n- Log: {log_path}\n- Output: {whitelist_path}\n- Plots: {plot_prefix}_*.png")

    return {
        "log": log_path,
        "whitelist": whitelist_path,
        "plot_prefix": plot_prefix
    }

def convert_txt_to_whitelist_mapping_df_from_path(whitelist_path):
    """Convert UMI-tools whitelist output to sequence mapping DataFrame.

    Reads a UMI-tools whitelist text file and creates a mapping DataFrame
    that shows the relationship between original sequences and their
    canonical (corrected) forms for downstream error correction.

    Args:
        whitelist_path (str or Path): Path to the whitelist text file generated
            by UMI-tools. Expected format: tab-separated with columns
            [canonical, collapsed, largest_count, counts].

    Returns:
        pd.DataFrame: Mapping DataFrame with columns:
            - "original": Original sequence (canonical or collapsed)
            - "canonical": Canonical sequence this original maps to

    Example:
        >>> mapping_df = convert_txt_to_whitelist_mapping_df_from_path("whitelist.txt")
        >>> print(mapping_df.head())
               original    canonical
        0  ATCGATCGATCG  ATCGATCGATCG
        1  ATCGATCGATCT  ATCGATCGATCG  
        2  ATCGATCGATTG  ATCGATCGATCG
        
    Note:
        Each canonical sequence maps to itself, and all sequences listed in the
        "collapsed" column map to their corresponding canonical sequence.
        Currently hardcoded to not perform reverse complement transformation.
    """
    whitelist_path = Path(whitelist_path)
    
    # Assume umi_tools output has canonical, collapsed, largest_count, counts
    wl_df = pd.read_csv(
        whitelist_path, sep="\t", header=None,
        names=["canonical", "collapsed", "largest_count", "counts"]
    )

    def maybe_revcomp(seq):
        return seq
        
    all_rows = []
    for _, row in wl_df.iterrows():
        canonical = maybe_revcomp(row["canonical"])
        all_rows.append({"original": canonical, "canonical": canonical})
        
        if pd.notna(row["collapsed"]):
            for c in row["collapsed"].split(","):
                all_rows.append({"original": maybe_revcomp(c), "canonical": canonical})

    mapping_df = pd.DataFrame(all_rows)
    return mapping_df
