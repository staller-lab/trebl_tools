import sys

sys.path = [
    p for p in sys.path if "/.local/lib" not in p
]  # Use conda env installation of duckdb
print(sys.path)

import duckdb  # For connecting to your DuckDB database
import pandas as pd  # For DataFrame manipulation
import numpy as np  # For numerical operations (e.g., np.round, np.isfinite)
import seaborn as sns  # For plotting (barplot, styling)
import matplotlib
import os
import matplotlib.pyplot as plt
from trebl_tools import preprocess
from trebl_tools import finder
import tempfile
import os
import re
import shutil
import pathlib
import dask.dataframe as dd
from pathlib import Path
import subprocess

os.environ["MPLBACKEND"] = "Agg"
from tqdm import tqdm
from trebl_tools.preprocess import time_it
from trebl_tools import error_correct
import atexit


class UMIDeduplicator:
    """Unified UMI Deduplication class combining simple counts and full UMI-tools pipeline.

    This class provides comprehensive UMI (Unique Molecular Identifier) deduplication
    functionality, supporting both simple counting methods and advanced UMI-tools
    directional deduplication workflows.

    Args:
        db_path (str): Path to the DuckDB database file.
        bc_objects (list): List of barcode objects with name and length attributes.
        step_name (str): Name identifier for this processing step.
        descriptor (str): Additional descriptor for table naming.
        step1_map_name (str): Name of the step1 mapping table in the database.
        fastq_path (str): Path to the input FASTQ file.
        output_path (str): Directory path for output files.
        refined_map_suffix (str): Suffix for the final refined mapping table (last step from map refiner).

    Attributes:
        db_path (str): Path to DuckDB database.
        con (duckdb.DuckDBPyConnection): Database connection object.
        bc_objects (list): Barcode objects.
        cols (list): List of barcode column names.
        step_name (str): Processing step identifier.
        table_prefix (str): Generated prefix for database tables.
        base (str): Base filename derived from FASTQ path.

    Example:
        >>> bc_objs = [barcode1, barcode2]
        >>> dedup = UMIDeduplicator(
        ...     db_path="data.db",
        ...     bc_objects=bc_objs,
        ...     step_name="trebl_experiment",
        ...     descriptor="test",
        ...     step1_map_name="step1_map",
        ...     fastq_path="sample.fastq.gz",
        ...     output_path="results/",
        ...     refined_map_suffix="quality_designed"
        ... )
        >>> merged_df = dedup.run_both_deduplications()
    """

    def __init__(
        self,
        db_path,
        bc_objects,
        step_name,
        descriptor,
        step1_map_name,
        fastq_path,
        output_path,
        refined_map_suffix,
    ):
        self.db_path = db_path
        self.con = duckdb.connect(self.db_path)
        self.bc_objects = bc_objects
        self.cols = [bc.name for bc in bc_objects]

        self.step_name = step_name
        self.descriptor = descriptor
        self.step1_map_name = step1_map_name
        self.output_path = output_path

        self.refined_map_suffix = refined_map_suffix

        self.table_prefix = self.step_name + "_" + "_".join(self.cols) + "_"
        if descriptor:
            self.table_prefix += f"{descriptor}_"

        # UMI-tools specific attributes
        if fastq_path:
            path = Path(fastq_path)
            base = path.name
            for suffix in [".fastq", ".gz", ".assembled"]:
                base = base.replace(suffix, "")
            self.base = base
            self.fastq_path = fastq_path

    def counts_per_umi(self):
        """Return DataFrame with simple counts of each UMI per barcode combination.

        Returns:
            pd.DataFrame: DataFrame with columns for barcode(s), UMI, and read counts.
                Sorted by barcode combinations and read counts (descending).

        Example:
            >>> counts_df = dedup.counts_per_umi()
            >>> print(counts_df.head())
               ADBC2  HawkBCs      UMI  reads
            0  AAAA      TTTT  ACGTACGT     15
            1  AAAA      TTTT  TGCATGCA      8
        """
        select_cols_sql = ", ".join(self.cols)  # e.g., "ADBC2, HawkBCs"

        query = f"""
            SELECT {select_cols_sql}, UMI, COUNT(*) AS reads
            FROM {self.table_prefix}{self.refined_map_suffix}
            GROUP BY {select_cols_sql}, UMI
            ORDER BY {select_cols_sql}, reads DESC
        """

        df_counts = self.con.execute(query).fetchdf()
        return df_counts

    # --------- Methods from SimpleUMIDeduplicator ---------
    def unique_umis_per_barcodes(self):
        """Count unique UMIs per barcode combination and save as new table.

        Creates a new table in the database containing the count of distinct UMIs
        for each barcode combination.

        Note:
            This method creates a new table with suffix '_umis_collapsed' and
            sets the alias_name attribute to 'count' for the count column.

        Example:
            >>> dedup.unique_umis_per_barcodes()
            # Creates table: trebl_experiment_ADBC2_HawkBCs_test_quality_designed_umis_collapsed
        """

        # Include individual columns in SELECT and GROUP BY
        select_cols_sql = ", ".join(self.cols)

        self.new_table_name = (
            f"{self.table_prefix}{self.refined_map_suffix}_umis_collapsed"
        )

        # print(f"Counting unique UMIs per barcode(s) and saving as {self.new_table_name}...")

        # Build the alias
        # self.alias_name = f"{'_'.join(self.cols)}_umis_unique"
        self.alias_name = "count"

        query = f"""
            CREATE OR REPLACE TABLE {self.new_table_name} AS
            SELECT {select_cols_sql}, 
                   COUNT(DISTINCT UMI) AS {self.alias_name}
            FROM {self.table_prefix}{self.refined_map_suffix}
            GROUP BY {select_cols_sql}
            ORDER BY {self.alias_name} DESC
        """

        self.con.execute(query)

    def merge_simple_with_step1_map(self, save=True):
        """Merge simple UMI counts with step1 mapping table.

        Performs an inner join between the UMI count table and the step1 mapping
        table on barcode columns. Filters out quality and designed columns.

        Args:
            save (bool, optional): Whether to save results to CSV file.
                Defaults to True.

        Returns:
            pd.DataFrame: Merged DataFrame with UMI counts and step1 mapping data.
                Includes an 'info' column with the step_name value.

        Note:
            Removes columns named "Designed" and columns containing "_qual".
            CSV filename format: {base}_{barcode_cols}_umis_unique_with_step1_map.csv

        Example:
            >>> merged_df = dedup.merge_simple_with_step1_map(save=True)
            Saved to results/sample_ADBC2_HawkBCs_umis_unique_with_step1_map.csv
        """
        select_cols_sql = ", ".join(self.cols)  # e.g., AD_BC, RPTR_BC

        join_condition = " AND ".join([f"a.{col} = b.{col}" for col in self.cols])

        query = f"""
            SELECT b.*, a."{self.alias_name}"
            FROM "{self.new_table_name}" AS a
            INNER JOIN "{self.step1_map_name}" AS b
            ON {join_condition}
        """
        merged_df = self.con.execute(query).fetchdf()
        merged_df["info"] = self.step_name

        filtered_cols = [
            col
            for col in merged_df.columns
            if col != "Designed" and "_qual" not in col and col != "count"
        ]
        merged_df = merged_df[filtered_cols]

        if self.output_path and save == True:
            # Ensure the output directory exists
            os.makedirs(self.output_path, exist_ok=True)

            # Build a safe filename
            filename = os.path.join(
                self.output_path,
                f"{self.base}_{'_'.join(self.cols)}_umis_unique_with_step1_map.csv",
            )
            file_path = os.path.join(self.output_path, filename)
            merged_df.to_csv(file_path, index=False)
            print(f"Saved to {file_path}")

        return merged_df

    @time_it
    def save_simple_deduplication(self):
        """Save simple UMI deduplication results to TSV files.

        Saves two files:
        1. Simple UMI counts per barcode combination
        2. Read counts per individual UMI

        Returns:
            pd.DataFrame: DataFrame containing the simple UMI counts.

        Note:
            Files are saved with tab separation (.tsv format).
            Output filenames: {base}_simple_umi_counts.tsv and {base}_reads_per_umi.tsv
            Execution time is automatically logged via @time_it decorator.

        Example:
            >>> result_df = dedup.save_simple_deduplication()
            Saved to results/sample_simple_umi_counts.tsv
            Saved to results/sample_reads_per_umi.tsv
            Done in 2.34 seconds.
        """

        query = f"""
            SELECT * FROM "{self.new_table_name}"
        """
        merged_df = self.con.execute(query).fetchdf()

        if self.output_path:
            # Ensure the output directory exists
            os.makedirs(self.output_path, exist_ok=True)

            if self.output_path:
                output_dir = self.output_path
            else:
                output_dir = ""

            counts_per_bc = os.path.join(output_dir, f"{self.base}")
            merged_df.to_csv(
                f"{counts_per_bc}_simple_umi_counts.tsv", index=False, sep="\t"
            )
            print(f"Saved to {counts_per_bc}_simple_umi_counts.tsv")

            counts_per_umi_df = self.counts_per_umi()
            counts_per_umi_df.to_csv(
                f"{counts_per_bc}_reads_per_umi.tsv", index=False, sep="\t"
            )
            print(f"Saved to {counts_per_bc}_reads_per_umi.tsv")

        return merged_df

    def run_simple_deduplication(self):
        """Execute the simple UMI deduplication workflow.

        Runs the basic UMI counting approach by calling unique_umis_per_barcodes().
        This method only counts distinct UMIs per barcode combination without
        considering UMI similarity or directional information.

        Note:
            This is the faster, simpler alternative to the full UMI-tools pipeline.
            Results can be retrieved using merge_simple_with_step1_map() or
            save_simple_deduplication().
        """
        self.unique_umis_per_barcodes()
        # return self.merge_simple_with_step1_map()

    def show_tables(self):
        """Display all tables in the connected DuckDB database.

        Returns:
            pd.DataFrame: DataFrame listing all table names in the database.

        Example:
            >>> tables_df = dedup.show_tables()
            >>> print(tables_df)
                     name
            0  step1_map
            1  trebl_experiment_ADBC2_HawkBCs_quality_designed
        """
        return self.con.execute("SHOW TABLES").fetchdf()

    # --------- Methods from UMIToolsDeduplicator ---------
    @preprocess.time_it
    def generate_fastq(self, suffix="_umi_extracted.fastq"):
        """Generate FASTQ file with UMIs in headers and concatenated barcodes as sequences.
        Used for UMI-tools directional deduplication workflow.

        Creates a FASTQ file where:
        - Read ID contains the UMI sequence
        - Read sequence is concatenated barcode columns
        - Quality scores are placeholder 'I' characters

        Args:
            suffix (str, optional): Suffix for output FASTQ filename.
                Defaults to "_umi_extracted.fastq".

        Note:
            Only includes reads with valid UMI and barcode sequences.

        Example:
            >>> dedup.generate_fastq("_custom_suffix.fastq")
            Generating FASTQ with UMIs in header and barcodes as sequence...
            Writing FASTQ to results/sample_custom_suffix.fastq (1000 reads)...
            FASTQ complete: results/sample_custom_suffix.fastq
            Done in 1.23 seconds.
        """
        print("Generating FASTQ with UMIs in header and barcodes as sequence...")

        # Concatenate barcode columns
        if len(self.cols) == 1:
            concat_expr = self.cols[0]
        else:
            concat_expr = " || ".join(self.cols)  # no separator

        query = f"""
            SELECT UMI, {concat_expr} AS barcode_seq
            FROM {self.table_prefix}{self.refined_map_suffix}
        """
        result = self.con.execute(query).fetchall()
        n_rows = len(result)

        # display(result)

        # Prepare output path
        if self.output_path:
            output_file = os.path.join(self.output_path, f"{self.base}{suffix}")
            os.makedirs(self.output_path, exist_ok=True)
        else:
            output_file = f"{self.base}{suffix}"

        print(f"Writing FASTQ to {output_file} ({n_rows} reads)...")

        with open(output_file, "w") as f:
            for umi, seq in tqdm(result, total=n_rows, desc="Writing FASTQ"):
                # Only keep rows with both a barcode and umi
                if (
                    umi is not None
                    and seq is not None
                    and len(seq) > 1
                    and len(umi) > 1
                ):
                    header = f"@{umi}"
                    plus = "+"
                    qual = "I" * len(seq)
                    f.write(f"{header}\n{seq}\n{plus}\n{qual}\n")

        self.umi_fastq = output_file
        print(f"FASTQ complete: {self.umi_fastq}")

    @preprocess.time_it
    def generate_barcode_fasta_and_index(self, suffix="_barcodes"):
        """Generate FASTA reference of unique barcodes and create bowtie2 index.
        Used for UMI-tools directional deduplication workflow to align reads to barcode sequences.

        Creates a FASTA file containing all unique concatenated barcode sequences
        and builds a bowtie2 index for alignment. The FASTA file is removed after
        indexing to save space.

        Args:
            suffix (str, optional): Suffix for FASTA filename. Defaults to "_barcodes".

        Note:
            Loads bowtie module from Savio.
            Uses 32 threads for index building (BOWTIE2_INDEX_BUILDER_THREADS).

        Example:
            >>> dedup.generate_barcode_fasta_and_index("_ref")
            Saving unique barcode(s) as reference file...
            Creating table of unique concatenated barcodes: trebl_experiment_ADBC2_unique_barcodes
            Writing FASTA to results/sample_ref.fa...
            Indexing FASTA with bowtie2-build, prefix: results/sample_ref_index
            Done in 5.67 seconds.
        """
        print("Saving unique barcode(s) as reference file...")

        if len(self.cols) == 1:
            concat_expr = self.cols[0]
        else:
            concat_expr = " || ".join(self.cols)  # no separator

        new_table_name = f"{self.table_prefix}unique_barcodes"

        query = f"""
            CREATE OR REPLACE TABLE {new_table_name} AS
            SELECT DISTINCT {concat_expr} AS barcode
            FROM {self.table_prefix}{self.refined_map_suffix}
        """

        print(f"Creating table of unique concatenated barcodes: {new_table_name}")
        self.con.execute(query)

        if self.output_path:
            # Ensure output_path ends with a slash safely
            output_file = os.path.join(self.output_path, f"{self.base}{suffix}.fa")

        else:
            output_file = os.path.join(f"{self.base}{suffix}.fa")

        print(f"Writing FASTA to {output_file}...")

        # Fetch all barcodes from DuckDB
        result = self.con.execute(f"SELECT barcode FROM {new_table_name}").fetchall()
        with open(output_file, "w") as f:
            for i, (barcode,) in enumerate(result, start=1):
                # Use the barcode itself as the header
                f.write(f">{barcode}\n{barcode}\n")

        # Index the fasta file too
        bowtie2_index_prefix = output_file.replace(".fa", "_index")
        self.bowtie2_index_prefix = bowtie2_index_prefix

        print(f"Indexing FASTA with bowtie2-build, prefix: {bowtie2_index_prefix}")

        # Module system fallback
        cmd = [
            "bash",
            "-c",
            f"module load bio/bowtie2/2.5.1-gcc-11.4.0 && "
            f"export BOWTIE2_INDEX_BUILDER_THREADS=32 && "
            f"bowtie2-build {output_file} {bowtie2_index_prefix}",
        ]
        subprocess.run(cmd, check=True)

        # Optionally remove the FASTA after indexing
        os.remove(output_file)

    @preprocess.time_it
    def align_sort_and_deduplicate_umis(self):
        """Execute the complete UMI-tools alignment and deduplication pipeline."""
        # --- Setup output directory ---
        output_dir = self.output_path or os.getcwd()
        os.makedirs(output_dir, exist_ok=True)

        # --- Create temporary working directory ---
        tmp_dir = tempfile.mkdtemp(prefix="trebl_pipeline_")
        atexit.register(lambda: shutil.rmtree(tmp_dir, ignore_errors=True))

        # --- Paths for intermediate and final files ---
        sam_file = os.path.join(tmp_dir, f"{self.base}_umi_extracted.sam")
        bam_file = os.path.join(tmp_dir, f"{self.base}_umi_extracted.bam")
        sorted_bam = os.path.join(tmp_dir, f"{self.base}_umi_extracted.sorted.bam")
        dedup_bam = os.path.join(output_dir, f"{self.base}_umi_deduplicated.bam")
        counts_per_bc = os.path.join(
            output_dir, f"{self.base}_directional_umi_counts.tsv"
        )

        bc_len = sum(bc.length for bc in self.bc_objects)
        L = min(32, bc_len)

        # --- Load modules for Bowtie2 and Samtools ---
        module_cmd = (
            "module load bio/bowtie2/2.5.1-gcc-11.4.0 "
            "bio/samtools/1.17-gcc-11.4.0 && "
        )

        # --- Bowtie2 alignment ---
        print("Aligning .FASTQ to reference .FA ...")
        subprocess.run(
            module_cmd + f"bowtie2 -p 32 -x {self.bowtie2_index_prefix} "
            f"-U {self.umi_fastq} -S {sam_file} --norc --end-to-end "
            f"--very-sensitive -N 0 -L {L} -k 1 --score-min L,0,0",
            shell=True,
            check=True,
            executable="/bin/bash",
        )

        # --- Convert SAM -> BAM ---
        print("Converting SAM -> BAM ...")
        subprocess.run(
            module_cmd + f"samtools view -b {sam_file} -o {bam_file}",
            shell=True,
            check=True,
            executable="/bin/bash",
        )

        # --- Sort BAM ---
        print("Sorting BAM ...")
        subprocess.run(
            module_cmd + f"samtools sort -@ 32 -o {sorted_bam} {bam_file}",
            shell=True,
            check=True,
            executable="/bin/bash",
        )

        # --- Index sorted BAM ---
        print("Indexing BAM ...")
        subprocess.run(
            module_cmd + f"samtools index {sorted_bam}",
            shell=True,
            check=True,
            executable="/bin/bash",
        )

        # --- Deduplicate UMIs ---
        print("Deduplicating UMIs ...")
        subprocess.run(
            [
                "umi_tools",
                "dedup",
                "-I",
                sorted_bam,
                "-S",
                dedup_bam,
                "--method=directional",
                "--per-contig",
                "--per-gene",
                # add this if needed:
                # "--umi-tag=RX",
            ],
            check=True,
        )

        # --- Index deduplicated BAM (NOW it exists) ---
        print("Indexing dedup BAM ...")
        subprocess.run(
            module_cmd + f"samtools index {dedup_bam}",
            shell=True,
            check=True,
            executable="/bin/bash",
        )

        # --- Count UMIs per barcode ---
        print("Saving final counts ...")
        subprocess.run(
            [
                "umi_tools",
                "count",
                "-I",
                dedup_bam,
                "-S",
                counts_per_bc,
                "--per-contig",
            ],
            check=True,
        )

        print("UMI workflow complete!")

    def run_umi_tools_deduplication(self):
        """Execute the complete UMI-tools deduplication workflow.

        Runs the full UMI-tools pipeline by calling:
        1. generate_fastq() - Create FASTQ with UMIs in headers
        2. generate_barcode_fasta_and_index() - Create reference and index
        3. align_sort_and_deduplicate_umis() - Complete alignment and deduplication

        Note:
            This is the comprehensive, slower method that uses UMI-tools
            directional deduplication algorithm. Results in more accurate
            deduplication compared to simple counting methods.

        Example:
            >>> dedup.run_umi_tools_deduplication()
            # Executes complete pipeline automatically
        """
        self.generate_fastq()
        self.generate_barcode_fasta_and_index()
        self.align_sort_and_deduplicate_umis()

    def merge_complex_with_step1_map(self, save=True):
        """Merge UMI-tools directional deduplication results with step1 mapping table.

        Reads the UMI-tools output TSV file and merges it with the step1 mapping
        table on barcode columns. Renames columns to follow naming conventions.

        Args:
            save (bool, optional): Whether to save merged results to CSV.
                Defaults to True.

        Returns:
            pd.DataFrame: Merged DataFrame with directional UMI counts and
                step1 mapping information.

        Note:
            Expects UMI-tools output file: {base}_directional_umi_counts.tsv

        Example:
            >>> merged_df = dedup.merge_complex_with_step1_map(save=True)
            Saved to results/sample_ADBC2_HawkBCs_umis_directional_deduplic_with_step1_map.csv
        """
        if self.output_path:
            output_dir = self.output_path
        else:
            output_dir = ""

        counts_file = os.path.join(
            output_dir, f"{self.base}_directional_umi_counts.tsv"
        )
        final_counts = pd.read_csv(counts_file, sep="\t")

        # Rename columns to match barcode naming convention
        final_counts = final_counts.rename(
            columns={
                "gene": "_".join(self.cols),
                "count": f"{'_'.join(self.cols)}_umis_directional_deduplic",
            }
        )

        # Read step1 map from DuckDB
        step1_df = self.con.execute(f'SELECT * FROM "{self.step1_map_name}"').fetchdf()

        # Merge on barcode columns
        merge_cols = [
            col
            for col in self.cols
            if col in final_counts.columns and col in step1_df.columns
        ]
        merged_df = step1_df.merge(final_counts, on=merge_cols, how="inner")

        # Optional: filter out unwanted columns like before
        filtered_cols = [
            col
            for col in merged_df.columns
            if col != "Designed" and "_qual" not in col and col != "count"
        ]
        merged_df = merged_df[filtered_cols]

        # Add info column
        merged_df["info"] = self.step_name

        # Save CSV if output_path specified
        if self.output_path and save == True:
            os.makedirs(self.output_path, exist_ok=True)
            filename = (
                f"{self.base}_{'_'.join(self.cols)}_umis_unique_with_step1_map.csv"
            )
            file_path = os.path.join(self.output_path, filename)
            merged_df.to_csv(file_path, index=False)
            print(f"Saved to {file_path}")

        return merged_df

    def run_both_deduplications(self):
        """Execute both simple and UMI-tools deduplication methods and merge results.

        Runs both deduplication approaches:
        1. Simple UMI counting (unique_umis_per_barcodes)
        2. UMI-tools directional deduplication (run_umi_tools_deduplication)

        Then merges both results with the step1 mapping table if available.

        Returns:
            pd.DataFrame or None: If step1_map_name is provided, returns merged
                DataFrame with both deduplication results. Otherwise, saves
                simple deduplication results only and returns None.

        Note:
            This method provides comprehensive UMI analysis by comparing both
            approaches. The merged output contains results from both methods
            for comparison and validation.
            Output filename: {base}_{barcode_cols}_umis_deduplic_with_step1_map.csv

        Example:
            >>> merged_df = dedup.run_both_deduplications()
            Starting simple deduplication.
            Finished simple deduplication.

            Starting UMI Tools directional deduplication.
            # ... UMI-tools pipeline output ...
            Finished UMI Tools directional deduplication.

            Saved to results/sample_ADBC2_HawkBCs_umis_deduplic_with_step1_map.csv
        """

        print("Starting simple deduplication.")
        self.unique_umis_per_barcodes()
        print("Finished simple deduplication.\n")

        print("Starting UMI Tools directional deduplication.")
        self.run_umi_tools_deduplication()
        print("Finished UMI Tools directional deduplication.\n")

        if self.step1_map_name:
            df1 = self.merge_simple_with_step1_map(save=False)
            df2 = self.merge_complex_with_step1_map(save=False)

            merged_df = pd.merge(df1, df2, how="outer")

            if self.output_path:
                os.makedirs(self.output_path, exist_ok=True)
                filename = os.path.join(
                    self.output_path,
                    f"{self.base}_{'_'.join(self.cols)}_umis_deduplic_with_step1_map.csv",
                )
                file_path = os.path.join(self.output_path, filename)
                merged_df.to_csv(file_path, index=False)
                print(f"Saved to {file_path}")

            return merged_df
        else:
            self.save_simple_deduplication()
            # Need to save the simple deduplication result
