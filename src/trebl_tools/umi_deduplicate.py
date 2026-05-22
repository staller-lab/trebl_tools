import sys

sys.path = [
    p for p in sys.path if "/.local/lib" not in p
]  # Use conda env installation of duckdb

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


# Get the directory of the current Python executable
conda_bin = os.path.dirname(sys.executable)
# Prepend it to PATH so subprocess calls see it
os.environ["PATH"] = conda_bin + os.pathsep + os.environ.get("PATH", "")


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
            FROM "{self.table_prefix}{self.refined_map_suffix}_umis_collapsed" AS a
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
            SELECT * FROM "{self.table_prefix}{self.refined_map_suffix}_umis_collapsed"
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
    def run_umi_tools_deduplication(self):
        """
        Run umi_tools count_tab on preprocessed DuckDB table instead of FASTQ/BAM.
        Exports TSV and produces normalized directional counts output with barcode columns.
        """
        import subprocess
        import pandas as pd
    
        os.makedirs(self.output_path, exist_ok=True)
    
        tsv_path = os.path.join(self.output_path, f"{self.base}_umi_count_input.tsv")
        concat_cols = " || ".join(self.cols)
    
        self.con.execute(f"""
            COPY (
                SELECT '_' || UMI AS umi_col,
                       {concat_cols} AS gene
                FROM {self.table_prefix}{self.refined_map_suffix}
                ORDER BY {concat_cols}
            ) TO '{tsv_path}'
            WITH (FORMAT 'csv', DELIMITER E'\t', HEADER FALSE)
        """)
        print(f"Exported TSV for count_tab to {tsv_path}")
    
        output_path = os.path.join(self.output_path, f"{self.base}_directional_umi_counts.tsv")
        cmd = [
            "umi_tools",
            "count_tab",
            "--stdin", tsv_path,
            "--stdout", output_path
        ]
        subprocess.run(cmd, check=True)
        self.count_tab_output = output_path
        print(f"Count_tab output saved to {output_path}")
    
        # Normalize output: split gene back into barcode columns
        directional_df = pd.read_csv(output_path, sep="\t")
    
        if "gene" in directional_df.columns:
            start = 0
            for bc in self.bc_objects:
                end = start + bc.length
                directional_df[bc.name] = directional_df["gene"].astype(str).str.slice(start, end)
                start = end
    
        # remaining_cols = [c for c in directional_df.columns if c not in ordered_cols]
        # directional_df = directional_df[ordered_cols + remaining_cols]
        directional_df = directional_df.drop(columns = "gene")
    
        directional_df.to_csv(output_path, sep="\t", index=False)
        print(f"Normalized directional output saved to {output_path}")
    
        if os.path.exists(tsv_path):
            os.remove(tsv_path)
            print(f"Deleted intermediate file {tsv_path}")

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
