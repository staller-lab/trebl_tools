import duckdb                # For connecting to your DuckDB database
import pandas as pd          # For DataFrame manipulation
import numpy as np           # For numerical operations (e.g., np.round, np.isfinite)
import seaborn as sns        # For plotting (barplot, styling)
import matplotlib
import os
import matplotlib.pyplot as plt
import tempfile
import os
import shutil
import pathlib
import dask.dataframe as dd

from trebl_tools import preprocess
from trebl_tools import finder

class ComplexityChecker:
    """Analyze barcode complexity and coverage between different processing steps.
    
    This class compares barcode sequences between the initial mapping table
    (designed sequences) and processed datasets to assess library coverage,
    identify missing sequences, and evaluate processing efficiency.
    
    Args:
        db_path (str): Path to the DuckDB database containing mapping tables.
        step_name (str): Identifier for the processing step being analyzed.
        step1_map_name (str): Name of the initial mapping table in the database
            containing the designed barcode sequences.
        step_suffix (str): Suffix of the processed table to compare against
            the initial mapping (e.g., "quality_designed").
        barcode_groups (list): List of barcode groups to analyze. Each group
            can be either a single Barcode object or a list of Barcode objects
            for multi-barcode analyses.
            
    Example:
        >>> # Single barcode analysis
        >>> bc1 = Barcode(name="AD_BC", preceder="ATCG", post="GCTA", length=20)
        >>> bc2 = Barcode(name="REP_BC", preceder="TTAA", post="GGCC", length=15)
        >>> 
        >>> # Multi-barcode combination analysis  
        >>> checker = ComplexityChecker(
        ...     db_path="experiment.duckdb",
        ...     step_name="step2", 
        ...     step1_map_name="step1_AD_ADBC_RTBC_quality_designed",
        ...     step_suffix="quality_designed",
        ...     barcode_groups=[bc1, bc2, [bc1, bc2]]  # Individual + combined
        ... )
        >>> 
        >>> # Analyze coverage
        >>> coverage_df = checker.count_overlap()
        >>> print(coverage_df)
              BC_type  map_unique  step2  seen_in_both  percent_of_map_seen
        0       AD_BC        1000    800           750                 75.0
        1      REP_BC         500    400           350                 70.0
        2  AD_BC,REP_BC       1500   1200          1100                 73.3
        
    Note:
        Requires that both the initial mapping table and processed tables
        exist in the connected DuckDB database. Barcode objects must have
        a 'name' attribute for SQL column identification.
    """
    
    def __init__(
        self,
        db_path,
        step_name,
        step1_map_name,
        step_suffix,
        barcode_groups
    ):
        """Initialize ComplexityChecker with database connection and analysis parameters.
        
        Args:
            db_path (str): Path to DuckDB database file.
            step_name (str): Processing step identifier.
            step1_map_name (str): Initial mapping table name.
            step_suffix (str): Processed table suffix.
            barcode_groups (list): List of barcode groups for analysis.
        """
        self.con = duckdb.connect(db_path)
        self.step_name = step_name
        self.step1_map_name = step1_map_name
        self.step_suffix = step_suffix
        self.barcode_groups = barcode_groups
        
    def show_tables(self):
        """Display all tables in the connected DuckDB database.
        
        Prints a list of all table names in the database for reference
        and debugging purposes.
        
        Example:
            >>> checker.show_tables()
            ('step1_map',)
            ('step2_AD_BC_quality_designed',)
            ('step2_REP_BC_quality_designed',)
        """
        tables = self.con.execute("SHOW TABLES").fetchall()
        for table in tables:
            print(table)
            
    def count_overlap_one_barcode_group(self, barcode_group):
        bc_columns = [barcode.name for barcode in barcode_group]
        table_prefix = self.step_name + "_" + "_".join(bc_columns) + "_"

        # Combine statistics for all barcodes
        map_column_list = ", ".join([f"m.{bc}" for bc in bc_columns])    
        step_column_list = ", ".join([f"s.{bc}" for bc in bc_columns])    
        where_clause = " AND ".join([f"m.{bc} = s.{bc}" for bc in bc_columns])

        pair_query = f"""
            SELECT
                '{",".join(bc_columns)}' AS BC_type,
                (SELECT COUNT(*) FROM (SELECT DISTINCT {map_column_list} FROM {self.step1_map_name} m) t) AS map_unique,
                (SELECT COUNT(*) FROM (SELECT DISTINCT {step_column_list} FROM {table_prefix}{self.step_suffix} s) t1) AS {self.step_name},
                (SELECT COUNT(*) FROM (SELECT DISTINCT {map_column_list} FROM {self.step1_map_name} m JOIN {table_prefix}{self.step_suffix} s ON {where_clause}) t2) AS seen_in_both
        """
        df = self.con.execute(pair_query).df()
        
    
        return df

    def count_overlap(self):
        """Calculate overlap statistics for a single barcode group.
        
        Compares unique barcode combinations between the initial mapping table
        and a processed table to determine coverage and overlap metrics.
        
        Args:
            barcode_group (list[Barcode]): List of Barcode objects to analyze
                as a group. For single barcode analysis, pass a list with one element.
                
        Returns:
            pd.DataFrame: DataFrame with a single row containing:
                - BC_type (str): Comma-separated list of barcode names in the group
                - map_unique (int): Number of unique barcode combinations in initial mapping
                - {step_name} (int): Number of unique combinations in processed table
                - seen_in_both (int): Number of combinations found in both tables
                
        Example:
            >>> bc_group = [ad_bc, rep_bc]  # Two-barcode combination
            >>> overlap_df = checker.count_overlap_one_barcode_group(bc_group)
            >>> print(overlap_df)
              BC_type  map_unique  step2  seen_in_both
            0  AD_BC,REP_BC    1500    1200          1100
            
        Note:
            Constructs table names using the pattern: {step_name}_{barcode_names}_{step_suffix}
            Uses DISTINCT counts to avoid duplicate sequences in the analysis.
        """
        results = []
        for barcode_group in self.barcode_groups:
            # If not iterable (or is a string), wrap in a list
            if not isinstance(barcode_group, (list, tuple, set)):
                barcode_group = [barcode_group]

            df = self.count_overlap_one_barcode_group(barcode_group)
            results.append(df)
    
        summary_df = pd.concat(results, ignore_index=True)
        summary_df["percent_of_map_seen"] = (
            100 * np.round(summary_df["seen_in_both"] / summary_df["map_unique"], 5)
        )
        return summary_df


