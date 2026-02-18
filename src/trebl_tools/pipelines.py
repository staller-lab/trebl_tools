import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path = [
    p for p in sys.path if "/.local/lib" not in p
]  # Use conda env installation of duckdb

import duckdb
import os
import pandas as pd
from trebl_tools import initial_map, map_refiner, complexity, finder, preprocess, error_correct, plotting, umi_deduplicate

import math
import seaborn as sns
import re
from tqdm import tqdm
import re
import sys  
import json
import subprocess
from pathlib import Path


# Defines the sequence of refinement operations for each step.
# Indexed by step name and whether error correction is enabled.
MAP_ORDERS = {
    "step1": {
        False: ['grouped','thresholded','barcode_exists','unique_target','quality','designed'],
        True: ['barcode_exists', 'quality', 'error_corrected','grouped','thresholded','unique_target','designed'],
    },
    "step2": {
        False: ['grouped','thresholded','barcode_exists', 'quality','designed'],
        True: ['barcode_exists', 'quality','error_corrected','grouped','thresholded','designed'],
    },
}


class TreblPipeline:
    """
    End-to-end TREBL analysis pipeline orchestrator.

    This class coordinates the complete TREBL workflow including initial barcode
    mapping, quality filtering, error correction, UMI deduplication, and activity
    score calculations. Supports both simple validation workflows and complex
    experimental analyses with multiple samples and replicates.

    Args:
        db_path (str): Path to the DuckDB database file. Will be created if
            it doesn't exist.
        design_file_path (str or None): Path to the design validation file. 
            If None, design-based filtering is skipped throughout the pipeline.
        error_correction (bool, optional): Whether to enable UMI-tools based
            error correction during refinement steps. Defaults to False.
        output_path (str or Path, optional): Output directory for results,
            figures, and intermediate files. If None, results are not saved
            to disk. Defaults to None.
        test_n_reads (int or bool, optional): If set to an integer, limits
            initial mapping to the first N reads for rapid testing/debugging.
            If False, processes all reads. Defaults to False.

    Example:
        >>> # Initialize pipeline for full experiment
        >>> pipeline = TreblPipeline(
        ...     db_path="experiment.duckdb",
        ...     design_file_path="designs.csv",
        ...     error_correction=True,
        ...     output_path="results/"
        ... )
        >>> 
        >>> # Run Step 1: Initial barcode mapping
        >>> step1_df = pipeline.run_step_1(
        ...     seq_file="step1_reads.fastq",
        ...     bc_objects=[ad_bc, rep_bc],
        ...     column_pairs=[("REP_BC", "AD")],
        ...     reads_threshold=10
        ... )
        >>> 
        >>> # Analyze full experiment with UMI deduplication
        >>> results = pipeline.trebl_experiment_analysis(
        ...     AD_seq_files=["ad_sample1.fastq", "ad_sample2.fastq"],
        ...     AD_bc_objects=[ad_bc, adbc_bc],
        ...     RT_seq_files=["rt_sample1.fastq", "rt_sample2.fastq"], 
        ...     RT_bc_objects=[rt_bc],
        ...     reverse_complement=True,
        ...     AD_umi_object=umi_obj,
        ...     RT_umi_object=umi_obj
        ... )
        >>> 
        >>> # Calculate final activity scores
        >>> mean_activity, summed_activity = pipeline.calculate_activity_scores(
        ...     step1_path="results/step1.csv",
        ...     AD_bc_objects=[ad_bc, adbc_bc],
        ...     RT_bc_objects=[rt_bc],
        ...     time_regex=r"(\d+)h",
        ...     rep_regex=r"rep(\d+)"
        ... )

    Note:
        The pipeline automatically handles table naming, intermediate file management,
        and result aggregation across multiple samples. All processing steps are
        logged and intermediate results are saved for debugging and validation.
    """

    def __init__(
        self,
        db_path,
        design_file_path,
        error_correction=False,
        output_path=None,
        test_n_reads=False
    ):
        """
        Initialize a TREBL pipeline run.

        Sets up database connection, output directories, and processing configuration.
        Creates necessary directory structure for results and figures if output_path
        is provided.

        Args:
            db_path (str): Path to the DuckDB database file. Will be created if
                it doesn't exist.
            design_file_path (str or None): Path to the design validation file. 
                If None, design-based filtering is skipped throughout the pipeline.
            error_correction (bool, optional): Whether to enable UMI-tools based
                error correction during refinement steps. Defaults to False.
            output_path (str or Path, optional): Output directory for results,
                figures, and intermediate files. If None, results are not saved
                to disk. Defaults to None.
            test_n_reads (int or bool, optional): If set to an integer, limits
                initial mapping to the first N reads for rapid testing/debugging.
                If False, processes all reads. Defaults to False.

        Note:
            When output_path is provided, creates subdirectories:
            - {output_path}/figures/ for all plots and visualizations
            - {output_path}/{sample_name}/ for per-sample intermediate files
        """
        self.con = duckdb.connect(db_path)
        self.db_path = db_path
        self.error_correction = error_correction
        self.design_file_path = design_file_path
        self.test_n_reads = test_n_reads
        
        if output_path:
            self.output_path = Path(output_path)
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_figures_path = self.output_path / "figures"
            self.output_figures_path.mkdir(exist_ok=True)
        else:
            self.output_path = None
            self.output_figures_path = None
            
    def _run_initial_mappers(self, mapper_specs):
        """
        Execute one or more InitialMapper instances with specified configurations.

        Internal helper method that creates and runs InitialMapper objects based
        on the provided specifications. Automatically handles test mode if enabled.

        Args:
            mapper_specs (list[dict]): List of mapper specifications. Each dict
                must contain:
                    - seq_file (str): FASTQ file path for barcode extraction
                    - step_name (str): Step name prefix for DuckDB table naming
                    - bc_objects (list): Barcode objects defining extraction patterns
                    - reverse_complement (bool): Whether to reverse complement reads
                    - design_file_path (str or None): Design validation file path

        Note:
            Respects the test_n_reads setting for rapid testing. Each mapper
            creates its own initial mapping table in the connected database.
        """
        for spec in mapper_specs:
            mapper = initial_map.InitialMapper(
                db_path=self.db_path,
                seq_file=spec["seq_file"],
                step_name=spec["step_name"],
                bc_objects=spec["bc_objects"],
                reverse_complement=spec["reverse_complement"],
                design_file_path=spec["design_file_path"],
            )
            if self.test_n_reads:
                mapper.create_test_map(test_n_reads=self.test_n_reads)
            else:
                mapper.create_map()

    def _run_refiners(self, refiners, plot_titles):
        """
        Execute MapRefiner instances and optionally generate loss plots.

        Internal helper method that runs refinement pipelines and creates
        standardized loss visualization plots if output path is configured.

        Args:
            refiners (list[MapRefiner]): List of MapRefiner instances to execute.
                Each refiner should be fully configured with processing parameters.
            plot_titles (list[str or None]): Titles for generated loss plots.
                Must match the length of refiners list. Use None to skip plotting
                for a particular refiner.

        Returns:
            list[MapRefiner]: The same refiner objects after execution, allowing
                access to results and generated tables.

        Note:
            Loss plots show processing efficiency across refinement steps and
            are saved automatically if output_figures_path is configured.
        """
        for refiner, plot_title in zip(refiners, plot_titles):
            refiner.refine_map_from_db()
    
            if self.output_figures_path and plot_titles:
                refiner.plot_loss(text_offset = -0.2)
                if plot_title:
                    plt.title(plot_title)
    
        return refiners

    def reads_distribution(self, seq_file, bc_objects, step_name, reverse_complement):
        """
        Visualize read-count distributions for threshold determination.

        Performs initial mapping and grouping for a single FASTQ file and generates
        a histogram showing the distribution of read counts per barcode combination.
        This analysis helps determine appropriate read count thresholds for filtering.

        Args:
            seq_file (str): Path to FASTQ file containing sequences for analysis.
            bc_objects (list[Barcode]): Barcode objects defining extraction patterns
                and expected lengths for the library.
            step_name (str): Step identifier used for DuckDB table naming and
                plot titles.
            reverse_complement (bool): Whether sequences should be reverse
                complemented prior to barcode extraction.

        Note:
            Creates grouped barcode combinations table and generates a read count
            histogram. Useful for determining reads_threshold parameters before
            running full refinement pipelines.

        Example:
            >>> pipeline.reads_distribution(
            ...     seq_file="test_reads.fastq",
            ...     bc_objects=[ad_bc, rep_bc], 
            ...     step_name="threshold_test",
            ...     reverse_complement=True
            ... )
            # Generates histogram showing read count distribution
        """
        
        self._run_initial_mappers([
            {
                "seq_file": seq_file,
                "step_name": step_name,
                "bc_objects": bc_objects,
                "reverse_complement": reverse_complement,
                "design_file_path": self.design_file_path,
            }
        ])

        # Create refiner
        refiner = map_refiner.MapRefiner(
            db_path=self.db_path,
            bc_objects=bc_objects,
            column_pairs=[],
            reads_threshold=1,
            map_order=["grouped"],
            step_name=step_name,
            output_figures_path=self.output_figures_path,
            design_file = self.design_file_path
        )

        refiners = self._run_refiners([refiner], plot_titles=[f"{step_name} grouped"])
        
    def step1_reads_distribution(self, seq_file, bc_objects, reverse_complement, step_suffix=""):
        """
        Convenience wrapper for Step 1 read distributions.

        Args:
            seq_file (str): FASTQ file path.
            bc_objects (list): Barcode objects.
            reverse_complement (bool): Whether to reverse complement reads.
            step_suffix (str, optional): Suffix appended to the step name.
        """
        step_name = "step1" + step_suffix
        self.reads_distribution(seq_file, bc_objects, step_name, reverse_complement)

    def step2_reads_distribution(self, 
                                 AD_seq_file,
        AD_bc_objects,
        RT_seq_file,
        RT_bc_objects,
        reverse_complement,
        step_suffix=""):
        """
        Convenience wrapper for Step 2 AD and RT libraries.

        Args:
            AD_seq_file (str): FASTQ file for AD reads.
            AD_bc_objects (list): AD barcode objects.
            RT_seq_file (str): FASTQ file for RT reads.
            RT_bc_objects (list): RT barcode objects.
            reverse_complement (bool): Whether to reverse complement reads.
            step_suffix (str, optional): Suffix appended to the step name.
        """

        step_name = "step2" + step_suffix
        
        self.reads_distribution(AD_seq_file, AD_bc_objects, step_name, reverse_complement)
        self.reads_distribution(RT_seq_file, RT_bc_objects, step_name, reverse_complement)

    def trebl_experiment_reads_distribution(self, AD_seq_files,
        AD_bc_objects,
        RT_seq_files,
        RT_bc_objects,
        reverse_complement,
        step_suffix=""):
        """
        Convenience wrapper for TREBL Experiment AD and RT libraries.

        Args:
            AD_seq_file (str): FASTQ file for AD reads.
            AD_bc_objects (list): AD barcode objects.
            RT_seq_file (str): FASTQ file for RT reads.
            RT_bc_objects (list): RT barcode objects.
            reverse_complement (bool): Whether to reverse complement reads.
            step_suffix (str, optional): Suffix appended to the step name.
        """

        step_name_base = "trebl_experiment" + step_suffix
        
        for file_path in AD_seq_files:
            base_name = os.path.basename(file_path)
            name_only = base_name.split('.')[0]
            step_name = f"{step_name_base}_{name_only}"
            self.reads_distribution(file_path, AD_bc_objects, step_name, reverse_complement)
            
        for file_path in RT_seq_files:
            base_name = os.path.basename(file_path)
            name_only = base_name.split('.')[0]
            step_name = f"{step_name_base}_{name_only}"
            self.reads_distribution(file_path, RT_bc_objects, step_name, reverse_complement)
            
    def run_step_1(
        self,
        seq_file,
        bc_objects,
        column_pairs,
        reads_threshold,
        reverse_complement=False,
        step_suffix = ""
    ):
        """
        Run TREBL Step 1: map reads to designed AD barcodes.

        Performs initial mapping, refinement, optional error correction,
        and outputs the final designed map.

        Args:
            seq_file (str): FASTQ file containing Step 1 reads.
            bc_objects (list): Barcode objects defining extraction rules.
            column_pairs (list[tuple]): Column-pair checks used
                to remove barcode collisions. 

                Each tuple has the form:

                    (key_column(s), target_column(s))

                where `key_column(s)` and `target_column(s)` can be either
                a single column name (str) or a tuple/list of column names.

                The constraint enforces that, for each **key**, at least
                90% of its reads must map to a single **target**. If a key
                value maps to multiple targets without one reaching the 90%
                threshold, the key is considered ambiguous and discarded.

                **Single-column example**:

                    column_pairs = [("RPTR_BC", "AD")]

                This ensures that each reporter barcode (``RPTR_BC``) maps
                predominantly to a single AD. If a reporter barcode has
                reads mapping to multiple ADs without one exceeding 90%,
                that reporter barcode is removed.

                **Multi-column example**:

                    column_pairs = [
                        (("RPTR_BC"), ("Hawk_BC", "AD_BC"))
                    ]

                This checks that each RPTR BC
                maps to a single Hawkins AD barcode and AD barcode combination 
                with ≥90% of reads. Ambiguous key combinations are removed.

                Multiple constraints may be applied sequentially:

                    column_pairs = [
                        ("AD_BC", "AD"),
                        ("RPTR_BC", "AD"),
                        (("AD_BC", "RPTR_BC"), ("AD", "SampleID")),
                    ]

            reads_threshold (int): Minimum number of reads required for a
                barcode or barcode pair to be kept.
            reverse_complement (bool, optional): Whether reads should be
                reverse complemented prior to barcode extraction.
                Defaults to False.
            step_suffix (str, optional): Suffix appended to the step name
                for DuckDB table naming. This is useful when you want
                to distinguish multiple runs or subsets of the same step.

                Example:
                    "_spike_in" — if processing spike-in samples separately
                    from your main dataset.
                    "_new_data" — for data from new sequencing run.

        Returns:
            pd.DataFrame: Final designed Step 1 mapping after all refinement
            steps have been applied.
        """


        step_name = "step1" + step_suffix

        
        # Initial mapping
        self._run_initial_mappers([
            {
                "seq_file": seq_file,
                "step_name": step_name,
                "bc_objects": bc_objects,
                "reverse_complement": reverse_complement,
                "design_file_path": self.design_file_path,
            }
        ])

        manual_ec_threshold = True
        if self.error_correction and reads_threshold == 1:
                manual_ec_threshold = False

        # Create refiner
        refiner = map_refiner.MapRefiner(
            db_path=self.db_path,
            bc_objects=bc_objects,
            column_pairs=column_pairs,
            reads_threshold=reads_threshold,
            map_order=MAP_ORDERS["step1"][self.error_correction],
            step_name=step_name,
            output_figures_path=self.output_figures_path,
            manual_ec_threshold = manual_ec_threshold,
            design_file = self.design_file_path
        )
    
        # Use _run_refiners to handle refinement and plotting
        refiners = self._run_refiners([refiner], plot_titles=["Step 1"])
    
        # Get the dataframe
        df = refiners[0].get_map_df('designed')
        if self.output_path:
            df.to_csv(self.output_path / f"{step_name}.csv", index=False)

        self.step1_df = df
        
        return df

    def run_step_2(
        self,
        AD_seq_file,
        AD_bc_objects,
        RT_seq_file,
        RT_bc_objects,
        reverse_complement,
        reads_threshold_AD,
        reads_threshold_RT,
        step1_map_csv_path,  # Updated argument to accept CSV path
        step_suffix="",
    ):
        """
        Run TREBL Step 2: Analyze intermediate complexity.

        Performs separate refinement for AD and RT libraries and computes
        overlap with the Step 1 map.

        Args:
            AD_seq_file (str): FASTQ file containing AD reads.
            AD_bc_objects (list): Barcode objects for ADs.
            RT_seq_file (str): FASTQ file containing RT reads.
            RT_bc_objects (list): Barcode objects for reporters.
            reverse_complement (bool): Whether to reverse complement reads.
            reads_threshold_AD (int): Minimum reads per AD barcode.
            reads_threshold_RT (int): Minimum reads per RT barcode.
            step1_map_csv_path (str): Path to the Step 1 map CSV file.
            step_suffix (str, optional): Suffix appended to the step name.

        Returns:
            dict: Dictionary with keys:
                - "AD_step2": AD DataFrame
                - "RT_step2": RT DataFrame
                - "step1_overlap": Overlap statistics
        """
        step_name = "step2" + step_suffix

        # Initial mapping
        self._run_initial_mappers([
            {
                "seq_file": AD_seq_file,
                "step_name": step_name,
                "bc_objects": AD_bc_objects,
                "reverse_complement": reverse_complement,
                "design_file_path": self.design_file_path,  # AD uses design
            },
            {
                "seq_file": RT_seq_file,
                "step_name": step_name,
                "bc_objects": RT_bc_objects,
                "reverse_complement": reverse_complement,
                "design_file_path": None,  # RT skips design
            },
        ])

        refiners = []

        for bc_objs, reads_threshold in zip(
            (AD_bc_objects, RT_bc_objects),
            (reads_threshold_AD, reads_threshold_RT)
        ):
            step2_map_order = MAP_ORDERS["step2"][self.error_correction].copy()

            refiners.append(
                map_refiner.MapRefiner(
                    db_path=self.db_path,
                    bc_objects=bc_objs,
                    column_pairs=[],
                    reads_threshold=reads_threshold,
                    map_order=step2_map_order,
                    step_name=step_name,
                    design_file=self.design_file_path if bc_objs is AD_bc_objects else None,
                    output_figures_path=self.output_figures_path,
                )
            )

        # Run refiners
        refiners = self._run_refiners(
            refiners,
            plot_titles=["AD Step 2", "RT Step 2"],
        )

        # Load Step 1 map CSV into DuckDB
        step1_map_name = "step1_map_temp"
        self.con.execute(f"DROP TABLE IF EXISTS {step1_map_name}")
        self.con.execute(f"""
            CREATE TABLE {step1_map_name} AS
            SELECT * FROM read_csv_auto('{step1_map_csv_path}')
        """)

        # Compute AD–reporter complexity and overlap
        checker = complexity.ComplexityChecker(
            db_path=self.db_path,
            step_name=step_name,
            step1_map_name=step1_map_name,  # Use the temporary table name
            step_suffix="designed",
            barcode_groups=[AD_bc_objects, RT_bc_objects],
        )

        overlap = checker.count_overlap()

        AD_df = refiners[0].get_map_df('designed')
        RT_df = refiners[1].get_map_df('designed')

        # Save CSVs
        if self.output_path:
            AD_df.to_csv(self.output_path / f"{step_name}_AD.csv", index=False)
            RT_df.to_csv(self.output_path / f"{step_name}_RT.csv", index=False)

        return {
            "AD_step2": AD_df,
            "RT_step2": RT_df,
            "step1_overlap": overlap
        }
    
    def _duckdb_safe_name(self, base_name):
        """
        Convert a filename to a DuckDB-safe table name.
        
        Replaces periods, hyphens, and spaces with underscores, removes special
        characters, and ensures the name doesn't start with a digit.
        
        Args:
            base_name (str): Original filename or identifier.
            
        Returns:
            str: DuckDB-safe table name.
        """
        name = base_name.replace(".", "_")    # replace periods with underscores
        name = base_name.replace("-", "_").replace(" ", "_")
        name = re.sub(r'[^0-9a-zA-Z_]', '_', name)
        if re.match(r'^\d', name):
            name = f"_{name}"
        return name
    

    def _run_trebl_experiment_helper(
        self,
        seq_files,
        bc_objects,
        reverse_complement,
        reads_threshold=1,
        umi_object=None,
        step_name_suffix="",
        count_col_name=None,
        gene_col_name=None,
        concat_gene=False,
        umi_deduplication='both'
    ):
        """
        Core TREBL experiment runner.

        Handles both UMI and non-UMI workflows. If a UMI object is provided,
        UMI-based deduplication is applied. The workflow supports two deduplication
        modes: 'simple' (single-step deduplication) or 'both' (simple + directional/complex deduplication).

        Args:
            seq_files (list[str]): FASTQ file paths to process.
            bc_objects (list): Barcode objects describing library barcodes.
            reverse_complement (bool): Whether to reverse complement reads before processing.
            reads_threshold (int, optional): Minimum reads per barcode/UMI to retain. Defaults to 1.
            umi_object (optional): UMI configuration object. If provided, triggers UMI-based workflow.
            step_name_suffix (str, optional): Suffix appended to DuckDB step names. Defaults to "".
            count_col_name (str, optional): Column name for UMI counts in final merged results.
            gene_col_name (str, optional): Column name for gene/barcode identifiers.
            concat_gene (bool, optional): If True, concatenates barcode columns to form gene identifiers. Defaults to False.
            umi_deduplication (str, optional): Deduplication mode. Options:
                - 'simple': Only simple UMI deduplication.
                - 'both' (default): Runs both simple and directional/complex deduplication.

        Returns:
            pd.DataFrame: Aggregated experiment results.
                - If UMI workflow: Merged DataFrame containing both "simple" and "complex" UMI counts.
                - If non-UMI workflow: Barcode-level counts for the experiment.

        Notes:
            - When `umi_object` is provided, directional/complex UMI counts are stored in 
            "{output_path}/{sample}_directional_umi_counts.tsv" and simple UMI counts in 
            "{output_path}/{sample}_simple_umi_counts.tsv".
            - Deduplicator also outputs reads-per-UMI summary files for quality control.
            - Non-UMI workflow skips UMI deduplication and uses barcode grouping and thresholding.
            - Error correction steps (if enabled) are applied before deduplication.
        """
        
        step_name_prefix = "trebl_experiment_" + step_name_suffix
        
        results = []
        simple_results = []
    
        for file_path in seq_files:
            base_name = os.path.basename(file_path)
            name_only = base_name.split('.')[0]
            name_only = self._duckdb_safe_name(name_only)
            step_name = f"{step_name_prefix}{name_only}"
    
            output_dir = self.output_path / step_name
            output_dir.mkdir(parents=True, exist_ok=True)
    
            design_file_path = self.design_file_path if "AD" in [bc.name for bc in bc_objects] else None
    
            # Initial mapping
            mapper_kwargs = dict(
                db_path=self.db_path,
                step_name=step_name,
                seq_file=file_path,
                bc_objects=bc_objects,
                reverse_complement=reverse_complement,
                design_file_path=design_file_path
            )
            if umi_object:
                mapper_kwargs["umi_object"] = umi_object
    
            mapper = initial_map.InitialMapper(**mapper_kwargs)
            mapper.create_map()
    
            # Refinement
            manual_ec_threshold = not (self.error_correction and reads_threshold == 1)
            map_order = ["quality", "error_corrected"] if self.error_correction else ["quality"]
            if not umi_object:
                # non-UMI workflow adds additional steps
                map_order = map_order + (["grouped", "thresholded", "designed"] if not self.error_correction else ["grouped", "thresholded", "designed"])
    
            refiner = map_refiner.MapRefiner(
                db_path=self.db_path,
                bc_objects=bc_objects,
                column_pairs=[],
                reads_threshold=reads_threshold,
                map_order=map_order,
                step_name=step_name,
                output_figures_path=output_dir,
                manual_ec_threshold=manual_ec_threshold,
            )
            refiner.refine_map_from_db()
            refiner.plot_loss()
            if self.error_correction:
                refiner.plot_error_correction()
    
            if umi_object:
                # Deduplication
                refined_map_suffix = "error_corrected" if self.error_correction else "quality"
                deduplicator = umi_deduplicate.UMIDeduplicator(
                    db_path=self.db_path,
                    bc_objects=bc_objects,
                    step_name=step_name,
                    descriptor="",
                    step1_map_name=None,
                    fastq_path=file_path,
                    output_path=output_dir,
                    refined_map_suffix=refined_map_suffix,
                )
                
                if umi_deduplication == 'simple':
                    deduplicator.run_simple_deduplication()
                    deduplicator.save_simple_deduplication()
                else:
                    deduplicator.run_both_deduplications()

                    # Load results
                    complex_df = pd.read_csv(output_dir / f"{name_only}_directional_umi_counts.tsv", sep="\t")
                    complex_df["name"] = name_only
                    results.append(complex_df)

                simple_df = pd.read_csv(output_dir / f"{name_only}_simple_umi_counts.tsv", sep="\t")
                simple_df["name"] = name_only
                simple_results.append(simple_df)

                # Reads per UMI
                one_file_reads_per_UMI = deduplicator.counts_per_umi()
                one_file_reads_per_UMI["name"] = name_only
                one_file_reads_per_UMI.to_csv(output_dir / f"{name_only}_reads_per_umi.tsv", sep="\t")

        # Merge results
        if umi_object:
            # Merge results safely for both simple and complex dedup
            
            # Complex DF may be empty if only simple dedup was run
            if results:
                complex_df = pd.concat(results, ignore_index=True).rename(
                    columns={"count": count_col_name, "gene": gene_col_name}
                )
            else:
                complex_df = pd.DataFrame(columns=[gene_col_name, count_col_name, "name"])
            
            # Simple DF should always exist
            simple_df = pd.concat(simple_results, ignore_index=True).rename(
                columns={"count": count_col_name}
            )
            
            # Concatenate barcode columns if requested
            if concat_gene:
                concat_cols = [bc.name for bc in bc_objects]
                simple_df[gene_col_name] = simple_df[concat_cols].agg("".join, axis=1)
            
            # Merge simple and complex DFs
            if not complex_df.empty:
                merged = pd.merge(
                    complex_df,
                    simple_df,
                    on=[gene_col_name, "name"],
                    suffixes=("_complex", "_simple"),
                    how="outer"
                )
            else:
                # If no complex DF, just return simple_df with renamed columns
                merged = simple_df.copy()
                merged = merged.rename(columns={count_col_name: f"{count_col_name}_simple"})
            
            return merged
        else:
             # Non-UMI workflow: grab all relevant tables
            tables = refiner.show_tables()
            first_bc_name = bc_objects[0].name
            for table in tables:
                if step_name_prefix in table[0] and first_bc_name in table[0]:
                    df = refiner.get_map_df(table[0])
                    df["sample"] = table[0][len(step_name_prefix):]
                    results.append(df)

            return pd.concat(results, ignore_index=True)

    def trebl_experiment_analysis(
        self,
        AD_seq_files,
        AD_bc_objects,
        RT_seq_files,
        RT_bc_objects,
        reverse_complement,
        step1_map_csv_path=None,
        AD_umi_object=None,
        RT_umi_object=None,
        reads_threshold_AD=1,
        reads_threshold_RT=1,
        step_name_suffix="",
        umi_deduplication='both'
    ):
        """
        Run TREBL experiment analysis for both AD and RT libraries.

        The workflow automatically selects between a UMI or non-UMI
        pipeline depending on whether a UMI object is provided. If UMI deduplication
        is enabled, results from simple and/or directional/complex deduplication are merged.

        Args:
            AD_seq_files (list[str]): Paths to FASTQ files for AD library reads.
            AD_bc_objects (list): Barcode objects for AD library extraction.
            RT_seq_files (list[str]): Paths to FASTQ files for reporter (RT) reads.
            RT_bc_objects (list): Barcode objects for reporter library extraction.
            reverse_complement (bool): Whether reads should be reverse complemented prior to barcode extraction.
            step1_map_csv_path (str, optional): Path to Step 1 map CSV for computing overlap plots.
            AD_umi_object (optional): UMI object for AD library. If provided, triggers UMI deduplication.
            RT_umi_object (optional): UMI object for RT library. If provided, triggers UMI deduplication.
            reads_threshold_AD (int, optional): Minimum reads per AD barcode to retain. Defaults to 1.
            reads_threshold_RT (int, optional): Minimum reads per RT barcode to retain. Defaults to 1.
            step_name_suffix (str, optional): Suffix for DuckDB table and output names.
            umi_deduplication (str, optional): Deduplication mode for both AD and RT libraries.
                Options:
                    - 'simple': Only simple UMI deduplication is applied.
                    - 'both' (default): Both simple and directional/complex deduplication are performed.

        Returns:
            dict: Dictionary containing final TREBL experiment results:
                - "AD_results" (pd.DataFrame): AD library results with merged UMI counts if UMI workflow.
                - "RT_results" (pd.DataFrame): RT library results with merged UMI counts if UMI workflow.

        Notes:
            - UMI-based workflows produce two count tables per sample: 
            simple UMI counts and directional/complex UMI counts. These are merged in the final output.
            - Non-UMI workflows return barcode counts after grouping, thresholding, and optional error correction.
            - If `output_path` is set, results are saved as CSV:
                "{output_path}/AD_trebl_experiment_results.csv" and
                "{output_path}/RT_trebl_experiment_results.csv".
            - Barcode quality/loss plots are generated for both AD and RT libraries.
        """
        
        step_name_prefix = "trebl_experiment_" + step_name_suffix
        
        experiments = {
            "AD": {
                "seq_files": AD_seq_files,
                "bc_objects": AD_bc_objects,
                "umi_object": AD_umi_object,
                "reads_threshold": reads_threshold_AD,
                "count_col_name": "AD_umi_count",
                "gene_col_name": "AD_ADBC_concat",
                "concat_gene": True,
                "output_file": "ADBC_trebl_experiment_results.csv",
                "umi_deduplication" : umi_deduplication
            },
            "RT": {
                "seq_files": RT_seq_files,
                "bc_objects": RT_bc_objects,
                "umi_object": RT_umi_object,
                "reads_threshold": reads_threshold_RT,
                "count_col_name": "RTBC_umi_count",
                "gene_col_name": "RTBC",
                "concat_gene": False,
                "output_file": "RTBC_trebl_experiment_results.csv",
                "umi_deduplication": umi_deduplication
            },
        }

        results = {}

        for name, spec in experiments.items():
            df = self._run_trebl_experiment_helper(
                seq_files=spec["seq_files"],
                bc_objects=spec["bc_objects"],
                reverse_complement=reverse_complement,
                reads_threshold=spec["reads_threshold"],
                umi_object=spec["umi_object"],
                count_col_name=spec.get("count_col_name"),
                gene_col_name=spec.get("gene_col_name"),
                concat_gene=spec.get("concat_gene", False),
                step_name_suffix=step_name_suffix,
                umi_deduplication=spec["umi_deduplication"]
            )

            if self.output_path:
                df.to_csv(self.output_path / f"{name}_trebl_experiment_results.csv", index=False)

            results[f"{name}_results"] = df

        if step1_map_csv_path:
            # Load Step 1 map CSV into DuckDB as a temporary table
            step1_map_name = "step1_map_temp"
            self.con.execute(f"DROP TABLE IF EXISTS {step1_map_name}")
            self.con.execute(f"""
                CREATE TABLE {step1_map_name} AS
                SELECT * FROM read_csv_auto('{step1_map_csv_path}')
            """)

            # Plot loss for AD and RT
            self.plot_trebl_experiment_loss(
                AD_bc_objects, 
                step1_map_name, 
                step_name_prefix=step_name_prefix
            )
            self.plot_trebl_experiment_loss(
                RT_bc_objects, 
                step1_map_name, 
                step_name_prefix=step_name_prefix
            )
        
        return results

    def plot_trebl_experiment_loss(self, bc_objects, step1_map_name=None, step_name_prefix="trebl_experiment_"):
        """
        Plot barcode quality and mapping loss for a TREBL experiment.

        Generates bar plots showing:
            1. Total number of reads in each initial mapping table.
            2. Number of reads passing barcode quality checks (`_qual` columns).
            3. Number of reads that match the Step 1 map for all barcodes.

        Plots are saved as a PNG file in `self.output_figures_path` and 
        returned as Matplotlib figure and axes objects.

        Args:
            bc_objects (list): List of barcode objects to evaluate.
            step1_map_name (str): Name of the Step 1 DuckDB table used to
                calculate overlap with mapped barcodes.
            step_name_prefix (str, optional): Prefix used to identify TREBL
                experiment tables in DuckDB. Defaults to "trebl_experiment_".

        Returns:
            tuple: (fig, axes)
                - fig (matplotlib.figure.Figure): Figure object containing all subplots.
                - axes (list[matplotlib.axes._subplots.AxesSubplot]): Flattened list of subplot axes.
                  Returns None, None if no matching tables are found.

        """
        
        # Connect to DuckDB
        con = duckdb.connect(self.db_path)
    
        # Get all tables matching step_name_prefix and bc_object names
        tables = con.execute("SHOW TABLES").fetchall()
        bc_names = [bc.name for bc in bc_objects]
    
        result_prefixes = [
            table[0] for table in tables
            if step_name_prefix in table[0] 
            and any(bc_name in table[0] for bc_name in bc_names)
            and 'initial' in table[0]
        ]
        
                
        if not result_prefixes:
            print("No matching tables found.")
            return None, None
    
        num_plots = len(result_prefixes)
        max_cols = 5
        ncols = min(max_cols, num_plots)
        nrows = math.ceil(num_plots / ncols)
    
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*3), dpi=300, sharey=True, sharex=True)
        axes = axes.flatten() if num_plots > 1 else [axes]
        
        bc_names_str = "_".join([bc.name for bc in bc_objects])

        for i, file_name in tqdm(enumerate(result_prefixes), total=len(result_prefixes), desc="Plotting BCs"):
            # Total count
            total_count = con.execute(f'SELECT COUNT(*) FROM "{file_name}"').fetchone()[0]
    
            # Count rows where all BCs are qualified (if the _qual column exists)
            qual_cols = [f'"{bc.name}_qual"' for bc in bc_objects]

            if qual_cols:
                qual_conditions = " AND ".join(qual_cols)
                both_true = con.execute(f'SELECT COUNT(*) FROM "{file_name}" WHERE {qual_conditions}').fetchone()[0]
            else:
                both_true = 0
    
            # Count rows present in Step1 map for all BCs
            join_conditions = " AND ".join([f'm."{bc.name}" = s."{bc.name}"' for bc in bc_objects])

            if step1_map_name== None:
                step1_count = 0
            else:
                step1_count = con.execute(f'''
                    SELECT COUNT(*) 
                    FROM "{file_name}" AS m
                    JOIN "{step1_map_name}" AS s
                    ON {join_conditions}
                ''').fetchone()[0]
        
            # Prepare DataFrame for plotting
            plot_counts = pd.DataFrame({
                "Category": ["Total", "BC\nQual", "In\nStep1"],
                "Count": [total_count, both_true, step1_count]
            })
    
            ax = axes[i]
            sns.barplot(x="Category", y="Count", data=plot_counts, ax=ax, palette=["gray", "green", "blue"])
            for container in ax.containers:
                ax.bar_label(container, fmt='%d', label_type='edge', fontsize='small', padding=2)

            title_regex = f'{step_name_prefix}(.*)_{bc_names_str}_initial'
            match = re.search(title_regex, file_name)
            group_name = match.group(1) if match else file_name
            ax.set_title(str(group_name), fontsize='medium', y = 1.1)
            ax.set_xlabel("")
            ax.set_ylabel("")
    
        # Hide unused axes
        for j in range(num_plots, len(axes)):
            axes[j].set_visible(False)
    
        fig.supylabel("Count")
        sns.despine()
        plt.tight_layout(pad=1)

        fig.suptitle("Trebl Experiment Loss")
    
        save_path = os.path.join(self.output_figures_path, f"{step_name_prefix}{bc_names_str}_bc_quality_loss.png")
        plt.savefig(save_path)

        con.close()

    
        return fig, axes

    def calculate_activity_scores(
        self,
        step1_path,
        AD_bc_objects,
        RT_bc_objects,
        time_regex,
        rep_regex
    ):
        """
        Compute activity scores from AD and RT experiment results using Step 1 mapping.

        Calculates two complementary activity metrics by combining AD library UMI counts
        (input) with RT library UMI counts (output) using the Step 1 barcode-to-gene
        mapping as the linking table. Returns a consolidated table with both per-barcode
        averaged activity and summed activity calculations.

        Args:
            step1_path (str): Path to the Step 1 mapping CSV file containing the
                canonical barcode-to-gene relationships.
            AD_bc_objects (list[Barcode]): List of AD barcode objects with `.name`
                attributes corresponding to columns in the Step 1 mapping.
            RT_bc_objects (list[Barcode]): List of RT barcode objects with `.name`
                attributes for reporter identification.
            time_regex (str): Regular expression to extract timepoint values from
                sample names. Should contain a capture group for the time value.
                Example: r"(\d+)h" for "2h", r"t(\d+)" for "t24".
            rep_regex (str): Regular expression to extract replicate identifiers
                from sample names. Should contain a capture group for replicate number.
                Example: r"rep(\d+)" for "rep1", r"r(\d+)" for "r2".

        Returns:
            pd.DataFrame: Consolidated activity scores table with hierarchical columns:
                - Index: (AD, replicate) multi-index
                - Columns: Multi-level structure with:
                    - Level 0: timepoint values (e.g., 0, 2, 24)
                    - Level 1: metric types ('mean_activity', 'std_activity', 'summed_activity')
                
                Example structure:
                                        0                    2                    24
                            mean  std  sum      mean  std  sum      mean  std  sum
                AD    rep                                                            
                GENE1  1      0.45  0.12  0.52    1.23  0.34  1.15    2.10  0.67  1.98
                    2      0.48  0.15  0.55    1.18  0.29  1.12    2.05  0.71  1.95

        Note:
            **Activity Score Calculations:**

            1. **Averaged Activity (Per-Barcode)**:
                - For each barcode: activity = RT_UMIs / AD_UMIs
                - Mean and std calculated per (AD, time, rep) group
                - Results show variability across barcodes within each gene

            2. **Summed Activity (Per-Gene)**:
                - For each (AD, time, rep): sum all AD_UMIs and RT_UMIs across barcodes
                - activity = sum(RT_UMIs) / sum(AD_UMIs) 
                - Results show overall gene-level activity

            **Output Files (if output_path configured):**
            - `unaggregated_activities.csv`: Raw per-barcode activity scores
            - `consolidated_activity_scores.csv`: Combined table with all metrics

            **Data Integration Workflow:**
            1. Load AD and RT experiment results from previous analysis
            2. Extract time/replicate metadata using provided regex patterns
            3. Merge with Step 1 mapping to link AD and RT measurements
            4. Calculate per-barcode activity ratios
            5. Aggregate using both averaging and summing strategies
            6. Consolidate into single multi-level DataFrame

        Example:
            >>> activity_df = pipeline.calculate_activity_scores(
            ...     step1_path="results/step1.csv",
            ...     AD_bc_objects=[ad_bc, adbc_bc], 
            ...     RT_bc_objects=[rt_bc],
            ...     time_regex=r"(\d+)h",          # Extract "24" from "24h"
            ...     rep_regex=r"rep(\d+)"          # Extract "1" from "rep1"  
            ... )
            >>> print("Activity table shape:", activity_df.shape)
            >>> print("Available timepoints:", activity_df.columns.get_level_values(0).unique())
            >>> print("Available metrics:", activity_df.columns.get_level_values(1).unique())
            >>> 
            >>> # Access specific metric for all timepoints
            >>> mean_activities = activity_df.xs('mean_activity', level=1, axis=1)
            >>> summed_activities = activity_df.xs('summed_activity', level=1, axis=1)

        Raises:
            ValueError: If regex patterns fail to match sample names
            FileNotFoundError: If required input CSV files are not found
            KeyError: If expected columns are missing from input files
        """

        import sys
        import pandas as pd

        def extract_with_regex(series, regex, group=1, column_name=""):
            """Extract integer values using regex with validation."""
            try:
                extracted = series.str.extract(regex).iloc[:, group - 1]
                if extracted.isnull().any():
                    raise ValueError(
                        f"Regex failed to match all values in column '{column_name}'."
                    )
                return extracted.astype(int)
            except Exception:
                print(f"Error extracting '{column_name}' with regex '{regex}'.")
                print(series.head(10))
                sys.exit(1)

        # Load experiment results
        ad_results_path = self.output_path / "AD_trebl_experiment_results.csv"
        rt_results_path = self.output_path / "RT_trebl_experiment_results.csv"

        ad_column_names = [bc.name for bc in AD_bc_objects]
        rt_column_names = [bc.name for bc in RT_bc_objects]

        ad_bc_results = pd.read_csv(ad_results_path, index_col=0)
        ad_bc_results["time"] = extract_with_regex(
            ad_bc_results["name"], time_regex, column_name="time"
        )
        ad_bc_results["rep"] = extract_with_regex(
            ad_bc_results["name"], rep_regex, column_name="rep"
        )
        ad_bc_results = ad_bc_results[
            ad_column_names + ["time", "rep", "AD_umi_count_simple"]
        ].reset_index(drop=True)

        rt_bc_results = pd.read_csv(rt_results_path)
        rt_bc_results["time"] = extract_with_regex(
            rt_bc_results["name"], time_regex, column_name="time"
        )
        rt_bc_results["rep"] = extract_with_regex(
            rt_bc_results["name"], rep_regex, column_name="rep"
        )
        rt_bc_results = rt_bc_results[
            rt_column_names + ["time", "rep", "RTBC_umi_count_simple"]
        ].reset_index(drop=True)

        # Load Step 1 mapping
        step1_map = pd.read_csv(step1_path)
        step1_map = step1_map[ad_column_names + rt_column_names + ["AD"]]

        step1_map_with_ad = pd.merge(step1_map, ad_bc_results)
        step1_map_with_rt = pd.merge(step1_map, rt_bc_results)
        step1_map_with_ad_rt = pd.merge(
            step1_map_with_ad, step1_map_with_rt, how="outer"
        )

        step1_map_with_ad_rt["AD_umi_count_simple"] = (
            step1_map_with_ad_rt["AD_umi_count_simple"].fillna(0)
        )
        step1_map_with_ad_rt["RTBC_umi_count_simple"] = (
            step1_map_with_ad_rt["RTBC_umi_count_simple"].fillna(0)
        )

        # Per-barcode activity
        step1_map_with_ad_rt["activity"] = (
            step1_map_with_ad_rt["RTBC_umi_count_simple"]
            / step1_map_with_ad_rt["AD_umi_count_simple"]
        )

        if self.output_path:
            step1_map_with_ad_rt.to_csv(
                self.output_path / "bc_activities.csv", index=False
            )

        # Calculate both averaged and summed activities
        # Averaged activity across barcodes
        grouped_avg = (
            step1_map_with_ad_rt
            .groupby(["AD", "time", "rep"])["activity"]
            .agg(["mean", "std"])
            .reset_index()
        )

        # Summed activity across barcodes
        grouped_sum = (
            step1_map_with_ad_rt
            .groupby(["AD", "time", "rep"])
            .agg(
                summed_AD_UMIs=("AD_umi_count_simple", "sum"),
                summed_RT_UMIs=("RTBC_umi_count_simple", "sum"),
            )
            .reset_index()
        )
        
        grouped_sum["summed_activity"] = (
            grouped_sum["summed_RT_UMIs"] / grouped_sum["summed_AD_UMIs"]
        )

        # Pivot each metric separately
        pivoted_mean = grouped_avg.pivot(
            index=["AD", "rep"], columns="time", values="mean"
        )
        pivoted_std = grouped_avg.pivot(
            index=["AD", "rep"], columns="time", values="std"
        )
        pivoted_summed = grouped_sum.pivot(
            index=["AD", "rep"], columns="time", values="summed_activity"
        )

        # Create consolidated multi-level column DataFrame
        # First, ensure all DataFrames have the same columns (timepoints)
        all_timepoints = sorted(set(pivoted_mean.columns) | set(pivoted_std.columns) | set(pivoted_summed.columns))
        
        # Reindex all DataFrames to have the same columns
        pivoted_mean = pivoted_mean.reindex(columns=all_timepoints)
        pivoted_std = pivoted_std.reindex(columns=all_timepoints)
        pivoted_summed = pivoted_summed.reindex(columns=all_timepoints)
        
        # Create multi-level columns
        consolidated_data = {}
        for timepoint in all_timepoints:
            consolidated_data[(timepoint, 'bc_activity_avg')] = pivoted_mean[timepoint]
            consolidated_data[(timepoint, 'bc_activity_std')] = pivoted_std[timepoint]
            consolidated_data[(timepoint, 'pooled_activity')] = pivoted_summed[timepoint]
        
        # Create the consolidated DataFrame
        consolidated_df = pd.DataFrame(consolidated_data)
        
        # Sort columns by timepoint first, then by metric
        consolidated_df = consolidated_df.sort_index(axis=1)
        
        # Set proper column names
        consolidated_df.columns.names = ['timepoint', 'metric']
        
        if self.output_path:
            consolidated_df.to_csv(self.output_path / "AD_activities.csv")

        return consolidated_df
