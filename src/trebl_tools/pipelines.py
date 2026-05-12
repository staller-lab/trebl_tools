import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path = [
    p for p in sys.path if "/.local/lib" not in p
]  # Use conda env installation of duckdb

import duckdb
import os
import pandas as pd
from trebl_tools import (
    initial_map,
    map_refiner,
    complexity,
    finder,
    preprocess,
    error_correct,
    plotting,
    umi_deduplicate,
)

import math
import seaborn as sns
import re
from tqdm import tqdm
import re
import sys
import json
import subprocess
from pathlib import Path

from IPython.display import display

# Defines the sequence of refinement operations for each step.
# Indexed by step name and whether error correction is enabled.
MAP_ORDERS = {
    "step1": {
        False: [
            "grouped",
            "thresholded",
            "barcode_exists",
            "unique_target",
            "quality",
            "designed",
        ],
        True: [
            "barcode_exists",
            "quality",
            "error_corrected",
            "grouped",
            "thresholded",
            "unique_target",
            "designed",
        ],
    },
    "step2": {
        False: ["grouped", "thresholded", "barcode_exists", "quality", "designed"],
        True: [
            "barcode_exists",
            "quality",
            "error_corrected",
            "grouped",
            "thresholded",
            "designed",
        ],
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
        test_n_reads=False,
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
                refiner.plot_loss(text_offset=-0.2)
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

        self._run_initial_mappers(
            [
                {
                    "seq_file": seq_file,
                    "step_name": step_name,
                    "bc_objects": bc_objects,
                    "reverse_complement": reverse_complement,
                    "design_file_path": self.design_file_path,
                }
            ]
        )

        # Create refiner
        refiner = map_refiner.MapRefiner(
            db_path=self.db_path,
            bc_objects=bc_objects,
            column_pairs=[],
            reads_threshold=1,
            map_order=["grouped"],
            step_name=step_name,
            output_figures_path=self.output_figures_path,
            design_file=self.design_file_path,
        )

        refiners = self._run_refiners([refiner], plot_titles=[f"{step_name} grouped"])

    def step1_reads_distribution(
        self, seq_file, bc_objects, reverse_complement, step_suffix=""
    ):
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

    def step2_reads_distribution(
        self,
        AD_seq_file,
        AD_bc_objects,
        RT_seq_file,
        RT_bc_objects,
        reverse_complement,
        step_suffix="",
    ):
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

        self.reads_distribution(
            AD_seq_file, AD_bc_objects, step_name, reverse_complement
        )
        self.reads_distribution(
            RT_seq_file, RT_bc_objects, step_name, reverse_complement
        )

    def run_step_1(
        self,
        seq_file,
        bc_objects,
        column_pairs,
        reads_threshold,
        reverse_complement=False,
        step_suffix="",
        min_fraction_major_target = 0.9,  
        custom_map_order = None
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

                **Single-column example**::

                    column_pairs = [("RPTR_BC", "AD")]

                This ensures that each reporter barcode (``RPTR_BC``) maps
                predominantly to a single AD. If a reporter barcode has
                reads mapping to multiple ADs without one exceeding 90%,
                that reporter barcode is removed.

                **Multi-column example**::

                    column_pairs = [
                        (("RPTR_BC"), ("Hawk_BC", "AD_BC"))
                    ]

                This checks that each RPTR BC
                maps to a single Hawkins AD barcode and AD barcode combination
                with ≥90% of reads. Ambiguous key combinations are removed.

                Multiple constraints may be applied sequentially::

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
            min_fraction_major_target (float): Fraction threshold for target filtering.
            custom_map_order (list): List of strings with a custom intermediate order 
                if preferred.

        Returns:
            pd.DataFrame: Final designed Step 1 mapping after all refinement
            steps have been applied.
        """

        step_name = "step1" + step_suffix

        # Initial mapping
        self._run_initial_mappers(
            [
                {
                    "seq_file": seq_file,
                    "step_name": step_name,
                    "bc_objects": bc_objects,
                    "reverse_complement": reverse_complement,
                    "design_file_path": self.design_file_path,
                }
            ]
        )

        manual_ec_threshold = True
        if self.error_correction and reads_threshold == 1:
            manual_ec_threshold = False
        
        if custom_map_order:
            map_order = custom_map_order
        else:
            map_order = MAP_ORDERS["step1"][self.error_correction]
        
        # Create refiner
        refiner = map_refiner.MapRefiner(
            db_path=self.db_path,
            bc_objects=bc_objects,
            column_pairs=column_pairs,
            reads_threshold=reads_threshold,
            map_order=map_order,
            step_name=step_name,
            output_figures_path=self.output_figures_path,
            manual_ec_threshold=manual_ec_threshold,
            design_file=self.design_file_path,
            min_fraction_major_target = min_fraction_major_target
        )

        # Use _run_refiners to handle refinement and plotting
        refiners = self._run_refiners([refiner], plot_titles=["Step 1"])

        # Get the dataframe
        df = refiners[0].get_map_df("designed")
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
        self._run_initial_mappers(
            [
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
            ]
        )

        refiners = []

        for bc_objs, reads_threshold in zip(
            (AD_bc_objects, RT_bc_objects), (reads_threshold_AD, reads_threshold_RT)
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
                    design_file=(
                        self.design_file_path if bc_objs is AD_bc_objects else None
                    ),
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
        self.con.execute(
            f"""
            CREATE TABLE {step1_map_name} AS
            SELECT * FROM read_csv_auto('{step1_map_csv_path}')
        """
        )

        # Compute AD–reporter complexity and overlap
        checker = complexity.ComplexityChecker(
            db_path=self.db_path,
            step_name=step_name,
            step1_map_name=step1_map_name,  # Use the temporary table name
            step_suffix="designed",
            barcode_groups=[AD_bc_objects, RT_bc_objects],
        )

        overlap = checker.count_overlap()

        AD_df = refiners[0].get_map_df("designed")
        RT_df = refiners[1].get_map_df("designed")

        # Save CSVs
        if self.output_path:
            AD_df.to_csv(self.output_path / f"{step_name}_AD.csv", index=False)
            RT_df.to_csv(self.output_path / f"{step_name}_RT.csv", index=False)
            overlap.to_csv(self.output_path / f"{step_name}_AD_RT_overlap.csv", index=False)

        return {"AD_step2": AD_df, "RT_step2": RT_df, "step1_overlap": overlap}

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
        name = base_name.replace(".", "_")  # replace periods with underscores
        name = base_name.replace("-", "_").replace(" ", "_")
        name = re.sub(r"[^0-9a-zA-Z_]", "_", name)
        if re.match(r"^\d", name):
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
        umi_deduplication="both",
    ):
        """
        Memory-efficient TREBL experiment runner. Writes everything to disk in the
        original output paths, no DataFrames are returned.
        """
        step_name_prefix = "trebl_experiment_" + step_name_suffix
    
        # Prepare per-library final output CSVs
        merged_output_file = self.output_path / f"{bc_objects[0].name}_trebl_experiment_{step_name_suffix}results.csv"
        if merged_output_file.exists():
            merged_output_file.unlink()
    
        for file_path in seq_files:
            try:
                base_name = os.path.basename(file_path)
                name_only = base_name.split(".")[0]
                name_only = self._duckdb_safe_name(name_only)
                step_name = f"{step_name_prefix}{name_only}"
        
                # Each sample gets its own folder
                output_dir = self.output_path / step_name_prefix / step_name
                output_dir.mkdir(parents=True, exist_ok=True)
        
                design_file_path = (
                    self.design_file_path
                    if "AD" in [bc.name for bc in bc_objects]
                    else None
                )
        
                # Initial mapping
                mapper_kwargs = dict(
                    db_path=self.db_path,
                    step_name=step_name,
                    seq_file=file_path,
                    bc_objects=bc_objects,
                    reverse_complement=reverse_complement,
                    design_file_path=design_file_path,
                )
                if umi_object:
                    mapper_kwargs["umi_object"] = umi_object
        
                mapper = initial_map.InitialMapper(**mapper_kwargs)
                mapper.create_map()
        
                # Refinement
                manual_ec_threshold = not (self.error_correction and reads_threshold == 1)
                map_order = ["quality", "error_corrected"] if self.error_correction else ["quality"]
                if not umi_object:
                    map_order += ["grouped", "thresholded"]
        
                refiner = map_refiner.MapRefiner(
                    db_path=self.db_path,
                    bc_objects=bc_objects,
                    column_pairs=[],
                    reads_threshold=reads_threshold,
                    map_order=map_order,
                    step_name=step_name,
                    output_figures_path=output_dir,
                    manual_ec_threshold=manual_ec_threshold,
                    umi_object=umi_object,
                )
                refiner.refine_map_from_db()
                refiner.plot_loss()
                plt.title(f"{name_only} TREBL Experiment")
                plt.show()
                if self.error_correction:
                    refiner.plot_error_correction()
        
                if umi_object:
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
        
                    if umi_deduplication == "simple":
                        deduplicator.run_simple_deduplication()
                        deduplicator.save_simple_deduplication()
                    else:
                        deduplicator.run_both_deduplications()
        
                    # # Merge results immediately to final CSV
                    # for suffix in ["simple", "directional_umi_counts"]:
                    #     file_name = f"{name_only}_{suffix}.tsv"
                    #     file_path_full = output_dir / file_name
                    #     if file_path_full.exists():
                    #         df = pd.read_csv(file_path_full, sep="\t")
                    #         df["name"] = name_only
                    #         if suffix == "simple":
                    #             if concat_gene:
                    #                 concat_cols = [bc.name for bc in bc_objects]
                    #                 df[gene_col_name] = df[concat_cols].agg("".join, axis=1)
                    #             else:
                    #                 bc_col = bc_objects[0].name
                    #                 df[gene_col_name] = df.get(bc_col, df.get(gene_col_name))
        
                    #         # Append to final CSV
                    #         df.to_csv(
                    #             merged_output_file,
                    #             sep=",",
                    #             index=False,
                    #             mode="a",
                    #             header=not merged_output_file.exists(),
                    #         )
        
                    # Writes reads per UMI
                    one_file_reads_per_UMI = deduplicator.counts_per_umi()
                    one_file_reads_per_UMI["name"] = name_only
                    one_file_reads_per_UMI.to_csv(
                        output_dir / f"{name_only}_reads_per_umi.tsv", sep="\t", index=False
                    )
        
                else:
                    # Non-UMI workflow: write thresholded table per sample
                    tables = refiner.show_tables()
                    first_bc_name = bc_objects[0].name
                    for table in tables:
                        if step_name_prefix in table[0] and first_bc_name in table[0] and "thresholded" in table[0]:
                            df = refiner.get_map_df(table[0])
                            bc_col_names = [bc.name for bc in bc_objects]
                            keep_cols = [col for col in df.columns if col in bc_col_names or col == "count"]
                            df = df[keep_cols].copy()
                            df["sample"] = table[0][len(step_name_prefix):]
        
                            sample_file = output_dir / f"{name_only}_thresholded.tsv"
                            df.to_csv(sample_file, sep="\t", index=False)

            except Exception as e:
                print(f"\nERROR processing {file_path}")
                print(e)
                continue
    
    
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
        umi_deduplication="both",
    ):
        """
        Memory-efficient TREBL experiment analysis for AD and RT libraries.
        Saves all outputs to the same locations as the original pipeline.
        """
        experiments = {
            "AD": {
                "seq_files": AD_seq_files,
                "bc_objects": AD_bc_objects,
                "umi_object": AD_umi_object,
                "reads_threshold": reads_threshold_AD,
                "count_col_name": "AD_umi_count",
                "gene_col_name": "AD_ADBC_concat",
                "concat_gene": True,
                "umi_deduplication": umi_deduplication,
            },
            "RT": {
                "seq_files": RT_seq_files,
                "bc_objects": RT_bc_objects,
                "umi_object": RT_umi_object,
                "reads_threshold": reads_threshold_RT,
                "count_col_name": "RTBC_umi_count",
                "gene_col_name": "RTBC",
                "concat_gene": False,
                "umi_deduplication": umi_deduplication,
            },
        }
    
        for name, spec in experiments.items():
            self._run_trebl_experiment_helper(
                seq_files=spec["seq_files"],
                bc_objects=spec["bc_objects"],
                reverse_complement=reverse_complement,
                reads_threshold=spec["reads_threshold"],
                umi_object=spec["umi_object"],
                count_col_name=spec.get("count_col_name"),
                gene_col_name=spec.get("gene_col_name"),
                concat_gene=spec.get("concat_gene", False),
                step_name_suffix=step_name_suffix,
                umi_deduplication=spec["umi_deduplication"],
            )
    
        if step1_map_csv_path:
            step_name_prefix = "trebl_experiment_" + step_name_suffix
            self.plot_trebl_experiment_loss(
                AD_bc_objects, step1_map_csv_path, step_name_prefix=step_name_prefix
            )
            self.plot_trebl_experiment_loss(
                RT_bc_objects, step1_map_csv_path, step_name_prefix=step_name_prefix
            )

    def new_plot_trebl_experiment_loss(self, 
                                        bc_objects, 
                                        step1_map_csv_path, 
                                        step_name_prefix,
                                       table_type="AD",
                                        rep_regex="AD_(\d+)",
                                        time_regex="AD_\d+_(\d+)"):

        bc_names_str = "_".join([bc.name for bc in bc_objects])

        plotting.plot_trebl_experiment_loss(
            prefix = step_name_prefix,
            bc_objects = bc_objects,
            con = duckdb.connect(self.db_path),
            output_figures_path =  self.output_figures_path,
            table_prefix_with_descriptor = f"{step_name_prefix}{bc_names_str}_bc_",
            table_type = table_type,
            step1_table=step1_map_csv_path,
            rep_regex=rep_regex,
            time_regex=time_regex,
        )

    def plot_trebl_experiment_loss(
        self,
        bc_objects,
        step1_map_csv_path=None,
        step_name_prefix="trebl_experiment_",
    ):
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
            step1_map_csv_path (str, optional): Path to Step 1 map CSV for 
                computing overlap plots.
            step_name_prefix (str, optional): Prefix used to identify TREBL
                experiment tables in DuckDB. Defaults to ``"trebl_experiment_"``.

        Returns:
            tuple: (fig, axes)
                - fig (matplotlib.figure.Figure): Figure object containing all subplots.
                - axes (list[matplotlib.axes._subplots.AxesSubplot]): Flattened list of subplot axes.
                  Returns None, None if no matching tables are found.

        """

        # Connect to DuckDB
        con = duckdb.connect(self.db_path)

        step1_map_name = None
        if step1_map_csv_path is not None:
            step1_map_name = "step1_map_temp"
            con.execute(f"DROP TABLE IF EXISTS {step1_map_name}")
            con.execute(
                f"""
                CREATE TEMP TABLE {step1_map_name} AS
                SELECT * FROM read_csv_auto('{step1_map_csv_path}')
                """
            )

        # Get all tables matching step_name_prefix and bc_object names
        tables = con.execute("SHOW TABLES").fetchall()
        bc_names = [bc.name for bc in bc_objects]

        result_prefixes = [
            table[0]
            for table in tables
            if step_name_prefix in table[0]
            and any(bc_name in table[0] for bc_name in bc_names)
            and "initial" in table[0]
        ]

        if not result_prefixes:
            print("No matching tables found.")
            return None, None

        num_plots = len(result_prefixes)
        max_cols = 5
        ncols = min(max_cols, num_plots)
        nrows = math.ceil(num_plots / ncols)

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(ncols * 4, nrows * 3),
            dpi=300,
            sharey=True,
            sharex=True,
        )
        axes = axes.flatten() if num_plots > 1 else [axes]

        bc_names_str = "_".join([bc.name for bc in bc_objects])

        for i, file_name in tqdm(
            enumerate(result_prefixes), total=len(result_prefixes), desc="Plotting BCs"
        ):
            # Total count
            total_count = con.execute(f'SELECT COUNT(*) FROM "{file_name}"').fetchone()[
                0
            ]

            # Count rows where all BCs are qualified (if the _qual column exists)
            qual_cols = [f'"{bc.name}_qual"' for bc in bc_objects]

            if qual_cols:
                qual_conditions = " AND ".join(qual_cols)
                both_true = con.execute(
                    f'SELECT COUNT(*) FROM "{file_name}" WHERE {qual_conditions}'
                ).fetchone()[0]
            else:
                both_true = 0

            # Count rows present in Step1 map for all BCs
            join_conditions = " AND ".join(
                [f'm."{bc.name}" = s."{bc.name}"' for bc in bc_objects]
            )

            if step1_map_name == None:
                step1_count = 0
            else:
                step1_count = con.execute(
                    f"""
                    SELECT COUNT(*) 
                    FROM "{file_name}" AS m
                    JOIN "{step1_map_name}" AS s
                    ON {join_conditions}
                """
                ).fetchone()[0]

            # Prepare DataFrame for plotting
            plot_counts = pd.DataFrame(
                {
                    "Category": ["Total", "BC\nQual", "In\nStep1"],
                    "Count": [total_count, both_true, step1_count],
                }
            )

            ax = axes[i]
            sns.barplot(
                x="Category",
                y="Count",
                data=plot_counts,
                ax=ax,
                palette=["gray", "green", "blue"],
            )
            for container in ax.containers:
                ax.bar_label(
                    container, fmt="%d", label_type="edge", fontsize="small", padding=2
                )

            import textwrap
            
            title_regex = f"{step_name_prefix}(.*)_{bc_names_str}_initial"
            match = re.search(title_regex, file_name)
            group_name = match.group(1) if match else file_name
            
            wrapped_title = "\n".join(textwrap.wrap(group_name, width=20))
            
            ax.set_title(wrapped_title, fontsize="medium", y=1.01)
            ax.set_xlabel("")
            ax.set_ylabel("")

        # Hide unused axes
        for j in range(num_plots, len(axes)):
            axes[j].set_visible(False)

        fig.supylabel("Count")
        sns.despine()
        plt.tight_layout(pad=1)

        fig.suptitle("Trebl Experiment Loss", y = 1.02)

        save_path = os.path.join(
            self.output_figures_path,
            f"{step_name_prefix}{bc_names_str}_bc_quality_loss.png",
        )
        plt.savefig(save_path)

        con.close()

        return fig, axes

def calculate_activity_scores(
    self,
    step1_path,
    AD_bc_objects,
    RT_bc_objects,
    time_regex=r"AD_\d+_(\d+)",
    rep_regex=r"AD_(\d+)_",
):
    """
    Calculate per-barcode activity scores from TREBL experiment AD simple and
    directional counts.

    Activity is defined as:
        np.log10(count_directional / count_simple)

    Returns
    -------
    ad_activity_per_barcode_df : pd.DataFrame
        Per-barcode, per-sample merged AD simple/directional counts with
        activity scores.
    """
    from pathlib import Path
    import glob
    import os
    import numpy as np
    import pandas as pd

    def read_and_label(files, label):
        dfs = []
        for f in files:
            df = pd.read_csv(f, sep="\t")
            df["file_base"] = f"{label}_{os.path.splitext(os.path.basename(f))[0]}"
            dfs.append(df)
        if not dfs:
            return pd.DataFrame()
        return pd.concat(dfs, ignore_index=True)

    def compute_activity(simple_vals, directional_vals):
        return np.log10(directional_vals / simple_vals)

    if self.output_path is None:
        raise ValueError("self.output_path must be set to calculate activity scores.")

    ad_col_names = [bc.name for bc in AD_bc_objects]
    rt_col_names = [bc.name for bc in RT_bc_objects]

    if len(ad_col_names) < 2:
        raise ValueError("AD_bc_objects must contain at least two barcode objects.")

    step1_required_cols = ad_col_names + rt_col_names
    step1_map = pd.read_csv(step1_path)

    missing_step1_cols = [col for col in step1_required_cols if col not in step1_map.columns]
    if missing_step1_cols:
        raise ValueError(
            f"Step 1 map is missing required columns: {missing_step1_cols}"
        )

    step1_map = step1_map[step1_required_cols].drop_duplicates()

    ad_gene_cols = ad_col_names
    step1_map["gene"] = step1_map[ad_gene_cols].astype(str).agg("".join, axis=1)

    trebl_root = Path(self.output_path) / "trebl_experiment_"
    search_root = trebl_root if trebl_root.exists() else Path(self.output_path)

    simple_AD_files = glob.glob(
        str(search_root / "**" / "*AD*" / "*simple*.tsv"),
        recursive=True,
    )
    directional_AD_files = glob.glob(
        str(search_root / "**" / "*AD*" / "*directional*.tsv"),
        recursive=True,
    )

    if len(simple_AD_files) == 0:
        raise ValueError("No AD simple count files found.")
    if len(directional_AD_files) == 0:
        raise ValueError("No AD directional count files found.")

    simple_AD_df = read_and_label(simple_AD_files, label="AD_simple")
    directional_AD_df = read_and_label(directional_AD_files, label="AD_directional")

    missing_simple_cols = [col for col in ad_gene_cols if col not in simple_AD_df.columns]
    if missing_simple_cols:
        raise ValueError(
            f"AD simple files are missing required columns: {missing_simple_cols}"
        )

    simple_AD_df["gene"] = simple_AD_df[ad_gene_cols].astype(str).agg("".join, axis=1)
    simple_AD_df["rep"] = simple_AD_df["file_base"].str.extract(rep_regex).astype(int)
    simple_AD_df["time"] = simple_AD_df["file_base"].str.extract(time_regex).astype(int)

    if "gene" not in directional_AD_df.columns:
        missing_directional_cols = [col for col in ad_gene_cols if col not in directional_AD_df.columns]
        if missing_directional_cols:
            raise ValueError(
                "AD directional files must contain either a 'gene' column "
                f"or these barcode columns: {missing_directional_cols}"
            )
        directional_AD_df["gene"] = (
            directional_AD_df[ad_gene_cols].astype(str).agg("".join, axis=1)
        )

    directional_AD_df["rep"] = directional_AD_df["file_base"].str.extract(rep_regex).astype(int)
    directional_AD_df["time"] = directional_AD_df["file_base"].str.extract(time_regex).astype(int)

    ad_activity_per_barcode_df = pd.merge(
        simple_AD_df,
        directional_AD_df,
        on=["gene", "rep", "time"],
        suffixes=("_simple", "_directional"),
    )

    step1_merge_cols = ["gene"] + step1_required_cols
    ad_activity_per_barcode_df = pd.merge(
        ad_activity_per_barcode_df,
        step1_map[step1_merge_cols].drop_duplicates(),
        on="gene",
        how="left",
    )

    ad_activity_per_barcode_df["count_diff"] = (
        ad_activity_per_barcode_df["count_simple"]
        - ad_activity_per_barcode_df["count_directional"]
    )

    ad_activity_per_barcode_df["activity_score"] = compute_activity(
        ad_activity_per_barcode_df["count_simple"],
        ad_activity_per_barcode_df["count_directional"],
    )

    output_dir = Path(self.output_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    ad_activity_per_barcode_df.to_csv(
        output_dir / "AD_activity_scores_per_barcode.csv",
        index=False,
    )

    return ad_activity_per_barcode_df
