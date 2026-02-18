from pathlib import Path
import duckdb
from trebl_tools.preprocess import time_it   

class InitialMapper:
    """Extract and map barcodes from DNA sequences in FASTQ or TXT files.

    This class reads sequence files, optionally reverse-complements them, 
    extracts barcodes and UMIs based on specified patterns, and merges with 
    a design file to produce an initial mapping table stored in DuckDB.

    Args:
        db_path (str): Path to DuckDB database file where tables will be created.
        seq_file (str or list[str]): Path(s) to input FASTQ or sequence files.
            Supports both single files and lists of multiple files.
        bc_objects (list): List of barcode objects, each with attributes:
            - name (str): Column name for the barcode
            - preceder (str): Sequence pattern before the barcode
            - post (str): Sequence pattern after the barcode  
            - length (int): Expected barcode length
        step_name (str): Identifier used as prefix for database table names.
        reverse_complement (bool): Whether to reverse complement sequences before mapping.
        design_file_path (str, optional): Path to CSV file containing designed sequences
            for validation. Expected to have sequences in first column.
        umi_object (object, optional): UMI extraction object with same attributes
            as barcode objects. Defaults to None.

    Example:
        >>> # Define barcode objects
        >>> ad_bc = BarcodeObject(name="AD_BC", preceder="ATCG", post="GCTA", length=20)
        >>> rep_bc = BarcodeObject(name="REP_BC", preceder="TTAA", post="GGCC", length=15)
        >>> 
        >>> # Create mapper
        >>> mapper = InitialMapper(
        ...     db_path="experiment.duckdb",
        ...     seq_file=["sample1.fastq.gz", "sample2.fastq.gz"],
        ...     bc_objects=[ad_bc, rep_bc],
        ...     step_name="step1",
        ...     reverse_complement=True,
        ...     design_file_path="designed_sequences.csv"
        ... )
        >>> 
        >>> # Run full mapping pipeline
        >>> mapper.create_map()
        >>> 
        >>> # Preview results
        >>> preview_df = mapper.preview_map()
        >>> print(preview_df.head())

    Note:
        For FASTQ files, only sequence lines (every 4th line starting from line 2)
        are extracted. Quality flags are automatically generated for each extracted
        barcode based on expected length matching.
    """

    def __init__(self, 
                 db_path, 
                 seq_file, 
                 bc_objects, 
                 step_name, 
                 reverse_complement, 
                 design_file_path, 
                 umi_object=None):
        """Initialize the InitialMapper object and set up DuckDB connection.
        
        Args:
            db_path (str): Path to DuckDB database file.
            seq_file (str or list[str]): Input sequence file path(s).
            bc_objects (list): List of barcode extraction objects.
            step_name (str): Step identifier for table naming.
            reverse_complement (bool): Whether to reverse complement sequences.
            design_file_path (str): Path to design validation file.
            umi_object (object, optional): UMI extraction object.
        """        
        self.con = duckdb.connect(db_path)
        self.seq_files = seq_file if isinstance(seq_file, list) else [seq_file]
        self.design_file_path = str(Path(design_file_path).resolve()) if design_file_path else None
        self.bc_objects = bc_objects
        self.reverse_complement = reverse_complement
        self.umi_object = umi_object
        self.step_name = step_name

        self.cols = [bc_object.name for bc_object in bc_objects]
        cols_str = "_".join(self.cols)
        self.table_prefix = f"{step_name}_{cols_str}_"
        self.test_n_reads = None

    @time_it
    def read_fastq(self):
        """Read sequences from FASTQ/TXT files and create an initial sequence table.

        Processes input files to extract sequence data. For FASTQ files, extracts
        only the sequence lines (every 4th line starting from line 2). Creates
        a 'seq' table in the connected DuckDB database with all sequences.

        Note:
            If `self.test_n_reads` is set, limits the number of reads processed
            for testing purposes. 

        Example:
            >>> mapper.test_n_reads = 1000  # Optional: limit for testing
            >>> mapper.read_fastq()
            Reading 2 FASTQ/TXT file(s)...
            Done in 15.34 seconds.
        """    
        con = self.con
        if isinstance(self.seq_files, str):
            self.seq_files = [self.seq_files]

        print(f"Reading {len(self.seq_files)} FASTQ/TXT file(s)...")
        seq_queries = []
        for f in self.seq_files:
            f_path = str(Path(f).resolve())
            seq_queries.append(f"""
                SELECT column0 AS sequence
                FROM (
                    SELECT *, ROW_NUMBER() OVER () AS rn
                    FROM read_csv_auto('{f_path}', header=False)
                ) t
                WHERE rn % 4 = 2
            """)

        combined_query = " UNION ALL ".join(seq_queries)

        if self.test_n_reads:
            # Limit to test_n_reads
            combined_query = f"SELECT * FROM ({combined_query}) LIMIT {self.test_n_reads}"

        con.execute(f"CREATE OR REPLACE TABLE seq AS {combined_query};")

        
    @time_it
    def reverse_complement_fastq(self):
        """Replace sequence column with reverse complement of DNA sequences.
        
        Converts all sequences in the 'seq' table to their reverse complement
        using DuckDB's built-in string functions. Only processes sequences if
        self.reverse_complement is True.
            
        Example:
            >>> mapper.reverse_complement = True
            >>> mapper.reverse_complement_fastq()
            Reverse complement of sequences...
            Done in 8.21 seconds.
        """
        if not self.reverse_complement:
            return
        print("Reverse complement of sequences...")
        con = self.con
    
        # Create a new table with only the reverse complemented sequence
        con.execute("""
            CREATE TABLE seq_rc AS
            SELECT reverse(translate(sequence, 'ACGTacgt', 'TGCAtgca')) AS sequence
            FROM seq;
        """)
    
        # Drop the old table and rename
        con.execute("DROP TABLE seq;")
        con.execute("ALTER TABLE seq_rc RENAME TO seq;")

    def _extract_sequence_object(self, seq_obj):
        """Extract a barcode or UMI sequence into a new column and compute quality.

        Internal helper method that extracts sequences based on the provided
        sequence object's pattern specifications. Handles four extraction modes:
        1. Both preceder and post: Extract sequence between patterns
        2. Only preceder: Extract N bases after preceder pattern
        3. Only post: Extract N bases before post pattern  
        4. Neither: Extract last N bases from sequence

        Args:
            seq_obj (object): Barcode or UMI object with attributes:
                - name (str): Column name for extracted sequence
                - preceder (str or None): Pattern before target sequence
                - post (str or None): Pattern after target sequence
                - length (int): Expected length of target sequence

        Note:
            Creates two columns: {name} for the extracted sequence and
            {name}_qual as a boolean flag indicating if extraction succeeded
            (length matches expected). Uses non-greedy regex matching when
            both preceder and post patterns are specified.

        Example:
            >>> # Extract 20bp barcode between "ATCG" and "GCTA"
            >>> bc_obj = BarcodeObject(name="BC1", preceder="ATCG", post="GCTA", length=20)
            >>> mapper._extract_sequence_object(bc_obj)
            BC1: extracting between 'ATCG' and 'GCTA'
        """
        con = self.con

        col = seq_obj.name
        length = seq_obj.length
        con.execute(f"ALTER TABLE seq ADD COLUMN {col} VARCHAR;")

        if seq_obj.preceder and seq_obj.post:
            # Both exist → regex between preceder and post
            regex = f"{seq_obj.preceder}(.*?){seq_obj.post}"  # non-greedy
            print(f"{col}: extracting between '{seq_obj.preceder}' and '{seq_obj.post}'")
            con.execute(f"""
                UPDATE seq
                SET {col} = coalesce(regexp_extract(sequence::VARCHAR, '{regex}', 1), '')
            """)
    
        elif not seq_obj.preceder and not seq_obj.post:
            # Neither exists → take last N bases
            print(f"{col}: extracting last {length} bases")
            con.execute(f"""
                UPDATE seq
                SET {col} = CASE
                    WHEN length(sequence) >= {length}
                    THEN substr(sequence, length(sequence) - {length} + 1, {length})
                    ELSE ''
                END
            """)
    
        elif seq_obj.preceder:
            # Only preceder → take N bases after preceder
            regex = f"{seq_obj.preceder}(.{{{length}}})"
            print(f"{col}: extracting {length} bases after preceder '{seq_obj.preceder}'")
            con.execute(f"""
                UPDATE seq
                SET {col} = coalesce(regexp_extract(sequence::VARCHAR, '{regex}', 1), '')
            """)
    
        elif seq_obj.post:
            # Only post → take N bases before post
            regex = f"(.{{{length}}}){seq_obj.post}"
            print(f"{col}: extracting {length} bases before post '{seq_obj.post}'")
            con.execute(f"""
                UPDATE seq
                SET {col} = coalesce(regexp_extract(sequence::VARCHAR, '{regex}', 1), '')
            """)
    
        # Add quality flag column
        qual_col = f"{col}_qual"
        con.execute(f"ALTER TABLE seq ADD COLUMN {qual_col} BOOLEAN;")
        con.execute(f"""
            UPDATE seq
            SET {qual_col} = LENGTH({col}) = {length};
        """)
        
    @time_it
    def extract_barcodes(self):
        """Extract all barcode sequences and compute their quality flags.
        
        Processes each barcode object in self.bc_objects to extract the
        corresponding barcode sequences from the 'seq' table. Creates
        both the barcode columns and their associated quality flag columns.
        
        Note:
            Calls _extract_sequence_object() for each barcode. Quality flags
            are boolean columns indicating whether the extracted sequence
            matches the expected length. 
            
        Example:
            >>> mapper.extract_barcodes()
            Extracting 3 barcodes...
            AD_BC: extracting between 'ATCG' and 'GCTA'
            RT_BC: extracting 15 bases after preceder 'TTAA'
            UMI: extracting last 12 bases
            Done in 45.67 seconds.
        """
        print(f"Extracting {len(self.bc_objects)} barcodes...")
        for bc in self.bc_objects:
            self._extract_sequence_object(bc)

    @time_it
    def extract_umi(self):
        """Extract the UMI sequence and compute its quality flag.
        
        Processes the UMI object (if provided) to extract UMI sequences
        from the 'seq' table using the same pattern-matching approach
        as barcode extraction.
        
        Note:
            Only runs if self.umi_object is not None. Creates both UMI
            column and UMI quality flag column. 
            
        Example:
            >>> mapper.extract_umi()
            Extracting UMI...
            UMI: extracting 12 bases after preceder 'AAAA'
            Done in 12.34 seconds.
        """
        if self.umi_object:
            print("Extracting UMI...")
            self._extract_sequence_object(self.umi_object)
                
    @time_it
    def merge_design(self):        
        """Merge with design file and create final mapping table.
        
        Merges extracted barcode data with the design file (if provided) to
        create a 'Designed' column indicating whether sequences match the
        designed library. Creates the final initial mapping table and removes
        the raw sequence column to save space.
        
        Note:
            If design_file_path is provided and "AD" column exists, performs
            LEFT JOIN with design file. Otherwise, sets all sequences as designed (1).
            Creates final table with name pattern: {table_prefix}initial
            
        Example:
            >>> mapper.merge_design()
            Merging with design file...
            Done in 5.43 seconds.
            
        Raises:
            FileNotFoundError: If design_file_path is specified but file doesn't exist.
        """  
        con = self.con
        if self.design_file_path and "AD" in self.cols:
            print("Merging with design file...")
            con.execute(f"""
                CREATE OR REPLACE TABLE design AS
                SELECT CAST(column0 AS VARCHAR) AS AD
                FROM read_csv_auto('{self.design_file_path}', header=False)
            """)
            con.execute(f"""
                CREATE OR REPLACE TABLE {self.table_prefix}initial_tmp AS
                SELECT s.*, CASE WHEN d.AD IS NOT NULL THEN 1 ELSE 0 END AS Designed
                FROM seq s
                LEFT JOIN design d USING(AD);
            """)
        else:
            con.execute(f"""
                CREATE OR REPLACE TABLE {self.table_prefix}initial_tmp AS
                SELECT *, 1 AS Designed
                FROM seq;
            """)
        con.execute("DROP TABLE IF EXISTS seq;")

        # Drop the sequence column from the new table
        con.execute(f"ALTER TABLE {self.table_prefix}initial_tmp DROP COLUMN sequence;")
        con.execute(
            f"ALTER TABLE {self.table_prefix}initial_tmp RENAME TO {self.table_prefix}initial;"
        )
        
    def _run_pipeline(self):
        """Execute the complete mapping pipeline in sequence.
        
        Internal helper method that runs all pipeline steps in the correct order:
        1. read_fastq() - Load sequences from files
        2. reverse_complement_fastq() - Optionally reverse complement
        3. extract_barcodes() - Extract all barcode sequences  
        4. extract_umi() - Extract UMI sequences (if specified)
        5. merge_design() - Merge with design file and finalize table
        
        Note:
            This method is called by both create_map() and create_test_map().
            All individual steps are decorated with @time_it for performance monitoring.
        """
        self.read_fastq()
        self.reverse_complement_fastq()
        self.extract_barcodes()
        self.extract_umi()
        self.merge_design()
        print("Mapping complete.")

    def _table_exists(self, table_name: str) -> bool:
        """Check if a table exists in the connected DuckDB database.
        
        Args:
            table_name (str): Name of the table to check for existence.
            
        Returns:
            bool: True if table exists, False otherwise.
            
        Example:
            >>> exists = mapper._table_exists("step1_AD_BC_initial")
            >>> print(f"Table exists: {exists}")
        """
        return self.con.execute(
            """
            SELECT COUNT(*) > 0
            FROM information_schema.tables
            WHERE table_name = ?
            """,
            [table_name],
        ).fetchone()[0]

    
    def create_map(self):
        """Run the full mapping pipeline on all reads.
        
        Executes the complete barcode extraction pipeline on all sequences
        in the input files. Skips processing if the initial mapping table
        already exists in the database.
        
        Note:
            Creates table with name pattern: {step_name}_{barcode_cols}_initial
            Resets test_n_reads to None to ensure all reads are processed.
            All pipeline steps are automatically timed and logged.
            
        Example:
            >>> mapper.create_map()
            ✓ Initial map already exists: step1_AD_BC_REP_BC_initial — skipping
            
            # Or if table doesn't exist:
            Reading 1 FASTQ/TXT file(s)...
            Done in 30.45 seconds.
            ...
            Mapping complete.
        """   
        initial_table = f"{self.table_prefix}initial"
        
        if self._table_exists(initial_table):
            print(f"✓ Initial map already exists: {initial_table} — skipping")
            return
    
        self.test_n_reads = None  # ensure no test limit
        self._run_pipeline()

    
    def create_test_map(self, test_n_reads: int = 100):
        """Run the mapping pipeline on a limited number of reads for testing.
        
        Executes the barcode extraction pipeline on a subset of sequences
        for rapid testing and validation. Skips if initial mapping table
        already exists.
        
        Args:
            test_n_reads (int, optional): Maximum number of reads to process
                for testing. Defaults to 100.
                
        Note:
            Useful for validating extraction patterns and pipeline functionality
            before processing large datasets. Resets test limit after completion.
            
        Example:
            >>> mapper.create_test_map(test_n_reads=500)
            Reading 1 FASTQ/TXT file(s)...
            Done in 2.15 seconds.
            ...
            Mapping complete.
        """ 
        initial_table = f"{self.table_prefix}initial"

        if self._table_exists(initial_table):
            print(f"✓ Initial map already exists: {initial_table} — skipping test map")
            return
    
        self.test_n_reads = test_n_reads
        self._run_pipeline()
        self.test_n_reads = None
  

    def preview_map(self):
        """Preview the first 5 rows of the created mapping table.

        Returns a sample of the initial mapping table to verify extraction
        results and inspect data quality before downstream processing.

        Returns:
            pd.DataFrame: DataFrame containing the first 5 rows of the initial
                mapping table with all extracted barcodes, UMI (if applicable),
                quality flags, and design validation results.

        Example:
            >>> preview_df = mapper.preview_map()
            step1_AD_BC_REP_BC_initial
            Total rows: 1500000
            >>> print(preview_df)
               AD_BC  AD_BC_qual  REP_BC  REP_BC_qual  UMI  UMI_qual  Designed
            0  ATCG...    True      TTAA...     True      ...    True         1
            1  GCTA...    False     GGCC...     True      ...    True         0
            ...

        Note:
            Also prints the table name and total row count for quick reference.
            Requires that create_map() or create_test_map() has been run first.
        """
        print(f"{self.table_prefix}initial")

        # Count total rows
        total_rows = self.con.execute(f"SELECT COUNT(*) FROM {self.table_prefix}initial").fetchone()[0]
        print(f"Total rows: {total_rows}")
        
        return self.con.execute(f"SELECT * FROM {self.table_prefix}initial LIMIT 5;").df()
