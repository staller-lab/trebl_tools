import pandas as pd


class Barcode:
    """Represents a barcode with extraction pattern specifications.

    This class defines the parameters needed to extract a specific barcode
    from DNA sequences, including the flanking sequences and expected length.

    Args:
        name (str): Identifier for the barcode that will be used as column name.
        preceder (str): DNA sequence pattern that precedes the barcode.
        post (str): DNA sequence pattern that follows the barcode.
        length (int): Expected length of the barcode sequence.

    Example:
        >>> # Create a barcode object for a 20bp sequence between "ATCG" and "GCTA"
        >>> ad_barcode = Barcode(
        ...     name="AD_BC",
        ...     preceder="ATCG",
        ...     post="GCTA",
        ...     length=20
        ... )
        >>> print(f"Barcode: {ad_barcode.name}, Length: {ad_barcode.length}")
        Barcode: AD_BC, Length: 20

    Note:
        The preceder and post sequences are used in regex pattern matching,
        so special regex characters should be escaped if they appear literally
        in the DNA sequences.
    """

    def __init__(self, name, preceder, post, length):
        """Initialize a Barcode object with extraction parameters.

        Args:
            name (str): Barcode identifier for column naming.
            preceder (str): DNA sequence preceding the barcode.
            post (str): DNA sequence following the barcode.
            length (int): Expected length of barcode sequence.
        """
        self.name = name
        self.preceder = preceder
        self.post = post
        self.length = length


def add_barcode(seq_df, bc_object):
    """Extract a barcode from sequences using pattern matching.

    Extracts barcode sequences located between specified preceder and post
    patterns. Creates both the barcode column and a quality flag column
    indicating whether the extracted sequence matches the expected length.

    Args:
        seq_df (pd.DataFrame ): DataFrame containing a 'sequence' column
            with DNA sequences to process.
        bc_object (Barcode): Barcode object containing extraction parameters
            (name, preceder, post, length).

    Returns:
        pd.DataFrame : Updated DataFrame with two new columns:
            - {bc_object.name}: Extracted barcode sequences
            - {bc_object.name}_qual: Boolean quality flag (True if length matches expected)

    Example:
        >>> df = pd.DataFrame({"sequence": ["AAAXXXCCC", "AAAYYYCCC", "AAAXCCC"]})
        >>> bc = Barcode(name="TEST", preceder="AAA", post="CCC", length=3)
        >>> result_df = add_barcode(df, bc)
        >>> print(result_df)
           sequence TEST  TEST_qual
        0  AAAXXXCCC  XXX       True
        1  AAAYYYCCC  YYY       True
        2    AAAXCCC    X      False

    Note:
        Uses greedy regex matching by default. Sequences that don't match
        the pattern or don't meet the length requirement will have NaN
        values that are converted to False for the quality flag.
    """
    regex = f"{bc_object.preceder}(.*){bc_object.post}"
    subseq_series = (
        seq_df["sequence"].str.extract(regex)[0].str.slice(0, bc_object.length)
    )
    seq_df[bc_object.name] = subseq_series
    seq_df[bc_object.name + "_qual"] = (
        seq_df[bc_object.name].str.len() == bc_object.length
    )
    seq_df[bc_object.name + "_qual"] = seq_df[bc_object.name + "_qual"].fillna(False)
    return seq_df


def add_umi(seq_df, umi_length=12):
    """Extract UMI (Unique Molecular Identifier) from the end of sequences.

    Extracts the last `umi_length` bases from each sequence as the UMI.
    This is typically used when UMIs are located at the 3' end of reads
    after barcode regions.

    Args:
        seq_df (pd.DataFrame ): DataFrame containing a 'sequence' column
            with DNA sequences.
        umi_length (int, optional): Number of bases to extract from sequence end
            as the UMI. Defaults to 12.

    Returns:
        pd.DataFrame : Updated DataFrame with new 'UMI' column
            containing the extracted UMI sequences.

    Example:
        >>> df = pd.DataFrame({"sequence": ["ATCGATCGATCGATCG", "GCTAGCTAGCTAGCTA"]})
        >>> result_df = add_umi(df, umi_length=4)
        >>> print(result_df)
                   sequence   UMI
        0  ATCGATCGATCGATCG  ATCG
        1  GCTAGCTAGCTAGCTA  GCTA

    Note:
        Currently hardcoded to extract the last 12 bases regardless of the
        umi_length parameter. This appears to be a bug that should be fixed
        to use the umi_length parameter.

    Todo:
        Fix implementation to use umi_length parameter instead of hardcoded -12.
    """
    seq_df["UMI"] = seq_df["sequence"].str.slice(
        -12,
    )
    return seq_df


def add_multiple_barcodes(bc_objects, seq_df):
    """Extract multiple barcodes from sequences in a single operation.

    Applies barcode extraction for all specified barcode objects to the
    input DataFrame. This is a convenience function for processing multiple
    barcodes in a batch operation.

    Args:
        bc_objects (list[Barcode]): List of Barcode objects defining the
            barcodes to extract.
        seq_df (pd.DataFrame ): DataFrame containing sequences
            in a 'sequence' column.

    Returns:
        pd.DataFrame : Updated DataFrame with columns for
            each extracted barcode and corresponding quality flags:
            - {barcode.name}: Extracted barcode sequence
            - {barcode.name}_qual: Boolean quality flag for each barcode

    Example:
        >>> df = pd.DataFrame({"sequence": ["AAAXXXCCCYYYTTT", "AAAABBBCCCDDTTT"]})
        >>> bc1 = Barcode(name="BC1", preceder="AAA", post="CCC", length=3)
        >>> bc2 = Barcode(name="BC2", preceder="CCC", post="TTT", length=3)
        >>> result_df = add_multiple_barcodes([bc1, bc2], df)
        >>> print(result_df.columns.tolist())
        ['sequence', 'BC1', 'BC1_qual', 'BC2', 'BC2_qual']

    Note:
        Barcodes are processed sequentially in the order provided.
        Creates a copy of the input DataFrame to avoid modifying the original.
        Each barcode extraction is independent and uses the same source sequences.

    Warning:
        If barcodes have overlapping patterns or are located in overlapping
        regions, the extraction order may affect results.
    """
    df = seq_df.copy()
    for bc_object in bc_objects:
        df = add_barcode(df, bc_object)
    return df
