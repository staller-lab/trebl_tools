import matplotlib
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd

def plot_reads_histogram(map_df, save_path=None, ax=None, **kwargs):    
    """Plot a histogram of read counts on a log-log scale.

    Creates a histogram visualization of read count distributions using
    logarithmic scaling on both axes for better visualization of wide
    dynamic ranges typical in sequencing data.

    Args:
        map_df (pd.DataFrame): DataFrame containing at least a "count" column
            with read count values.
        save_path (str, optional): Path to save the figure. If None, figure
            is not saved. Defaults to None.
        ax (matplotlib.axes.Axes, optional): Existing axis to plot on. If None,
            a new figure is created. Defaults to None.
        **kwargs: Additional keyword arguments passed to sns.histplot().

    Returns:
        matplotlib.axes.Axes: The axis object with the histogram plotted.

    Example:
        >>> df = pd.DataFrame({'count': [10, 100, 1000, 50, 500]})
        >>> ax = plot_reads_histogram(df, save_path="read_dist.png")
        >>> plt.show()

    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
    
    sns.histplot(map_df["count"], log_scale=(True, True), ax=ax, **kwargs)
    
    # Ensure log scale is applied
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel("Read Count")
    
    sns.despine()
    
    if save_path:
        ax.get_figure().savefig(save_path, bbox_inches="tight")
    
    return ax
    
    
def plot_error_correction(output_figures_path, table_prefix_with_descriptor, plot=True):
    """Summarize and visualize error correction results from whitelist files.
    
    Processes whitelist files generated during error correction to create
    summary statistics and visualizations of barcode collapse and correction
    performance. Reads the first whitelist file found and generates comprehensive
    plots showing before/after sequence counts and correction statistics.
    
    Args:
        output_figures_path (str or Path): Directory containing whitelist files
            and where output plots will be saved.
        table_prefix_with_descriptor (str): Prefix for output filenames.
        plot (bool, optional): Whether to generate and save plots. 
            Defaults to True.
    
    Returns:
        pd.DataFrame: Summary table with columns:
            - canonical: Canonical barcode sequence
            - collapsed: Comma-separated list of collapsed sequences
            - largest_count: Read count of most abundant sequence in group
            - counts: Comma-separated counts of collapsed sequences
            - collapsed_list: List of collapsed sequences
            - num_merged: Number of sequences merged into canonical
            - rest_count: Sum of counts from non-canonical sequences
    
    Raises:
        FileNotFoundError: If no whitelist files (*_whitelist.txt) are found
            in the specified directory.
    
    Example:
        >>> summary = plot_error_correction(
        ...     "results/error_corrected/",
        ...     "step1_",
        ...     plot=True
        ... )
        >>> print(f"Processed {len(summary)} canonical barcodes")
        Saved barcode whitelist plots: results/error_corrected/step1_whitelist_summary.png
        
    Note:
        Expects whitelist files to have tab-separated format with columns:
        canonical, collapsed, largest_count, counts. Uses seaborn 'talk' context
        for plot styling.
    """
    import glob
    
    whitelist_dir = Path(output_figures_path or "./error_corrected")
    whitelist_files = sorted(glob.glob(str(whitelist_dir / "*_whitelist.txt")))
    print(whitelist_files)
    
    if not whitelist_files:
        raise FileNotFoundError(f"No whitelist files found in {whitelist_dir}")
    
    
    wl_df = pd.read_csv(whitelist_files[0], sep="\t", header=None,
                        names=["canonical", "collapsed", "largest_count", "counts"])
    
    # Process collapsed sequences
    wl_df["collapsed_list"] = wl_df["collapsed"].fillna("").apply(lambda x: x.split(",") if x else [])
    wl_df["num_merged"] = wl_df["collapsed_list"].apply(len)
    
    # Compute rest counts
    def parse_counts(row):
        if pd.isna(row["counts"]) or str(row["counts"]).strip() == "":
            return row["largest_count"], 0
        rest = sum(int(c) for c in str(row["counts"]).split(",") if c.strip())
        largest = row["largest_count"] - rest
        return largest, rest
    
    wl_df[["largest_count", "rest_count"]] = wl_df.apply(parse_counts, axis=1, result_type="expand")
    
    summary_df = wl_df
    
    #Plot distributions
    if plot:
        sns.set(style="white", context="talk")
    
        fig = plot_all_whitelists_from_summary(summary_df)
        if output_figures_path:
            plot_path = Path(output_figures_path) / f"{table_prefix_with_descriptor}whitelist_summary.png"
            fig.savefig(plot_path, bbox_inches="tight")
            print(f"Saved barcode whitelist plots: {plot_path}")
        plt.show()
    
    return summary_df
    
def plot_all_whitelists_from_summary(summary_df, n_cols=4, dpi=300):
    """Create comprehensive visualization of barcode error correction summary.
    
    Generates a multi-panel figure showing various aspects of barcode error
    correction performance including sequence reduction, group size distributions,
    and read count relationships between canonical and collapsed sequences.
    
    Args:
        summary_df (pd.DataFrame): DataFrame with error correction summary data.
            Must contain columns: ['canonical', 'num_merged', 'largest_count', 'rest_count'].
        n_cols (int, optional): Number of panels per row. Currently fixed at 4
            for the standard panel layout. Defaults to 4.
        dpi (int, optional): Figure resolution in dots per inch. Defaults to 300.
    
    Returns:
        matplotlib.figure.Figure: Figure object containing four panels:
            - Panel 1: Bar chart of sequence counts before/after correction
            - Panel 2: Histogram of group sizes (number of sequences per canonical)
            - Panel 3: Scatter plot of largest member count vs group size  
            - Panel 4: Log-log scatter of largest vs summed smaller member counts
    
    Example:
        >>> summary_df = pd.DataFrame({
        ...     'num_merged': [0, 2, 1, 3],
        ...     'largest_count': [100, 500, 200, 1000], 
        ...     'rest_count': [0, 50, 25, 150]
        ... })
        >>> fig = plot_all_whitelists_from_summary(summary_df)
        >>> plt.show()
        
    Note:
        Panel 4 includes a diagonal reference line (y=x) to show the relationship
        between canonical and collapsed sequence read counts. Uses seaborn 'Paired'
        color palette for consistent styling.
    """
    n_barcodes = 1
    n_rows = n_barcodes
    n_panels = 4  # same as summarize_whitelist
    
    fig, axs = plt.subplots(n_rows, n_panels, figsize=(5 * n_panels, 5 * n_rows), dpi=dpi)
    
    # Ensure axs is 2D array for consistent indexing
    if n_rows == 1:
        axs = axs.reshape(1, -1)
    
    i = 0
    df_bc = summary_df
    
    # Panel 1: Bar chart of original vs canonical sequences
    total_original_seqs = df_bc["num_merged"].sum() + len(df_bc)
    total_canonical_seqs = len(df_bc)
    axs[i, 0].bar(["Before", "After"], [total_original_seqs, total_canonical_seqs],
                  color=[sns.color_palette('Paired')[0], sns.color_palette('Paired')[1]])
    axs[i, 0].set_ylabel("Number of sequences")
    for x, y in zip(["Before", "After"], [total_original_seqs, total_canonical_seqs]):
        axs[i, 0].text(x, y + max([total_original_seqs, total_canonical_seqs])*0.02, f"{y:,}", 
                       ha='center', va='bottom', fontsize='medium', weight='bold')
    
    # Panel 2: Histogram of group sizes
    axs[i, 1].hist(df_bc["num_merged"] + 1, bins=30, edgecolor='black')
    axs[i, 1].set_xlabel("Group size")
    axs[i, 1].set_ylabel("Frequency")
    
    # Panel 3: Scatter - largest member vs group size
    axs[i, 2].scatter(df_bc["largest_count"], df_bc["num_merged"] + 1, alpha=0.6)
    axs[i, 2].set_xlabel("Reads of largest")
    axs[i, 2].set_ylabel("Group size")
    
    # Panel 4: Scatter - largest member vs sum of smaller members
    axs[i, 3].scatter(df_bc["largest_count"], df_bc["rest_count"], alpha=0.6)
    axs[i, 3].set_xlabel("Reads of largest")
    axs[i, 3].set_ylabel("Summed merged reads")
    axs[i, 3].set_xscale('log')
    axs[i, 3].set_yscale('log')
    
    # Plot y=x line for panel 4
    xlims = axs[i, 3].get_xlim()
    ylims = axs[i, 3].get_ylim()
    line_min = min(xlims[0], ylims[0])
    line_max = max(xlims[1], ylims[1])
    axs[i, 3].plot([line_min, line_max], [line_min, line_max], 'r--', alpha=0.7)
    axs[i, 3].set_xlim(xlims)
    
    # Add row-level label for the barcode
    fig.text(
        -0.02,
        (n_rows - i - 0.5) / n_rows,
        f"Concat",
        fontsize='large',
        weight='bold',
        rotation=90,
        va='center'
    )
    
    sns.despine()
    plt.tight_layout(pad=2)
    return fig

def plot_loss_helper(ax, palette, text_offset, show_background, default_map_order, output_figures_path, table_prefix_with_descriptor, df):
    """Generate horizontal bar plot showing read loss through processing pipeline steps.
    
    Creates a visualization showing both total reads and unique counts at each
    step of the processing pipeline, with optional background bars to show
    the difference between total and unique counts.
    
    Args:
        ax (matplotlib.axes.Axes or None): Existing axis to plot on. If None,
            a new figure is created.
        palette (str): Name of seaborn color palette to use for bars.
        text_offset (float): Vertical offset for read count labels to avoid overlap.
        show_background (bool): Whether to show background bars for total reads
            with reduced alpha transparency.
        default_map_order (list): Ordered list of pipeline step names for
            consistent color mapping across plots.
        output_figures_path (str): Directory path where plot will be saved.
        table_prefix_with_descriptor (str): Prefix for output filename.
        df (pd.DataFrame): DataFrame with columns:
            - description: Step descriptions for y-axis labels  
            - map: Pipeline step names for color mapping
            - total_reads: Total read counts per step
            - unique_count: Unique sequence counts per step
            
    Returns:
        tuple: (fig, ax) matplotlib Figure and Axes objects.
        
    Example:
        >>> df = pd.DataFrame({
        ...     'description': ['Step 1', 'Step 2'], 
        ...     'map': ['map1', 'map2'],
        ...     'total_reads': [10000, 8000],
        ...     'unique_count': [5000, 4000]
        ... })
        >>> fig, ax = plot_loss_helper(
        ...     None, 'Set2', 0.1, True, ['map1', 'map2'], 
        ...     'results/', 'sample_', df
        ... )
        
    Note:
        Automatically saves plot to output_figures_path with filename 
        "{table_prefix_with_descriptor}loss.png". 
    """
    # Initialize plot
    if ax is None:
        sns.set(style="white", context="talk")
        fig, ax = plt.subplots(figsize=(6, 4), dpi=300)

    # Create color mapping
    cmap = sns.color_palette(palette, n_colors=len(default_map_order))
    palette_dict = {name: cmap[i] for i, name in enumerate(default_map_order)}
    colors = [palette_dict.get(name, (0.5, 0.5, 0.5)) for name in df["map"]]

    # Bar plots for total reads (light) and unique counts (solid)
    if show_background:
        background_alpha = 0.5
    else:
        background_alpha = 0
    sns.barplot(x="total_reads", y="description", data=df, ax=ax, palette=colors, alpha=background_alpha)
    sns.barplot(x="unique_count", y="description", data=df, ax=ax, palette=colors, alpha = 1)
    #sns.scatterplot(x="unique_AD_count", y="description", data=df, ax=ax, palette=colors, zorder = 20)

    # Draw black line marking unique_AD_count — only show top edge
    # sns.barplot(
    #     x="unique_AD_count", y="description", data=df,
    #     ax=ax, color='none', zorder=20
    # )
        
    # Add count labels to bars
    for i, row in df.iterrows():
        # Unique counts (front bar)
        if pd.notna(row["unique_count"]):
            if row["unique_count"] != row["total_reads"]:
                ax.text(
                    row["unique_count"] * 1.01,
                    i-text_offset,
                    f"{int(row['unique_count']):,}",
                    va="center",
                    ha="left",
                    fontsize='small',
                    color="black"#, zorder = 20
                )
        if show_background:
            background_color = 'black'
        else:
            background_color = 'white'
        # Total reads (back bar)
        if pd.notna(row["total_reads"]):
            ax.text(
                row["total_reads"] * 1.01,
                i+text_offset,
                f"{int(row['total_reads']):,}",
                va="center",
                ha="left",
                fontsize='small',
                color=background_color#, zorder = 20
            )
    
    ax.set_xlabel("Read Count (Unique, Total)")
    ax.set_ylabel("Map Step")
    ax.set_yticklabels(df["map"])
    sns.despine(bottom=True)
    ax.set_xticks([])

    if output_figures_path:
        filename = os.path.join(output_figures_path, f"{table_prefix_with_descriptor}loss.png")
        plt.savefig(filename, bbox_inches="tight")

    return fig, ax
    