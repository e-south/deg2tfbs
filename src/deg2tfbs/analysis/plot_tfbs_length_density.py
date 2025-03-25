"""
--------------------------------------------------------------------------------
<deg2tfbs project>

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def plot_and_save_tfbs_length_analysis(mapping_csv: Path, plot_output: Path, summary_output: Path) -> None:
    """
    Loads the tf2tfbs_mapping CSV file, computes the length of each TFBS string,
    generates a density plot with hue based on the 'tf' column,
    computes summary statistics per TF, and saves both the plot and a CSV summary.

    Parameters:
      mapping_csv (Path): Path to the tf2tfbs_mapping.csv file.
      plot_output (Path): Path where the density plot (tfbs_length_density.png) will be saved.
      summary_output (Path): Path where the summary CSV (tfbs_length_summary.csv) will be saved.
    """
    # Load the mapping CSV (assuming a standard CSV format)
    try:
        df = pd.read_csv(mapping_csv)
    except Exception as e:
        raise RuntimeError(f"Error loading mapping CSV from {mapping_csv}: {e}")
    
    # Ensure the expected columns are present
    required_columns = ['tf', 'tfbs']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in {mapping_csv}")

    # Compute the length of each binding site (as the length of the string in 'tfbs')
    df['tfbs_length'] = df['tfbs'].apply(lambda x: len(str(x)) if pd.notnull(x) else 0)
    
    # Set seaborn style and palette
    sns.set_style("ticks")
    sns.set_palette("pastel")
    
    # Create the density plot
    plt.figure(figsize=(8, 6))
    # Plot a density plot (KDE) for each unique TF
    # This will overlay many density estimates (one per unique TF)
    unique_tfs = df['tf'].unique()
    for tf in unique_tfs:
        subset = df[df['tf'] == tf]
        # If there is not enough data to compute a density, skip plotting for that TF.
        if len(subset) < 2:
            continue
        sns.kdeplot(subset['tfbs_length'], label=tf, fill=True)
    
    plt.xlabel("TFBS Length (bp)")
    plt.ylabel("Density")
    plt.title("TFBS Length Density per TF")
    plt.legend(title="TF", bbox_to_anchor=(1.05, 1), loc="upper left")
    sns.despine()  # Remove top and right spines
    
    # Save the density plot
    plt.tight_layout()
    plt.savefig(plot_output, dpi=150)
    plt.close()
    
    # Compute summary statistics per TF
    summary = df.groupby("tf")['tfbs_length'].agg(
        n_sites="count",
        mean_length="mean",
        std_length="std",
        min_length="min",
        max_length="max"
    ).reset_index()
    
    # Sort summary by descending count of sites
    summary = summary.sort_values(by="n_sites", ascending=False)
    
    # Save summary CSV
    summary.to_csv(summary_output, index=False)
