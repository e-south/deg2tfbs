"""
--------------------------------------------------------------------------------
<deg2tfbs project>
src/deg2tfbs/analysis/plot_tfbs_counts.py

For a given TFBS source CSV, generate a bar plot of the frequency counts of each
transcription factor. The operation groups by the 'tf' column and counts the number
of unique 'tfbs' values. The plot title is set to the name of the parent directory 
of source_csv. 

If an optional set of regulators to exclude is provided (via exclude_regulators),
those rows are removed before computing counts.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
import logging
from typing import Union, Optional
import matplotlib.patches as mpatches

logger = logging.getLogger(__name__)

def plot_tfbs_counts(source_csv: Union[str, Path],
                     mapping_csv: Union[str, Path],
                     output_file: str,
                     figsize: tuple = (10, 6),
                     default_color: str = "#d3d3d3",
                     global_regulator_color: str = "#a9a9a9",
                     sigma_factor_color: str = "#696969",
                     log_scale: bool = False,
                     exclude_regulators: Optional[set] = None) -> None:
    """
    Generate a bar plot of TF counts for a given TFBS source CSV.
    """
    # Load source data.
    df_source = pd.read_csv(source_csv)
    if "tf" not in df_source.columns or "tfbs" not in df_source.columns:
        raise ValueError(f"'tf' or 'tfbs' column missing in {source_csv}")
    
    # If regulators are to be excluded, filter them out.
    if exclude_regulators:
        df_source = df_source[~df_source["tf"].isin(exclude_regulators)]
    
    # Count unique tfbs per tf.
    counts = df_source.groupby("tf")["tfbs"].nunique().sort_values(ascending=False)
    
    # Load mapping data.
    df_map = pd.read_csv(mapping_csv)
    if "tf" not in df_map.columns:
        raise ValueError(f"'tf' column missing in {mapping_csv}")
    df_map_unique = df_map.drop_duplicates(subset=["tf"])
    mapping_dict = df_map_unique.set_index("tf").to_dict(orient="index")
    
    # Determine bar colors.
    bar_colors = []
    for tf in counts.index:
        color = default_color
        info = mapping_dict.get(tf, {})
        if info.get("is_sigma_factor", "").strip().lower() == "yes":
            color = sigma_factor_color
        elif info.get("is_global_regulator", "").strip().lower() == "yes":
            color = global_regulator_color
        bar_colors.append(color)
    
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=figsize)
    counts.plot(kind="bar", color=bar_colors, ax=ax)
    
    if log_scale:
        ax.set_yscale("log")
    
    # Title from parent folder name.
    parent_name = Path(source_csv).parent.name
    ax.set_title(f"TFBS Counts for {parent_name}")
    ax.set_xlabel("Transcription Factor")
    ax.set_ylabel("Unique TFBS Count")
    ax.tick_params(axis='x', rotation=90)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Create legend (without frame).
    patches = [
        mpatches.Patch(color=default_color, label="General"),
        mpatches.Patch(color=global_regulator_color, label="Global Regulator"),
        mpatches.Patch(color=sigma_factor_color, label="Sigma Factor")
    ]
    ax.legend(handles=patches, loc="upper right", title="TF Category", frameon=False)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"TFBS counts plot saved to {output_file}")
