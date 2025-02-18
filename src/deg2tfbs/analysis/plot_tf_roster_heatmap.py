"""
--------------------------------------------------------------------------------
<deg2tfbs project>
src/deg2tfbs/analysis/plot_tf_roster_clustering.py

Generates a heatmap from TF roster data. Each roster (binary vector) is a row in the
heatmap (columns represent regulators). Rows are ordered based on a provided clustering
(e.g. Leiden clustering labels) with the reference source placed at the top.

Each row is colored based on its Leiden cluster. For each cell, a 0 is rendered 
in a light (pastel) version of the cluster’s color and a 1 is rendered in the full 
(primary) version.
  
Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.colors as mcolors
from pathlib import Path
import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by blending it with white.
    
    Input can be a matplotlib color string, hex string, or RGB tuple.
    """
    try:
        c = mcolors.cnames[color]
    except KeyError:
        c = color
    c = mcolors.to_rgb(c)
    white = np.array([1, 1, 1])
    c_light = np.array(c) + (white - np.array(c)) * amount
    return tuple(c_light)

def trim_label(label: str) -> str:
    """
    Trim the given label so that only the text after the second underscore is returned.
    If there are fewer than 2 underscores, return the original label.
    """
    parts = label.split("_")
    if len(parts) > 2:
        return "_".join(parts[2:])
    return label

def plot_heatmap(rosters: Dict[str, np.array],
                 regs_reference: List[str],
                 cluster_mapping: Dict[str, str],
                 output_file: str,
                 figsize: tuple = (20, 8),
                 xlabel: str = "Transcription Factor",
                 ylabel: str = "Enriched Experimental Condition",
                 title: str = "Sets of transcription factors associated with isolated DEGs",
                 cell_edge_color: str = "white",
                 reference_set: Optional[str] = None) -> None:
    """
    Generate a heatmap of TF rosters.

    Parameters:
      - rosters: Dictionary mapping source name to binary vector (np.array).
      - regs_reference: List of regulators (defines the columns).
      - cluster_mapping: Dictionary mapping source name to cluster label.
      - output_file: File path to save the plot.
      - figsize: Figure size.
      - xlabel: Label for the x-axis (default: "Transcription Factor").
      - ylabel: Label for the y-axis (default: "Enriched Experimental Condition").
      - title: Plot title (default: "Sets of transcription factors associated with isolated DEGs").
      - cell_edge_color: Color for the outline of each cell.
      - reference_set: Name of the reference source. If provided, its row is forced to be at the top.
    """
    # Ensure every roster has a cluster label.
    for key in rosters.keys():
        if key not in cluster_mapping:
            logger.warning(f"No cluster information for {key}; assigning default cluster '0'.")
            cluster_mapping[key] = '0'
    
    # Order sources: if reference_set is provided, place it first.
    def sort_key(k):
        if reference_set is not None and k == reference_set:
            return (-1, trim_label(k))
        return (int(cluster_mapping.get(k, '0')), trim_label(k))
    
    ordered_sources = sorted(rosters.keys(), key=sort_key)
    matrix = np.vstack([rosters[src] for src in ordered_sources])
    n_rows, n_cols = matrix.shape

    # Determine unique clusters and assign each a primary color using a colorblind-friendly palette.
    unique_clusters = sorted(set(cluster_mapping.values()), key=lambda x: int(x))
    palette = sns.color_palette("colorblind", len(unique_clusters))
    cluster_color_map = {cl: palette[i] for i, cl in enumerate(unique_clusters)}
    # Compute a pastel (light) version for each cluster.
    cluster_color_pair = {cl: {"light": lighten_color(cluster_color_map[cl], amount=0.7),
                               "dark": cluster_color_map[cl]}
                          for cl in unique_clusters}
    
    # Create a color matrix of shape (n_rows, n_cols).
    color_matrix = np.empty((n_rows, n_cols), dtype=object)
    for i, src in enumerate(ordered_sources):
        cl = cluster_mapping.get(src, '0')
        colors = cluster_color_pair[cl]
        for j in range(n_cols):
            color_matrix[i, j] = colors["dark"] if matrix[i, j] == 1 else colors["light"]
    
    # Plot the colored grid.
    fig, ax = plt.subplots(figsize=figsize)
    for i in range(n_rows):
        for j in range(n_cols):
            rect = plt.Rectangle([j, i], 1, 1, facecolor=color_matrix[i, j],
                                 edgecolor=cell_edge_color, linewidth=0.5)
            ax.add_patch(rect)
    
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.invert_yaxis()  # So that row 0 is at the top.
    ax.set_aspect("equal")  # Ensure each cell is square.
    
    # Set tick labels.
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(regs_reference, rotation=90, fontsize=8, fontweight="normal")
    
    y_labels = [trim_label(src) for src in ordered_sources]
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(y_labels, rotation=0, fontsize=8, fontweight="normal")
    
    # Set emboldened axis labels and title.
    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    ax.set_title(title, fontweight="bold")
    
    # Remove spines.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"TF roster heatmap saved to {output_file}")
