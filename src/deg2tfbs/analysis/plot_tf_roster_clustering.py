"""
--------------------------------------------------------------------------------
<deg2tfbs project>
src/deg2tfbs/analysis/plot_tf_roster_clustering.py

Provides functions to run a Scanpy clustering (UMAP + Leiden) workflow and
to generate a UMAP scatter plot of the TF roster data.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="anndata.utils")

import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull
import numpy as np
import logging
from typing import Tuple, List, Optional

logger = logging.getLogger(__name__)

def run_scanpy_clustering(matrix: np.array,
                          sample_names: List[str],
                          n_neighbors: int = 15,
                          min_dist: float = 0.1,
                          random_state: int = 0) -> anndata.AnnData:
    """
    Runs the clustering workflow using Scanpy.
    
    Preconditions:
      - matrix is a 2D numpy array.
      - len(sample_names) matches the number of rows in matrix.
    
    Returns:
      - anndata.AnnData with UMAP coordinates and Leiden clustering in obs["leiden"].
    """
    assert matrix.ndim == 2, "Input matrix must be 2D."
    assert matrix.shape[0] == len(sample_names), "Mismatch between matrix rows and sample_names length."
    
    adata = anndata.AnnData(X=matrix)
    adata.obs["sample"] = sample_names
    sc.pp.pca(adata, n_comps=10, random_state=random_state)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=random_state)
    sc.tl.umap(adata, min_dist=min_dist, random_state=random_state)
    sc.tl.leiden(adata, random_state=random_state)
    logger.info("Scanpy clustering complete.")
    return adata

def custom_umap_plot(ax, adata, title: str, alpha: float = 0.8,
                     base_size: int = 30, scale_factor: int = 2,
                     show_legend: bool = False,
                     reference_set: Optional[str] = None):
    """
    Draws a UMAP scatter plot on a given Axes.

    For each Leiden cluster (excluding reference samples), this function assigns 
    unique marker shapes to each sample (if possible) so that no shape is repeated 
    within that cluster. The marker sizes for the scatter points are scaled based 
    on the sample's reg_count. However, the legend is constructed separately with 
    fixed marker sizes.

    Preconditions:
      - adata.obs must contain 'leiden', 'sample', and 'reg_count'.
    """
    sns.set_style("ticks")
    X = adata.obsm["X_umap"]
    
    for col in ["leiden", "sample", "reg_count"]:
        if col not in adata.obs.columns:
            raise ValueError(f"{col} column not found in adata.obs")
    
    # Get unique Leiden clusters (sorted numerically).
    clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
    # Assign colors to clusters using a colorblind-friendly palette.
    cluster_colors = {cl: col for cl, col in zip(clusters, sns.color_palette("colorblind", len(clusters)))}
    
    # Define a list of marker shapes.
    markers = ['o', 's', '^', 'v', 'D', 'P', 'X', '*', '+', 'p']
    
    # Dictionary to hold custom legend handles.
    legend_dict = {}
    
    # For each cluster (non-reference), assign markers uniquely per sample.
    for cl in clusters:
        idx_cluster = (adata.obs["leiden"] == cl)
        if reference_set is not None:
            idx_cluster = idx_cluster & (adata.obs["sample"] != reference_set)
        if idx_cluster.sum() == 0:
            continue
        # Get unique sample names in this cluster.
        cluster_samples = adata.obs.loc[idx_cluster, "sample"].unique()
        cluster_marker_dict = {sample: markers[i % len(markers)] for i, sample in enumerate(cluster_samples)}
        
        # For each sample in this cluster, plot its points.
        for sample in cluster_samples:
            idx_sample = idx_cluster & (adata.obs["sample"] == sample)
            if idx_sample.sum() == 0:
                continue
            sizes = base_size + scale_factor * adata.obs.loc[idx_sample, "reg_count"].astype(float)
            ax.scatter(X[idx_sample, 0], X[idx_sample, 1],
                       color=cluster_colors[cl],
                       marker=cluster_marker_dict[sample],
                       alpha=alpha,
                       s=sizes)
            label = f"Cluster {cl}: {sample}"
            # Store a legend handle with a fixed marker size.
            if label not in legend_dict:
                legend_dict[label] = Line2D([], [], marker=cluster_marker_dict[sample],
                                              color=cluster_colors[cl],
                                              linestyle='', markersize=8)
    
    # Overplot reference points (if provided) in bright red.
    if reference_set is not None:
        idx_ref = adata.obs["sample"] == reference_set
        if idx_ref.sum() > 0:
            sizes_ref = base_size + scale_factor * adata.obs.loc[idx_ref, "reg_count"].astype(float)
            ax.scatter(X[idx_ref, 0], X[idx_ref, 1],
                       color="#ff0000", marker="o",
                       s=sizes_ref, edgecolor="black")
            legend_dict["Reference"] = Line2D([], [], marker="o", color="#ff0000",
                                              linestyle='', markersize=10)
    
    # Draw convex hulls for clusters (non-reference).
    for cl in clusters:
        idx = (adata.obs["leiden"] == cl)
        if reference_set is not None:
            idx = idx & (adata.obs["sample"] != reference_set)
        if idx.sum() < 3:
            continue  # Need at least 3 points.
        points = X[idx]
        try:
            hull = ConvexHull(points)
            hull_points = points[hull.vertices]
            hull_points = np.concatenate([hull_points, hull_points[0:1]], axis=0)
            ax.plot(hull_points[:, 0], hull_points[:, 1],
                    linestyle="--", color=cluster_colors[cl], linewidth=1)
        except Exception as e:
            logger.warning(f"Could not compute convex hull for cluster {cl}: {e}")
    
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(title)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if show_legend:
        handles = list(legend_dict.values())
        labels = list(legend_dict.keys())
        ax.legend(handles, labels, title="Cluster: Sample",
                  bbox_to_anchor=(1.05, 1), loc='upper left',
                  frameon=False, fontsize=7, title_fontsize=8)
    
    return ax

def plot_umap(adata, title: str, output_file: str,
              figsize: tuple = (12, 6), alpha: float = 0.8,
              base_size: int = 30, scale_factor: int = 2, dpi: int = 300,
              reference_set: Optional[str] = None) -> None:
    """
    Creates a UMAP plot using custom_umap_plot and saves it to output_file.
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    custom_umap_plot(ax, adata, title, alpha=alpha, base_size=base_size,
                     scale_factor=scale_factor, show_legend=True,
                     reference_set=reference_set)
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches="tight", dpi=dpi)
    plt.close()
    logger.info(f"UMAP plot saved to {output_file}")
