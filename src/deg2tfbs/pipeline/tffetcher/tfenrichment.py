"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tfenrichment.py

Author(s): Eric J. South
Dunlop Lab

This module performs enrichment analysis for transcription factors (TFs)
by comparing the number of differentially expressed genes (DEGs) targeted by each TF
to the total number of targets for that TF (from a full regulatory network).

For each TF:
  - Let a = number of DEGs (from the union of all DEG files) that are targets of the TF.
  - Let K = total number of targets for the TF (from the full network).
  - The enrichment score is defined as a/K.
  - A Fisher's exact test is performed using the following 2x2 contingency table:
  
                |   DEG    | Non-DEG 
    ------------|----------|---------
    Targets     |    a     |  K - a
    Non-targets |  M - a   |  N - K - (M - a)
    
    where:
      M = number of DEGs (present in the full network)
      N = total number of genes in the full network.
      
  - P-values are corrected using the Benjamini-Hochberg method.

The module outputs:
  1. A CSV summary (tf_enrichment_summary.csv) that includes for each TF:
       regulator, num_degs_regulated, total_targets, enrichment_score, p_value, fdr, is_significant,
       up_appearance_count, down_appearance_count, is_global_regulator, is_sigma_factor, is_nucleoid_regulator,
       regulator_type, a, K, M, N, and a boolean "topN" indicating whether the TF is among the top-N
       (lowest FDR among non-sigma, non-nucleoid candidates).
  2. A consolidated enrichment plot with two subplots sharing the same x-axis:
       - Top subplot: a scatter plot of enrichment scores.
         * x-axis: Transcription factors (TFs) ranked by p-value (from the bottom subplot).
         * y-axis: Enrichment score (DEG targets / total targets).
         * Dot size encodes total targets (K) for each TF.
         * Dot color is red for top-N regulators and gray for the others.
         * A legend shows that point size equates to total gene targets.
       - Bottom subplot: a barplot of negative FDR values.
         Top-N candidates are colored in darker pastel red (#ff5555) and the remaining ones in light gray (#d3d3d3).
         A vertical dashed line (in the same red hue) marks the top-N partition, with a rotated annotation "Top {N}".
       - X-tick labels are uniformly black and shared between subplots.
       
TF Target Enrichment in DEGs: Scores and BH-FDR Adjusted Significance
--------------------------------------------------------------------------------
"""

import logging
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter

# Set random seed for reproducibility.
np.random.seed(42)

# Import regulator classification utilities.
from deg2tfbs.pipeline.tffetcher.regulator_utils import (
    is_master_regulator,
    is_sigma_factor,
    is_nucleoid_regulator
)

from deg2tfbs.pipeline.utils import get_regulator_display_name

logger = logging.getLogger("tfenrichment")

def run_enrichment(merged_network, mapping_rows, all_degs, out_dir: Path, params: dict, nuc_reg: set):
    # BACKGROUND
    background_genes = set(merged_network.keys())
    N = len(background_genes)
    M = len(all_degs.intersection(background_genes))
    logger.info(f"[tfenrichment] Background gene count (N): {N}")
    logger.info(f"[tfenrichment] DEGs intersecting background (M): {M}")
    
    # Build regulator-background mapping.
    regulator_background = defaultdict(set)
    for gene, reg_dict in merged_network.items():
        for regulator in reg_dict:
            regulator_background[regulator].add(gene)
    
    # Count DEG associations.
    regulator_deg_count = defaultdict(int)
    regulator_up_count = defaultdict(int)
    regulator_down_count = defaultdict(int)
    for row in mapping_rows:
        reg = row["regulator"]
        regulator_deg_count[reg] += 1
        deg_source = row.get("deg_source", "")
        if "_up" in deg_source:
            regulator_up_count[reg] += 1
        elif "_down" in deg_source:
            regulator_down_count[reg] += 1
    
    results = []
    pvals = []
    for reg in regulator_background.keys():
        a = regulator_deg_count.get(reg, 0)
        K = len(regulator_background[reg])
        if K == 0:
            p_value = 1.0
            es = 0.0
        else:
            b = K - a
            c = M - a
            d = N - K - (M - a)
            table = [[a, b], [c, d]]
            try:
                _, p_value = fisher_exact(table, alternative='greater')
            except Exception as e:
                logger.error(f"Error computing Fisher's test for regulator {reg}: {e}")
                p_value = 1.0
            es = a / K if K > 0 else 0.0
        pvals.append(p_value)
        results.append({
            "regulator": reg,
            "num_degs_regulated": a,
            "total_targets": K,
            "enrichment_score": es,
            "p_value": p_value,
            "up_appearance_count": regulator_up_count.get(reg, 0),
            "down_appearance_count": regulator_down_count.get(reg, 0),
            "is_global_regulator": "yes" if is_master_regulator(reg) else "no",
            "is_sigma_factor": "yes" if is_sigma_factor(reg) else "no",
            "is_nucleoid_regulator": "yes" if is_nucleoid_regulator(reg, nuc_reg) else "no",
            "a": a,
            "K": K,
        })
    
    if pvals:
        _, fdrs, _, _ = multipletests(pvals, method="fdr_bh")
    else:
        fdrs = [1.0] * len(results)
    
    top_n = params.get("top_n", 25)
    for i, res in enumerate(results):
        res["fdr"] = fdrs[i] if i < len(fdrs) else 1.0
        res["M"] = M
        res["N"] = N
    df_enrich = pd.DataFrame(results)
    df_enrich.sort_values(by="fdr", ascending=True, inplace=True)
    
    # Ensure "topN" column exists; if not, create it.
    if "topN" not in df_enrich.columns:
        df_enrich["topN"] = False

    # Identify candidates for top-N (non-sigma, non-nucleoid with DEG counts)
    df_candidates = df_enrich[(df_enrich["num_degs_regulated"] > 0) &
                              (df_enrich["is_sigma_factor"] != "yes") &
                              (df_enrich["is_nucleoid_regulator"] != "yes")].copy()
    df_candidates.sort_values(by="fdr", ascending=True, inplace=True)
    df_candidates["topN"] = False
    if not df_candidates.empty:
        top_indices = df_candidates.head(top_n).index
        df_candidates.loc[top_indices, "topN"] = True
        logger.info("[tfenrichment] Top-N regulators (non-sigma, non-nucleoid):")
        logger.info(df_candidates.loc[top_indices, ["regulator", "num_degs_regulated", "total_targets", "enrichment_score", "fdr"]].to_string(index=False))
    df_enrich = df_enrich.merge(df_candidates[["regulator", "topN"]], on="regulator", how="left")
    if "topN" not in df_enrich.columns:
        df_enrich["topN"] = False
    else:
        df_enrich["topN"] = df_enrich["topN"].fillna(False)
    
    # For plotting, only use non-sigma, non-nucleoid candidates.
    df_plot = df_candidates.copy()
    if df_plot.empty:
        logger.warning("[tfenrichment] No non-sigma/non-nucleoid regulators with DEG counts found for plotting.")
        return df_enrich
    df_plot.reset_index(drop=True, inplace=True)
    
    # Allow user control over figure dimensions.
    plot_figsize = params.get("plot_figsize", (18, 10))
    
    # Plotting: create two subplots with shared x-axis.
    fig, (ax_top, ax_bottom) = plt.subplots(nrows=2, sharex=True, figsize=plot_figsize,
                                             gridspec_kw={"height_ratios": [3, 1]})
    # Provide outer margins.
    plt.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.12, hspace=0.05)
    x = np.arange(len(df_plot))
    
    # --- TOP SUBPLOT: Scatter plot for enrichment values ---
    # Each point represents a TF, where:
    #   - x-axis: TF index (ranked by p-value)
    #   - y-axis: Enrichment score (a/K)
    #   - Dot size encodes total targets (K)
    #   - Dot color: Red if in top-N; gray otherwise.
    sizes = df_plot["K"] * 10  # Adjust multiplier as needed.
    def get_point_color(topN):
        return "#ff5555" if topN else "#d3d3d3"
    point_colors = [get_point_color(top) for top in df_plot["topN"]]
    ax_top.scatter(x, df_plot["enrichment_score"],
                   s=sizes,
                   color=point_colors,
                   alpha=0.8)
    ax_top.set_ylabel("Enrichment Score (DEG targets / total targets)", fontsize=14)
    # Set main title and add a subtitle with a smaller font.
    ax_top.set_title("Prioritizing Transcription Factors", fontsize=16, pad=20)
    ax_top.text(0.5, 0.99, "DEG Enrichment per TF, BH-FDR Significance, and Number of Genes Regulated", 
                transform=ax_top.transAxes, ha="center", fontsize=14)
    ax_top.tick_params(axis="both", labelsize=12)
    # Ensure y-axis ticks show two decimal places
    ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    
    # --- Add custom legend for point size ---
    # Select representative K values (min, median, max)
    min_k = df_plot["K"].min()
    med_k = int(np.median(df_plot["K"]))
    max_k = df_plot["K"].max()
    # Create dummy scatter handles for the legend
    handles = [plt.scatter([], [], s=sz, color="gray", alpha=0.8) 
               for sz in [min_k*10, med_k*10, max_k*10]]
    labels = [f"{min_k} targets", f"{med_k} targets", f"{max_k} targets"]
    ax_top.legend(handles, labels, title="Total gene targets", loc="upper right", 
                  labelspacing=1.5, frameon=False)
    
    sns.despine(ax=ax_top, top=True, right=True)
    
    # --- BOTTOM SUBPLOT: Bar plot for negative FDR values ---
    def get_bar_color(topN):
        return "#ff5555" if topN else "#d3d3d3"
    bar_colors = [get_bar_color(top) for top in df_plot["topN"]]
    ax_bottom.bar(x, -df_plot["fdr"], width=0.6, color=bar_colors, alpha=0.8)
    ax_bottom.set_ylabel("BH-FDR corrected p-value", fontsize=14)
    ax_bottom.tick_params(axis="both", labelsize=12)
    ax_bottom.set_ylim(-1, 0)
    ax_bottom.set_yticks(np.linspace(-1, 0, 6))
    # Set y-axis tick labels with two decimal places (using absolute values)
    ax_bottom.set_yticklabels([f"{abs(y):.2f}" for y in np.linspace(-1, 0, 6)], fontsize=12)
    ax_bottom.set_xlabel("Transcription Factor", fontsize=14)
    
    # Provide a slight margin between the y-axis and the first x tick.
    ax_bottom.margins(x=0.05)
    
    # Draw vertical dashed line to demarcate top-N partition.
    if len(df_plot) >= top_n:
        partition_x = top_n - 0.5
        ax_bottom.axvline(x=partition_x, color="#ff5555", linestyle="--", linewidth=1)
        ax_bottom.text(partition_x - 0.2, -0.95, f"Top {top_n}", color="black",
                       fontsize=14, ha="right", va="bottom", rotation=0)
    
    # Set x-axis tick labels to show original capitalization.
    ax_bottom.set_xticks(x)
    corrected_labels = [get_regulator_display_name(reg) for reg in df_plot["regulator"]]
    ax_bottom.set_xticklabels(corrected_labels, rotation=90, fontsize=12, color="black")

    sns.despine(ax=ax_bottom, top=True, right=True)
    
    plt.tight_layout(pad=0.5)
    consolidated_plot_path = out_dir / "plots" / "tf_consolidated_enrichment.png"
    (out_dir / "plots").mkdir(parents=True, exist_ok=True)
    plt.savefig(consolidated_plot_path, dpi=600)
    plt.close()
    logger.info(f"[tfenrichment] Consolidated enrichment plot saved: {consolidated_plot_path}")
    
    return df_enrich

if __name__ == "__main__":
    # This module is intended to be run as part of the pipeline.
    pass
