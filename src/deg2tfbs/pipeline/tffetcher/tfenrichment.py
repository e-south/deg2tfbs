#!/usr/bin/env python
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
  - A Fisher’s exact test is performed using the following 2×2 contingency table:
  
             |   DEG    | Non-DEG 
    ---------|----------|---------
    Targets  |    a     |  K - a
    Non-targets |  M - a |  N - K - (M - a)
    
    where:
      M = number of DEGs (present in the full network)
      N = total number of genes in the full network.
      
  - P-values are corrected using the Benjamini–Hochberg method.

The module outputs:
  1. A CSV summary (tf_enrichment_summary.csv) that includes for each TF:
       regulator, num_degs_regulated, total_targets, enrichment_score, p_value, fdr, is_significant,
       up_appearance_count, down_appearance_count, is_global_regulator, is_sigma_factor, is_nucleoid_regulator,
       regulator_type, a, K, M, N, and a boolean "topN" indicating whether the TF is among the top-N
       (lowest FDR among non-sigma, non-nucleoid candidates).
  2. A consolidated enrichment plot with two subplots sharing the same x-axis:
       - Top subplot: a barplot of enrichment scores with y-axis label 
         "Enrichment Score (DEG targets / total targets)". Each bar is annotated with a two‐line text:
         the first line shows the number of DEG targets (a), a horizontal line ("────") is drawn, and the second line
         shows the total targets (K). A composite annotation is placed in the top right explaining these values.
       - Bottom subplot: a barplot of negative FDR values. Top-N candidates are colored in darker pastel red (#ff5555)
         and the remaining ones in light gray (#d3d3d3). A vertical dashed line (in the same red hue) marks the top-N partition,
         with a rotated annotation "Top {N}" (without a colon) positioned to the left of the line.
       - X-tick labels are uniformly black, with a slight margin between the y-axis and the first x tick.
       
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
from matplotlib.patches import Patch


# Set random seed for reproducibility.
np.random.seed(42)

# Import regulator classification utilities.
from deg2tfbs.pipeline.tffetcher.regulator_utils import (
    is_master_regulator,
    is_sigma_factor,
    is_nucleoid_regulator
)

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
    
    # Plotting: create two subplots.
    fig, (ax_top, ax_bottom) = plt.subplots(nrows=2, sharex=True, figsize=plot_figsize,
                                            gridspec_kw={"height_ratios": [3, 1]})
    # Provide outer margins.
    plt.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.12, hspace=0.05)
    x = np.arange(len(df_plot))
    bar_width = 0.6

    # TOP SUBPLOT: Plot barplot of enrichment scores.
    bars = ax_top.bar(x, df_plot["enrichment_score"], width=bar_width, color="#7f7f7f", alpha=0.8)
    ax_top.set_ylabel("Enrichment Score (DEG targets / total targets)", fontsize=14)
    ax_top.set_title("TF Target Enrichment in DEGs: Scores and BH-FDR Adjusted Significance", fontsize=16)
    # Increase margin between title and plot by shifting title upward.
    ax_top.title.set_position([0.5, 1.20])
    ax_top.tick_params(axis="both", labelsize=12)
    # Annotate each bar with a two-line annotation: top line = DEG targets count, bottom line = total targets.
    for i, row in df_plot.iterrows():
        es = row["enrichment_score"]
        annotation = f"{row['a']}\n─\n{row['K']}"
        ax_top.text(x[i], es + 0.01, annotation, ha="center", va="bottom", fontsize=10, color="gray")
    # Composite legend (custom annotation) placed in the top right, vertically centered.
    composite_text = "DEG targets for TF\n────────────\nTotal number of targets for TF"
    ax_top.text(0.88, 0.5, composite_text, transform=ax_top.transAxes, ha="center", va="center",
                fontsize=14, color="black", bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))
    sns.despine(ax=ax_top, top=True, right=True)
    
    # BOTTOM SUBPLOT: Plot negative FDR values.
    def get_bar_color(topN):
        return "#ff5555" if topN else "#d3d3d3"
    bar_colors = [get_bar_color(top) for top in df_plot["topN"]]
    ax_bottom.bar(x, -df_plot["fdr"], width=bar_width, color=bar_colors, alpha=0.8)
    ax_bottom.set_ylabel("BH-FDR corrected p-value", fontsize=14)
    ax_bottom.tick_params(axis="both", labelsize=12)
    ax_bottom.set_ylim(-1, 0)
    ax_bottom.set_yticks(np.linspace(-1, 0, 6))
    ax_bottom.set_yticklabels([f"{abs(y):.2f}" for y in np.linspace(-1, 0, 6)], fontsize=12)
    ax_bottom.set_xlabel("Transcription Factor", fontsize=14)
    
    # Provide a slight margin between the y-axis and the first x tick.
    ax_bottom.margins(x=0.05)
    
    # Draw vertical dashed line to demarcate top-N partition.
    if len(df_plot) >= top_n:
        partition_x = top_n - 0.5
        ax_bottom.axvline(x=partition_x, color="#ff5555", linestyle="--", linewidth=1)
        # Annotate the vertical line with rotated "Top {N}" (without colon) positioned to the left.
        ax_bottom.text(partition_x - 0.2, -0.95, f"Top {top_n}", color="#ff5555",
                       fontsize=14, ha="right", va="bottom", rotation=90)
    
    # Set x-axis tick labels uniformly in black.
    ax_bottom.set_xticks(x)
    ax_bottom.set_xticklabels(df_plot["regulator"], rotation=90, fontsize=12, color="black")
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
