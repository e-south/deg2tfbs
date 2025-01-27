"""
--------------------------------------------------------------------------------
<deg2tfbs project>
schmidt.py

Module for reading and analyzing data from Schmidt et al., Table S6, which contains 
global absolute abundance estimations from both datasets including functional 
annotations using cluster of orthologues groups (COG).

"The quantitative and condition-dependent Escherichia coli proteome"
DOI: 10.1038/nbt.3418
Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_schmidt_data(config_data: dict) -> pd.DataFrame:
    """
    Loads Schmidt data from the YAML config, e.g.:

      dataset_key: "schmidt"
      sheet_name: "Schmidt et al."
      usecols: "C:AC"
      header: 2

    Replaces inf/-inf with NaN. Returns a DataFrame.
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    return df

def schmidt_comparisons(df, comparisons, threshold=2.0, save_plots=False, plot_dir=None):
    """
    For each comparison label in 'comparisons', we do:
      - log2FC = log2( df[condition] + 1 ) - log2( df[reference] + 1 )
      - average = ( condition + reference ) / 2
      - classify up/down using threshold

    Returns a combined up/down DataFrame across all comparisons.
    """
    all_up = []
    all_down = []

    for label, (condition, reference) in comparisons.items():
        # Avoid zero with small offset
        offset = 1.0
        df_sub = df[[condition, reference]].dropna().copy()
        df_sub["log2FC"] = np.log2(df_sub[condition] + offset) - np.log2(df_sub[reference] + offset)
        df_sub["avg_exp"] = (df_sub[condition] + df_sub[reference]) / 2

        df_sub["color"] = np.where(
            df_sub["log2FC"] >= threshold, "red",
            np.where(df_sub["log2FC"] <= -threshold, "green", "gray")
        )
        up = df_sub[df_sub["log2FC"] >= threshold]
        down = df_sub[df_sub["log2FC"] <= -threshold]
        all_up.append(up)
        all_down.append(down)

        # Optional plot
        if save_plots and plot_dir is not None:
            plt.figure(figsize=(6,5))
            plt.scatter(df_sub["avg_exp"], df_sub["log2FC"], c=df_sub["color"], alpha=0.25, edgecolors="none")
            plt.axhline(threshold, color="gray", linestyle="--")
            plt.axhline(-threshold, color="gray", linestyle="--")
            plt.xscale("log")
            plt.title(f"Schmidt comparison: {label} ({condition} vs {reference})")
            plt.xlabel("Average Protein Copies/Cell")
            plt.ylabel("Log2 Fold Change")
            sns.despine()
            plot_file = plot_dir / f"schmidt_{label}.png"
            plt.savefig(plot_file, dpi=150)
            plt.close()

    up_all = pd.concat(all_up, ignore_index=True).drop_duplicates()
    down_all = pd.concat(all_down, ignore_index=True).drop_duplicates()
    return up_all, down_all

def run_schmidt_pipeline(full_config: dict):
    """
    Expects in YAML:
      schmidt:
        data:
          dataset_key: "schmidt"
          sheet_name: "Schmidt et al."
          usecols: "C:AC"
          header: 2
        thresholds:
          log2_fc_threshold: 2.0
          comparisons:
            "Acetate_vs_Glucose": [ "Acetate", "Glucose" ]
            "Xylose_vs_Glucose": [ "Xylose", "Glucose" ]
            ...
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """
    config_schmidt = full_config["schmidt"]
    df = read_schmidt_data(config_schmidt["data"])

    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_schmidt["output"]["csv_subdir"]
    plot_dir = batch_dir / config_schmidt["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    threshold = config_schmidt["thresholds"]["log2_fc_threshold"]
    comparisons = config_schmidt["thresholds"]["comparisons"]

    up_all, down_all = schmidt_comparisons(
        df,
        comparisons=comparisons,
        threshold=threshold,
        save_plots=True,
        plot_dir=plot_dir
    )

    up_all.to_csv(csv_dir / "DEGs_upregulated_all_schmidt.csv", index=False)
    down_all.to_csv(csv_dir / "DEGs_downregulated_all_schmidt.csv", index=False)

    print(f"[Schmidt Pipeline] Completed. Found: {len(up_all)} up, {len(down_all)} down total.")
