"""
--------------------------------------------------------------------------------
<deg2tfbs project>
radzikowski.py

Module for reading and analyzing data from Radzikowski et al.. Radzikowski et al. 
developed and experimentally verified a model, in which persistence is established 
through a system-level feedback: Strong perturbations of metabolic homeostasis 
cause metabolic fluxes to collapse, prohibiting adjustments toward restoring 
homeostasis.

"Bacterial persistence is an active σS stress response to metabolic flux limitation"
DOI: 10.15252/msb.20166998

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

def read_radzikowski_data(config_data: dict) -> pd.DataFrame:
    """
    Loads Radzikowski data using the user-provided 'dataset_key', etc. from YAML.

    Since the actual file columns include 'glucose growth', 'persisters, 0.5h', 
    'starved, 0.5h', etc., we rename them here to match the keys 
    used in the YAML 'comparisons':

       "Persisters_0.5h"
       "GlucoseGrowth"
       "Starved_0.5h"
       ...

    Example config_data might be:
      dataset_key: "radzikowski"
      sheet_name: "Sheet1"
      usecols: "A:S"
      header: 2
    """
    # 1) Load dataset
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )

    # 2) Rename columns so they match the "Persisters_0.5h", "GlucoseGrowth", etc. 
    #    that the pipeline references in the YAML comparisons.
    column_renames = {
        "Accession": "Accession",
        "Gene": "Gene",
        "glucose growth": "GlucoseGrowth",
        "persisters, 0.5h": "Persisters_0.5h",
        "starved, 0.5h": "Starved_0.5h",
        "glucose growth": "GlucoseGrowth",
        "persisters, 0.5h": "Persisters_0.5h",
        "persisters, 1h": "Persisters_1h",
        "starved, 0.5h": "Starved_0.5h",
        "starved, 1h": "Starved_1h",
    }

    # Apply partial renaming. Columns not listed remain unchanged.
    df.rename(columns=column_renames, inplace=True)

    return df

def radzikowski_comparisons(df, comparisons, threshold=2.0, save_plots=False, plot_dir=None):
    """
    For each label in `comparisons`, do:
      log2FC = log2(numerator + eps) - log2(denominator + eps)
      avg_exp = (numerator + denominator) / 2
    Classify up/down by threshold.

    Returns aggregated up/down DataFrames.
    """
    all_up = []
    all_down = []

    eps = 1e-5
    for label, (numerator, denominator) in comparisons.items():
        # Ensure we only handle rows that have both columns
        if numerator not in df.columns or denominator not in df.columns:
            print(f"Skipping comparison '{label}' because columns '{numerator}' or '{denominator}' not in DataFrame.")
            continue

        df_sub = df[[numerator, denominator]].dropna().copy()
        df_sub[f"log2FC_{label}"] = np.log2(df_sub[numerator] + eps) - np.log2(df_sub[denominator] + eps)
        df_sub[f"avg_exp_{label}"] = (df_sub[numerator] + df_sub[denominator]) / 2

        df_sub["color"] = np.where(
            df_sub[f"log2FC_{label}"] >= threshold, "red",
            np.where(df_sub[f"log2FC_{label}"] <= -threshold, "green", "gray")
        )

        up = df_sub[df_sub[f"log2FC_{label}"] >= threshold]
        down = df_sub[df_sub[f"log2FC_{label}"] <= -threshold]
        all_up.append(up)
        all_down.append(down)

        # Optional plot
        if save_plots and plot_dir is not None:
            plt.figure(figsize=(6,5))
            plt.scatter(df_sub[f"avg_exp_{label}"], df_sub[f"log2FC_{label}"],
                        c=df_sub["color"], alpha=0.25, edgecolors="none")
            plt.axhline(threshold, color="gray", linestyle="--")
            plt.axhline(-threshold, color="gray", linestyle="--")
            plt.xscale("log")
            plt.title(f"Radzikowski et al.\n log2({label})")
            plt.xlabel("Average Abundance")
            plt.ylabel("Log2 Fold Change")
            sns.despine()
            plot_file = plot_dir / f"radzikowski_{label}.png"
            plt.savefig(plot_file, dpi=150)
            plt.close()

    up_all = pd.concat(all_up, ignore_index=True).drop_duplicates() if all_up else pd.DataFrame()
    down_all = pd.concat(all_down, ignore_index=True).drop_duplicates() if all_down else pd.DataFrame()
    return up_all, down_all

def run_radzikowski_pipeline(full_config: dict):
    """
    Expects in YAML something like:
      radzikowski:
        data:
          dataset_key: "radzikowski"
          sheet_name: "Sheet1"
          usecols: "A:S"
          header: 2
        thresholds:
          log2_fc_threshold: 2.0
          comparisons:
            Persisters_0.5h_Over_GlucoseGrowth: ["Persisters_0.5h", "GlucoseGrowth"]
            Starved_0.5h_Over_GlucoseGrowth: ["Starved_0.5h", "GlucoseGrowth"]
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """
    config_radz = full_config["radzikowski"]
    df = read_radzikowski_data(config_radz["data"])

    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_radz["output"]["csv_subdir"]
    plot_dir = batch_dir / config_radz["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    threshold = config_radz["thresholds"]["log2_fc_threshold"]
    comparisons = config_radz["thresholds"]["comparisons"]

    up_all, down_all = radzikowski_comparisons(
        df,
        comparisons=comparisons,
        threshold=threshold,
        save_plots=True,
        plot_dir=plot_dir
    )

    # Only save if we found any results
    if not up_all.empty:
        up_all.to_csv(csv_dir / "DEGs_upregulated_all_radzikowski.csv", index=False)
    if not down_all.empty:
        down_all.to_csv(csv_dir / "DEGs_downregulated_all_radzikowski.csv", index=False)

    print(f"[Radzikowski Pipeline] Completed. Found: {len(up_all)} up, {len(down_all)} down total.")
