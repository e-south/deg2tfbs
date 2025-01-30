"""
--------------------------------------------------------------------------------
<deg2tfbs project>
radzikowski.py

Module for loading and analyzing data from Radzikowski et al., which collected E. coli
proteomic datasets across different growth conditions to identify differentially expressed
in the context of growth limiation and persistence.

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and saves an MA plot (average expression vs. log2 Fold Change).

"Bacterial persistence is an active σS stress response to metabolic flux limitation"
DOI: 10.15252/msb.20166998

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import yaml
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_radzikowski_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Radzikowski dataset using keys from the YAML config.
    
    E.g. config:
      {
        dataset_key: "radzikowski"
        sheet_name: "Sheet1"
        usecols: "A:S"
        header: 2
      }
      
    Note: Since the actual file columns include 'glucose growth', 
          'persisters, 0.5h', 'starved, 0.5h', etc., we rename 
          them here to match the keys used in the YAML 'comparisons':
            "Persisters_0.5h"
            "GlucoseGrowth"
            "Starved_0.5h"
            ...
    """
    # Load dataset
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )

    # Rename columns so they match the "Persisters_0.5h", "GlucoseGrowth", etc. 
    # that the pipeline references in the YAML comparisons.
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

        # Calculate log2FC and average expression
        df_sub = df[["Gene", numerator, denominator]].dropna().copy()
        df_sub[f"log2FC_{label}"] = np.log2(df_sub[numerator] + eps) - np.log2(df_sub[denominator] + eps)
        df_sub[f"avg_exp_{label}"] = (df_sub[numerator] + df_sub[denominator]) / 2

        # Up- and down-regulated genes are colored red and green, respectively
        df_sub["color"] = np.where(
            df_sub[f"log2FC_{label}"] >= threshold, "red",
            np.where(df_sub[f"log2FC_{label}"] <= -threshold, "green", "gray")
        )
        
        # Classify up/down based on threshold
        up = df_sub[df_sub[f"log2FC_{label}"] >= threshold].copy()
        up["comparison"] = label
        down = df_sub[df_sub[f"log2FC_{label}"] <= -threshold].copy()
        down["comparison"] = label
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

    required_columns = ["gene", "source", "thresholds", "comparison"]
    
    # Tidy up the DataFrames of up- and down-regulated genes
    up_clean = up_all.copy().rename(columns={"Gene": "gene"})
    up_clean["source"] = "radzikowski"
    up_clean["thresholds"] = threshold
    up_clean = up_clean[required_columns]

    down_clean = down_all.copy().rename(columns={"Gene": "gene"})
    down_clean["source"] = "radzikowski"
    down_clean["thresholds"] = threshold
    down_clean = down_clean[required_columns]

    up_clean.to_csv(csv_dir / "radzikowski_upregulated_degs.csv", index=False)
    down_clean.to_csv(csv_dir / "radzikowski_downregulated_degs.csv", index=False)

    print(f"[Radzikowski et al. Pipeline] Completed. Identified DEGs across {len(comparisons)} condition pairs at log2 ≥ {threshold}: {len(up_all)} up, {len(down_all)} down.")  

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_radzikowski_pipeline(full_config)
