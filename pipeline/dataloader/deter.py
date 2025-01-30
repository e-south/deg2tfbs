"""
--------------------------------------------------------------------------------
<deg2tfbs project>
deter.py

Module for reading and analyzing data described in Deter et al., which generated
RNA-seq data on both antibiotic-treated and -untreated populations emerging from 
stationary phase.

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and saves an MA plot (average expression vs. log2 Fold Change).

"Antibiotic tolerance is associated with a broad and complex transcriptional 
response in E. coli"
DOI: 10.1038/s41598-021-85509-7

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os

import yaml
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_deter_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Deter dataset using keys from the YAML config.
    
    E.g. config:
      {
        "dataset_key": "deter",
        "header": 0,
        "usecols": [
          "Name","gene","Stationary rep1","Stationary rep2","Stationary rep3",
          "3h fresh media rep1","3h fresh media rep2","3h fresh media rep3",
          "3h amp media rep1","3h amp media rep2","3h amp media rep3"
        ]
      }
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),  # Not used for CSV
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def deter_ma_plot(
    df_comp: pd.DataFrame,
    target_avg_col: str,
    ref_avg_col: str,
    threshold: float,
    plot_path: Path
):
    """
    Creates and saves an MA plot comparing target vs. reference columns.

    Args:
      df_comp (pd.DataFrame): Must have 'log2_fc' and 'avg_expr' columns.
      target_avg_col (str): e.g. "3h amp media_avg"
      ref_avg_col (str): e.g. "3h fresh media_avg"
      threshold (float): log2 fold-change threshold
      plot_path (Path): Where to save the figure.
    """
    sns.set_style("ticks")

    # Up- and down-regulated genes are colored red and green, respectively
    colors = np.where(
        df_comp["log2_fc"] >= threshold, "red",
        np.where(df_comp["log2_fc"] <= -threshold, "green", "gray")
    )

    plt.figure(figsize=(6,5))
    plt.scatter(
        df_comp["avg_expr"],
        df_comp["log2_fc"],
        c=colors,
        alpha=0.25,
        edgecolors="none"
    )
    plt.xscale("log")
    plt.axhline(threshold, color="gray", linestyle="--")
    plt.axhline(-threshold, color="gray", linestyle="--")

    # Extract labels from column names
    target_label = target_avg_col.replace("_avg", "")
    ref_label    = ref_avg_col.replace("_avg", "")
    plt.title(f"Deter et al.\n{target_label} vs. {ref_label}")
    plt.xlabel("Average Expression")
    plt.ylabel("Log2 Fold Change")
    sns.despine()

    plt.savefig(plot_path, dpi=150)
    plt.close()


def run_deter_pipeline(full_config: dict):
    """
    Reads Deter et al. data from CSV, 
    averages replicates for multiple comparisons,
    computes log2 FC, generating up/down CSVs, 
    and saves an MA plot for each comparison.

    Example YAML:
      deter:
        data:
          dataset_key: "deter"
          header: 0
          usecols:
            - "Name"
            - "gene"
            - "Stationary rep1"
            - "Stationary rep2"
            - "Stationary rep3"
            - "3h fresh media rep1"
            - "3h fresh media rep2"
            - "3h fresh media rep3"
            - "3h amp media rep1"
            - "3h amp media rep2"
            - "3h amp media rep3"
        thresholds:
          log2_fc_threshold: 2.0
          # comparisons is an array of dicts, each with target + reference
          comparisons:
            - { target: "3h amp media", reference: "3h fresh media" }
            - { target: "Stationary", reference: "3h fresh media" }
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """

    config_deter = full_config.get("deter", None)
    if config_deter is None:
        print("[Deter et al. Pipeline] No 'deter' config found. Skipping.")
        return

    # Load data
    df = read_deter_data(config_deter["data"])

    # Prepare output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_deter["output"]["csv_subdir"]
    plot_dir = batch_dir / config_deter["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Retrieve threshold + comparisons
    threshold = config_deter["thresholds"].get("log2_fc_threshold", 2.0)
    comparisons = config_deter["thresholds"].get("comparisons", [])

    # Helper to average replicates for a condition
    def average_condition(df_local, cond_prefix):
        # e.g. cond_prefix="3h amp media" => columns "3h amp media rep1", "3h amp media rep2", ...
        replicate_cols = [col for col in df_local.columns if col.startswith(cond_prefix)]
        if not replicate_cols:
            raise KeyError(f"No replicate columns found for '{cond_prefix}'")
        return df_local[replicate_cols].mean(axis=1)

    all_up = []
    all_down = []

    # For each comparison in the config
    for comparison in comparisons:
        target_cond = comparison["target"]
        ref_cond    = comparison["reference"]

        # Create average columns if not already present
        # e.g. "3h amp media_avg", "3h fresh media_avg"
        target_avg_col = f"{target_cond}_avg"
        ref_avg_col    = f"{ref_cond}_avg"

        if target_avg_col not in df.columns:
            df[target_avg_col] = average_condition(df, target_cond)
        if ref_avg_col not in df.columns:
            df[ref_avg_col]    = average_condition(df, ref_cond)

        # Compute log2 fold change
        eps = 1e-5
        fc_col = f"log2_fc_{target_cond}_vs_{ref_cond}"
        df[fc_col] = np.log2(
            (df[target_avg_col] + eps) / (df[ref_avg_col] + eps)
        )

        # Average expression
        avg_expr_col = f"avg_expr_{target_cond}_vs_{ref_cond}"
        df[avg_expr_col] = (df[target_avg_col] + df[ref_avg_col]) / 2.0

        # Mark up/down regulated genes
        up = df[df[fc_col] >= threshold].copy()
        down = df[df[fc_col] <= -threshold].copy()

        # Create an MA plot for this comparison
        df_comp = pd.DataFrame({
            "log2_fc": df[fc_col],
            "avg_expr": df[avg_expr_col]
        })

        plot_fname = f"deter_ma_{target_cond}_vs_{ref_cond}.png"
        deter_ma_plot(
            df_comp=df_comp,
            target_avg_col=target_avg_col,
            ref_avg_col=ref_avg_col,
            threshold=threshold,
            plot_path=plot_dir / plot_fname
        )

        required_columns = ["gene", "source", "thresholds", "comparison"]
    
        # Define comparison    
        comparison_str = f"{target_cond}_versus_{ref_cond}"
        
        # Tidy up the DataFrames of up- and down-regulated genes
        up_clean = up.copy()
        up_clean["source"] = "deter"
        up_clean["thresholds"] = threshold
        up_clean["comparison"] = comparison_str
        up_clean = up_clean[required_columns]

        down_clean = down.copy() 
        down_clean["source"] = "deter"
        down_clean["thresholds"] = threshold
        down_clean["comparison"] = comparison_str
        down_clean = down_clean[required_columns]

        # Collect for final summary
        all_up.append(up_clean)
        all_down.append(down_clean)

    # Combine all up/down across comparisons
    df_up = pd.concat(all_up, ignore_index=True).drop_duplicates()
    df_down = pd.concat(all_down, ignore_index=True).drop_duplicates()

    df_up.to_csv(csv_dir / "deter_upregulated_degs.csv", index=False)
    df_down.to_csv(csv_dir / "deter_downregulated_degs.csv", index=False)

    print(f"[Deter et al. Pipeline] Completed. Identified DEGs across {len(config_deter['thresholds']['comparisons'])} condition pair at log2 ≥ {threshold}: {len(df_up)} up, {len(df_down)} down.")  

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_deter_pipeline(full_config)