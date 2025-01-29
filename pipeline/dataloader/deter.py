"""
--------------------------------------------------------------------------------
<deg2tfbs project>
deter.py

Module for reading and analyzing data described in Deter et al.
"Antibiotic tolerance is associated with a broad and complex 
transcriptional response in E. coli"

DOI: 10.1038/s41598-021-85509-7

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_deter_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the Deter dataset (CSV) according to keys from the YAML config.

    Example config_data might be:
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
    Creates and saves an MA-like plot comparing target vs. reference columns.

    Args:
      df_comp (pd.DataFrame): Must have 'log2_fc' and 'avg_expr' columns.
      target_avg_col (str): e.g. "3h amp media_avg"
      ref_avg_col (str): e.g. "3h fresh media_avg"
      threshold (float): log2 fold-change threshold
      plot_path (Path): Where to save the figure.
    """
    sns.set_style("ticks")

    # Color up/down
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

    # Create human-readable condition names:
    # E.g. "3h amp media_avg" => "3h amp media"
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
    Orchestrates reading Deter et al. data from CSV, 
    averaging replicates for multiple comparisons,
    computing log2 FC, generating up/down CSVs, 
    and saving an MA plot for each comparison.

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
        print("[Deter Pipeline] No 'deter' config found. Skipping.")
        return

    # 1) Load data
    df = read_deter_data(config_deter["data"])

    # 2) Prepare output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_deter["output"]["csv_subdir"]
    plot_dir = batch_dir / config_deter["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # 3) Retrieve threshold + comparisons
    threshold = config_deter["thresholds"].get("log2_fc_threshold", 2.0)
    comparisons = config_deter["thresholds"].get("comparisons", [])

    # 4) Helper to average replicates for a condition
    def average_condition(df_local, cond_prefix):
        # e.g. cond_prefix="3h amp media" => columns "3h amp media rep1", "3h amp media rep2", ...
        replicate_cols = [col for col in df_local.columns if col.startswith(cond_prefix)]
        if not replicate_cols:
            raise KeyError(f"No replicate columns found for '{cond_prefix}'")
        return df_local[replicate_cols].mean(axis=1)

    all_up = []
    all_down = []

    # 5) For each comparison in the config
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

        # 5.1) Compute log2 fold change
        eps = 1e-5
        fc_col = f"log2_fc_{target_cond}_vs_{ref_cond}"
        df[fc_col] = np.log2(
            (df[target_avg_col] + eps) / (df[ref_avg_col] + eps)
        )

        # 5.2) Average expression
        avg_expr_col = f"avg_expr_{target_cond}_vs_{ref_cond}"
        df[avg_expr_col] = (df[target_avg_col] + df[ref_avg_col]) / 2.0

        # 5.3) Mark up/down
        up = df[df[fc_col] >= threshold].copy()
        down = df[df[fc_col] <= -threshold].copy()

        # Collect for final summary
        all_up.append(up)
        all_down.append(down)

        # 5.4) Save partial CSV
        up_csv_name = f"deter_up_{target_cond}_vs_{ref_cond}.csv"
        down_csv_name = f"deter_down_{target_cond}_vs_{ref_cond}.csv"
        up.to_csv(csv_dir / up_csv_name, index=False)
        down.to_csv(csv_dir / down_csv_name, index=False)

        # 5.5) Create an MA plot
        # We'll create a small df with 'log2_fc' => df[fc_col], 'avg_expr' => df[avg_expr_col]
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

        print(f"[Deter Pipeline] Completed one comparison: {target_cond} vs {ref_cond}")
        print(f"   Up = {len(up)}, Down = {len(down)}")

    # 6) Combine all up/down across comparisons
    df_up = pd.concat(all_up, ignore_index=True).drop_duplicates()
    df_down = pd.concat(all_down, ignore_index=True).drop_duplicates()

    df_up.to_csv(csv_dir / "deter_DEGs_upregulated_all.csv", index=False)
    df_down.to_csv(csv_dir / "deter_DEGs_downregulated_all.csv", index=False)

    print(f"[Deter Pipeline] Summaries: Up = {len(df_up)} total, Down = {len(df_down)} total.")
