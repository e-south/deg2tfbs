"""
--------------------------------------------------------------------------------
<deg2tfbs project>
rajacharya.py

Module for reading and analyzing data described in Rajacharya et al., which 
investigated parent vs. recombinant strains (induced at various time points) via 
proteomics to track expression burden.

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and plots a volcano figure (MinusLog10PAd vs. log2 Fold Change).

"Proteomics and metabolic burden analysis to understand the impact of recombinant 
protein production in E. coli"
DOI: 10.1038/s41598-024-63148-y

Module Author(s): Eric J. South
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


def read_rajacharya_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the Rajacharya dataset using keys from the YAML config.

    E.g. config:
      {
        "dataset_key": "rajacharya",
        "sheet_name": "Sheet 1",
        "usecols": ["Accession","genes","logFC","PValue"],
        "header": 3
      }
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )

    # Rename columns to something consistent
    df.rename(columns={
        "Accession": "Accession",
        "genes": "Gene",
        "logFC": "Log2FoldChange",
        "PValue": "PValue"
    }, inplace=True)

    # Add a column for -log10(PValue)
    df["MinusLog10PValue"] = -np.log10(df["PValue"].replace(0, np.nan))  # avoid log10(0)

    return df


def rajacharya_volcano_plot(df, log2_fc_col, minus_log10p_col, threshold, output_path):
    """
    Creates and saves a volcano plot:
      x-axis: Log2FoldChange
      y-axis: -log10(PValue)
    """
    # Up- and down-regulated genes are colored red and green, respectively
    df["color"] = np.where(
        df[log2_fc_col] >= threshold, "red",
        np.where(df[log2_fc_col] <= -threshold, "green", "gray")
    )

    plt.figure(figsize=(7,6))
    sns.set_style("ticks")

    # Plot the volcano
    plt.scatter(
        df[log2_fc_col],
        df[minus_log10p_col],
        c=df["color"],
        alpha=0.25,
        edgecolors="none"
    )
    # Threshold lines
    plt.axvline(x=threshold, color='gray', linestyle='--')
    plt.axvline(x=-threshold, color='gray', linestyle='--')

    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-log10(PValue)")
    plt.title("Rajacharya et al.\nVolcano Plot")
    sns.despine()

    plt.savefig(output_path, dpi=150)
    plt.close()


def run_rajacharya_pipeline(full_config: dict):
    """
    Orchestrates:
      1) Reading the Rajacharya dataset from config
      2) Plotting a volcano plot (log2 FC vs. -log10 PValue)
      3) Identifying up/down regulated genes at a user-defined threshold
      4) Saving the up/down CSVs and the volcano PNG

    Example config snippet:
      rajacharya:
        data:
          dataset_key: "rajacharya"
          sheet_name: "Sheet 1"
          usecols: ["Accession","genes","logFC","PValue"]
          header: 3
        thresholds:
          log2_fc_threshold: 2.0
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """
    config_raj = full_config.get("rajacharya", None)
    if config_raj is None:
        print("[Rajacharya Pipeline] No 'rajacharya' config found. Skipping.")
        return

    # Load data
    df = read_rajacharya_data(config_raj["data"])

    # Build output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_raj["output"]["csv_subdir"]
    plot_dir = batch_dir / config_raj["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Retrieve threshold
    threshold = config_raj["thresholds"].get("log2_fc_threshold", 2.0)

    # Make volcano plot
    log2_fc_col = "Log2FoldChange"
    minus_log10p_col = "MinusLog10PValue"
    plot_fname = "rajacharya_volcano.png"
    rajacharya_volcano_plot(
        df=df,
        log2_fc_col=log2_fc_col,
        minus_log10p_col=minus_log10p_col,
        threshold=threshold,
        output_path=plot_dir / plot_fname
    )

    # Filter up/down
    up = df[df[log2_fc_col] >= threshold].copy()
    down = df[df[log2_fc_col] <= -threshold].copy()

    required_columns = ["gene", "source", "thresholds", "comparison"]
    
    # Define comparison
    target_condition = "acyl-ACP_reductase_overproduction"
    reference_condition = "control"
    comparison_str = f"{target_condition}_vs_{reference_condition}"
    
    # Tidy up the DataFrames of up- and down-regulated genes
    up_clean = up.copy()
    up_clean["gene"] = up_clean["Gene"]
    up_clean["source"] = "rajacharya"
    up_clean["thresholds"] = threshold
    up_clean["comparison"] = comparison_str
    up_clean = up_clean[required_columns]

    down_clean = down.copy()
    down_clean["gene"] = down_clean["Gene"] 
    down_clean["source"] = "rajacharya"
    down_clean["thresholds"] = threshold
    down_clean["comparison"] = comparison_str
    down_clean = down_clean[required_columns]

    # Save CSV
    up_csv = csv_dir / "rajacharya_upregulated.csv"
    down_csv = csv_dir / "rajacharya_downregulated.csv"
    up_clean.to_csv(up_csv, index=False)
    down_clean.to_csv(down_csv, index=False)

    print(f"[Rajacharya et al. Pipeline] Completed. Identified DEGs across 1 condition pair at log2 ≥ {threshold}: {len(up)} up, {len(down)} down.")  

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_rajacharya_pipeline(full_config)