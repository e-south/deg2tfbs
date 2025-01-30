"""
--------------------------------------------------------------------------------
<deg2tfbs project>
jovanovic.py

Module for reading and analyzing data from Jovanovic et al., which generated transcriptomic
readouts in *E. coli* strains with and without Psp-inducing protein IV secretin stress.

The module identifies extreme up- and down-regulated genes using IQR-based outlier 
detection applied to the log 'Fold regulation' distribution from a microarray dataset.

"Induction and Function of the Phage Shock Protein Extracytoplasmic Stress 
Response in Escherichia coli"
DOI: 10.1074/jbc.M602323200

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
from pathlib import Path

import yaml
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from deg2tfbs.pipeline.degfetcher.utils import load_dataset


def read_jovanovic_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the Jovanovic dataset using keys from the YAML config.

    E.g. config:
      {
        "dataset_key": "jovanovic",
        "sheet_name": "Sheet1",
        "usecols": ["Gene","Product","Fold regulation"],
        "header": 2
      }
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols", None),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def jovanovic_distribution_plot(
    df: pd.DataFrame,
    log2_fc_col: str,
    up_idx,
    down_idx,
    output_path: Path
):
    """
    Plots a vertical distribution of log2 fold-change values on the y-axis.

    Args:
        df (pd.DataFrame): Must contain a column [log2_fc_col].
        log2_fc_col (str): e.g. "log2_fold_reg"
        up_idx (Index): Pandas index of upregulated outliers
        down_idx (Index): Pandas index of downregulated outliers
        output_path (Path): Where to save the PNG.
    """
    sns.set_style("ticks")

    # Up- and down-regulated genes are colored red and green, respectively
    color_map = ["gray"] * len(df)  # default
    for i in up_idx:
        color_map[i] = "red"
    for i in down_idx:
        color_map[i] = "green"

    # Jitter the x-values for better visualization
    plt.figure(figsize=(5,7))

    xvals = np.random.uniform(-0.2, 0.2, size=len(df))  # jitter
    plt.scatter(
        xvals,
        df[log2_fc_col],
        c=color_map,
        alpha=0.5,
        edgecolors="none"
    )

    plt.title("Jovanovic et al.\nDistribution of log2(Fold Regulation)")
    plt.ylabel("log2(Fold Regulation)")
    plt.xticks([], [])  # remove x-axis ticks

    sns.despine(top=True, right=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def run_jovanovic_pipeline(full_config: dict):
    """
    1) Reading the Jovanovic dataset from .xls
    2) Ignoring zero rows in "Fold regulation"
    3) Converting to log2 scale
    4) Identifying outliers (extreme up/down) via boxplot IQR approach
    5) Plotting distribution on a vertical log2(FoldReg) axis
    6) Saving CSVs + figure

    YAML snippet example:
      jovanovic:
        data:
          dataset_key: "jovanovic"
          sheet_name: "Sheet1"
          usecols:
            - "Gene"
            - "Product"
            - "Fold regulation"
          header: 2
        thresholds:
          fc_col: "Fold regulation"
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """
    config_jov = full_config.get("jovanovic", None)
    if config_jov is None:
        print("[Jovanovic Pipeline] No 'jovanovic' config found. Skipping.")
        return

    # Load data
    df = read_jovanovic_data(config_jov["data"])

    # Prepare output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_jov["output"]["csv_subdir"]
    plot_dir = batch_dir / config_jov["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Retrieve column name
    fc_col = config_jov["thresholds"].get("fc_col", "Fold regulation")

    # Ignore zero rows + transform to log2
    #    A value of 1 => no change, >1 => up, <1 => down
    #    skip zeros because log2(0) is undefined
    df = df[df[fc_col] > 0].copy()
    df["log2_fold_reg"] = np.log2(df[fc_col])

    # Boxplot IQR approach for "extreme" up/down
    #    Q1, Q3 => IQR => define whiskers => outliers
    q1 = df["log2_fold_reg"].quantile(0.25)
    q3 = df["log2_fold_reg"].quantile(0.75)
    iqr = q3 - q1

    up_cutoff = q3 + 1.5 * iqr
    down_cutoff = q1 - 1.5 * iqr

    up_df = df[df["log2_fold_reg"] > up_cutoff].copy()
    down_df = df[df["log2_fold_reg"] < down_cutoff].copy()

    required_columns = ["gene", "source", "thresholds", "comparison"]
    
    # Define comparison
    target_condition = "IV_secretin_overproduction"
    reference_condition = "control"
    comparison_str = f"{target_condition}_versus_{reference_condition}"
    
    # Tidy up the DataFrames of up- and down-regulated genes
    up_clean = up_df.copy()
    up_clean["gene"] = up_clean["Gene"]
    up_clean["source"] = "jovanovic"
    up_clean["thresholds"] = "IQR"
    up_clean["comparison"] = comparison_str
    up_clean = up_clean[required_columns]

    down_clean = down_df.copy()
    down_clean["gene"] = down_clean["Gene"] 
    down_clean["source"] = "jovanovic"
    down_clean["thresholds"] = "IQR"
    down_clean["comparison"] = comparison_str
    down_clean = down_clean[required_columns]

    # Save CSV for up/down
    up_out = csv_dir / "jovanovic_upregulated_degs.csv"
    down_out = csv_dir / "jovanovic_downregulated_degs.csv"
    up_clean.to_csv(up_out, index=False)
    down_clean.to_csv(down_out, index=False)

    # Plot the distribution
    up_idx = up_df.index
    down_idx = down_df.index

    plot_fname = "jovanovic_fold_distribution_iqr.png"
    jovanovic_distribution_plot(
        df=df,
        log2_fc_col="log2_fold_reg",
        up_idx=up_idx,
        down_idx=down_idx,
        output_path=plot_dir / plot_fname
    )
    
    print(f"[Jovanovic et al. Pipeline] Completed. Identified DEGs across 1 condition pair at with IQR-based outlier detection: {len(up_df)} up, {len(down_df)} down.")
    
if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_jovanovic_pipeline(full_config)