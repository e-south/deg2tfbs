"""
--------------------------------------------------------------------------------
<deg2tfbs project>
jovanovic.py

Module for reading and analyzing data described in Jovanovic et al., 
where "Fold regulation" represents expression in MG1655 + PspA 
versus MG1655 + vector. A value of 1 => no change, >1 => upregulated, 
0 => not detected or zero expression. We ignore zero rows, 
then convert to log2 scale and identify outliers using the boxplot IQR method.

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


def read_jovanovic_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the Jovanovic dataset from an .xls file.
    Example config:
      {
        "dataset_key": "jovanovic",
        "sheet_name": "Sheet1",
        "usecols": ["Gene","Product","Fold regulation"],
        "header": 2
      }
    (Because columns start on row 3 => header=2 in Pandas)
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),  # e.g. "Sheet1"
        usecols=config_data.get("usecols", None),
        header=config_data.get("header", 0),       # e.g. 2
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
    Genes that are "extreme up" (index in up_idx) are colored red,
    while "extreme down" (index in down_idx) are colored green,
    and the rest are gray.

    Args:
        df (pd.DataFrame): Must contain a column [log2_fc_col].
        log2_fc_col (str): e.g. "log2_fold_reg"
        up_idx (Index): Pandas index of upregulated outliers
        down_idx (Index): Pandas index of downregulated outliers
        output_path (Path): Where to save the PNG.
    """
    sns.set_style("ticks")

    # Define colors
    color_map = ["gray"] * len(df)  # default
    for i in up_idx:
        color_map[i] = "red"
    for i in down_idx:
        color_map[i] = "green"

    # Make a vertical scatter
    plt.figure(figsize=(5,7))

    # We'll do a swarmplot or simple scatter with jitter in x
    # to show distribution along y-axis
    xvals = np.random.uniform(-0.2, 0.2, size=len(df))  # jitter
    plt.scatter(
        xvals,
        df[log2_fc_col],
        c=color_map,
        alpha=0.7,
        edgecolors="none"
    )

    plt.title("Jovanovic et al.\nDistribution of log2(Fold Regulation)")
    plt.ylabel("log2(Fold Regulation)")
    plt.xticks([], [])  # remove x-axis ticks

    sns.despine(left=True, bottom=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def run_jovanovic_pipeline(full_config: dict):
    """
    Orchestrates:
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

    # 1) Load data
    df = read_jovanovic_data(config_jov["data"])

    # 2) Prepare output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_jov["output"]["csv_subdir"]
    plot_dir = batch_dir / config_jov["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # 3) Retrieve column name
    fc_col = config_jov["thresholds"].get("fc_col", "Fold regulation")

    # 4) Ignore zero rows + transform to log2
    #    A value of 1 => no change, >1 => up, <1 => down
    #    skip zeros because log2(0) is undefined
    df = df[df[fc_col] > 0].copy()
    df["log2_fold_reg"] = np.log2(df[fc_col])

    # 5) Boxplot IQR approach for "extreme" up/down
    #    Q1, Q3 => IQR => define whiskers => outliers
    q1 = df["log2_fold_reg"].quantile(0.25)
    q3 = df["log2_fold_reg"].quantile(0.75)
    iqr = q3 - q1

    up_cutoff = q3 + 1.5 * iqr
    down_cutoff = q1 - 1.5 * iqr

    up_df = df[df["log2_fold_reg"] > up_cutoff].copy()
    down_df = df[df["log2_fold_reg"] < down_cutoff].copy()

    # 6) Save CSV for up/down
    up_out = csv_dir / "jovanovic_upregulated_iqr.csv"
    down_out = csv_dir / "jovanovic_downregulated_iqr.csv"
    up_df.to_csv(up_out, index=False)
    down_df.to_csv(down_out, index=False)

    print(f"[Jovanovic Pipeline] Up (IQR outliers) = {len(up_df)}, Down (IQR outliers) = {len(down_df)}")

    # 7) Make distribution plot (vertical)
    #   highlight outliers in red/green, rest in gray
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

    print("[Jovanovic Pipeline] Completed with IQR-based outlier detection.")
