"""
--------------------------------------------------------------------------------
<deg2tfbs project>
kim.py

Module for reading Kim et al. dataset, which applied RNA‐seq to capture early,
middle, and late stages of heat stress (2 min–40 h) in E. coli. 

The module identifies extreme up- and down-regulated genes using IQR-based outlier 
detection applied to the '1h' column, which is in reference to the -30min column.

"Heat-responsive and time-resolved transcriptome and metabolome analyses of 
Escherichia coli uncover thermo-tolerant mechanisms"
DOI: 10.1038/s41598-020-74606-8

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

from deg2tfbs.pipeline.utils import load_dataset


def read_kim_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Kim dataset using keys from the YAML config.
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def kim_distribution_plot(
    df: pd.DataFrame,
    log2_fc_col: str,
    up_idx,
    down_idx,
    output_path: Path
):
    """
    Creates a jittered distribution plot of log2 fold changes with IQR outliers colored.
    """
    sns.set_style("ticks")
    plt.figure(figsize=(6,5))


    # Create color mapping
    color_map = ["gray"] * len(df)
    for i in up_idx:
        color_map[i] = "red"
    for i in down_idx:
        color_map[i] = "green"

    xvals = np.random.uniform(-0.2, 0.2, size=len(df))
    plt.scatter(xvals, df[log2_fc_col], c=color_map, alpha=0.5, edgecolors="none")

    plt.title("Kim et al.\n log2(42°C Heat Stress / Control (1h)); Filtered via +- 1.5*IQR")
    plt.ylabel("log2 Fold Change")
    plt.xlabel("Gene")
    plt.xticks([], [])

    sns.despine(top=True, right=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def run_kim_pipeline(full_config: dict):
    """
    Reads Kim data, applies IQR-based outlier detection on '1h' log2FC:
      up => > Q3 + 1.5*IQR
      down => < Q1 - 1.5*IQR
    """
    config_kim = full_config['datasets'].get("kim", None)
    if config_kim is None:
        print("[Kim Pipeline] No 'kim' config found. Skipping.")
        return

    # Load and clean data
    df = read_kim_data(config_kim["data"])
    assert "Gene" in df.columns, "Kim data must have a 'Gene' column"
    assert "1h" in df.columns, "Kim data must have a '1h' column for log2 FC"
    
    # Handle duplicates by averaging
    df = df.groupby("Gene", as_index=False).agg({"1h": "mean"}).dropna()

    # Build paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_kim["output"]["csv_subdir"]
    plot_dir = batch_dir / config_kim["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # IQR detection logic
    log2_fc = df["1h"]
    q1, q3 = log2_fc.quantile([0.25, 0.75])
    iqr = q3 - q1
    
    up_bound = q3 + 1.5 * iqr
    down_bound = q1 - 1.5 * iqr

    up_df = df[log2_fc > up_bound]
    down_df = df[log2_fc < down_bound]

    # Generate plot
    kim_distribution_plot(
        df=df,
        log2_fc_col="1h",
        up_idx=up_df.index,
        down_idx=down_df.index,
        output_path=plot_dir / "kim_fold_distribution_iqr.png"
    )

    # Build final data
    up_clean = pd.DataFrame({
        "gene": up_df["Gene"],
        "source": "kim",
        "comparison": "42C_versus_control"
    }).drop_duplicates()

    down_clean = pd.DataFrame({
        "gene": down_df["Gene"],
        "source": "kim",
        "comparison": "42C_versus_control"
    }).drop_duplicates()

    # Save CSVs
    up_csv = csv_dir / "kim_upregulated_degs.csv"
    down_csv = csv_dir / "kim_downregulated_degs.csv"
    up_clean.to_csv(up_csv, index=False)
    down_clean.to_csv(down_csv, index=False)

    print(f"[Kim et al. Pipeline] Completed. Identified DEGs across 1 condition pair with IQR-based outlier detection: {len(up_clean)} up, {len(down_clean)} down.")  

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_kim_pipeline(full_config)