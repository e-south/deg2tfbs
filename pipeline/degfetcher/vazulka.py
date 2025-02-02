"""
--------------------------------------------------------------------------------
<deg2tfbs project>
vazulka.py

Module for reading data from Vazulka et al., which performed RNA-seq to characterize
gene expression responses to Fab production in E. coli fed-batch processes.

The module identifies extreme up- and down-regulated genes using IQR-based outlier 
detection applied to the 'log2FoldChange' column.

"RNA-seq reveals multifaceted gene expression response to Fab production 
in Escherichia coli fed-batch processes with particular focus on ribosome stalling"
DOI: 10.1186/s12934-023-02278-w

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


def read_vazulka_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the Vazulka dataset using keys from the YAML config.

    E.g. config:
      {
        "dataset_key": "vazulka",
        "sheet_name": "Sheet1",
        "usecols": [...],
        "header": 0
      }
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def vazulka_distribution_plot(
    df: pd.DataFrame,
    log2_fc_col: str,
    up_idx,
    down_idx,
    output_path: Path
):
    """
    Plots a vertical distribution of log2 fold-change values with IQR outliers colored.
    
    Args:
        df: DataFrame containing the fold change data
        log2_fc_col: Column name with log2 fold change values
        up_idx: Pandas index of upregulated outliers
        down_idx: Pandas index of downregulated outliers
        output_path: Path to save the PNG
    """
    sns.set_style("ticks")
    plt.figure(figsize=(5,7))

    # Create color mapping
    color_map = ["gray"] * len(df)
    for i in up_idx:
        color_map[i] = "red"
    for i in down_idx:
        color_map[i] = "green"

    # Jitter x-values for visualization clarity
    xvals = np.random.uniform(-0.2, 0.2, size=len(df))
    
    plt.scatter(
        xvals,
        df[log2_fc_col],
        c=color_map,
        alpha=0.5,
        edgecolors="none"
    )
    plt.title("Vazulka et al. (TableS1)\nFab Overexpression (2 hours) versus control")
    plt.ylabel("log2(Fold Change)")
    plt.xticks([], [])

    sns.despine(top=True, right=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    

def run_vazulka_pipeline(full_config: dict):
    """
    Reads the dataset, 
    applies IQR-based outlier detection on 'log2FoldChange':
      up => log2FC > Q3 + 1.5*IQR
      down => log2FC < Q1 - 1.5*IQR

    Then saves CSV with columns ["gene","source","comparison"].
    """
    config_vazulka = full_config['datasets'].get("vazulka", None)
    if config_vazulka is None:
        print("[Vazulka Pipeline] No 'vazulka' config found. Skipping.")
        return

    # Load
    df = read_vazulka_data(config_vazulka["data"])
    assert "Gene" in df.columns, "Expected 'Gene' in Vazulka data"
    assert "log2FoldChange" in df.columns, "Expected 'log2FoldChange' in Vazulka data"

    # Build paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_vazulka["output"]["csv_subdir"]
    plot_dir = batch_dir / config_vazulka["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Aggregate duplicate genes by averaging log2FoldChange before IQR analysis
    df = df.groupby("Gene", as_index=False).agg({
        "log2FoldChange": "mean"
    })

    # Remove any remaining missing values
    df = df.dropna(subset=["log2FoldChange"])

    # IQR detection logic
    log2_fc = df["log2FoldChange"]
    q1, q3 = log2_fc.quantile([0.25, 0.75])
    iqr = q3 - q1
    
    up_bound = q3 + 1.5 * iqr
    down_bound = q1 - 1.5 * iqr

    up_df = df[log2_fc > up_bound]
    down_df = df[log2_fc < down_bound]

    # Generate plot
    vazulka_distribution_plot(
        df=df,
        log2_fc_col="log2FoldChange",
        up_idx=up_df.index,
        down_idx=down_df.index,
        output_path=plot_dir / "vazulka_fold_distribution_iqr.png"
    )

    # Build final data
    up_clean = pd.DataFrame({
        "gene": up_df["Gene"],
        "source": "vazulka",
        "comparison": "fab_production_2h_versus_control"
    }).drop_duplicates()

    down_clean = pd.DataFrame({
        "gene": down_df["Gene"],
        "source": "vazulka",
        "comparison": "fab_production_2h_versus_control"
    }).drop_duplicates()

    # Save
    up_csv = csv_dir / "vazulka_upregulated_degs.csv"
    down_csv = csv_dir / "vazulka_downregulated_degs.csv"
    up_clean.to_csv(up_csv, index=False)
    down_clean.to_csv(down_csv, index=False)

    print(f"[Vazulka et al. Pipeline] Completed. Identified DEGs across 1 condition pair with IQR-based outlier detection: {len(up_clean)} up, {len(down_clean)} down.")
    
if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_vazulka_pipeline(full_config)