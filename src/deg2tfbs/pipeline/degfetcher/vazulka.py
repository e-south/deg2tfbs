"""
--------------------------------------------------------------------------------
<deg2tfbs project>
vazulka.py

Module for reading data from Vazulka et al., which performed RNA-seq to characterize
gene expression responses to Fab production in E. coli fed-batch processes.

The module identifies up- and down-regulated genes using a simplified filtering approach
applied to the 'log2FoldChange' column: positive values indicate up-regulation while 
negative values indicate down-regulation.

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
    Plots a vertical distribution of log₂ fold-change values with threshold-based 
    upregulated (red) and downregulated (green) genes highlighted.
    
    Args:
        df: DataFrame containing the fold change data.
        log2_fc_col: Column name with log₂ fold-change values.
        up_idx: Pandas index of genes with log₂FoldChange > 0.
        down_idx: Pandas index of genes with log₂FoldChange < 0.
        output_path: Path to save the PNG.
    """
    sns.set_style("ticks")
    plt.figure(figsize=(6, 5))

    # Create color mapping: default gray; mark upregulated red and downregulated green
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
    plt.title("Vazulka et al. (TableS1)\nlog₂(fab overexpression (2 h) / control; Filtered via positive/negative threshold")
    plt.ylabel("log₂(Fold Change)")
    plt.xlabel("Gene")
    plt.xticks([], [])

    sns.despine(top=True, right=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    

def run_vazulka_pipeline(full_config: dict):
    """
    Reads the dataset, aggregates duplicate genes, applies threshold-based filtering on 
    'log2FoldChange' (positive values as upregulated, negative values as downregulated),
    and then saves CSVs with columns ["gene", "source", "comparison"].
    """
    config_vazulka = full_config['datasets'].get("vazulka", None)
    if config_vazulka is None:
        print("[Vazulka Pipeline] No 'vazulka' config found. Skipping.")
        return

    # Load the Vazulka dataset
    df = read_vazulka_data(config_vazulka["data"])
    assert "Gene" in df.columns, "Expected 'Gene' in Vazulka data"
    assert "log2FoldChange" in df.columns, "Expected 'log2FoldChange' in Vazulka data"

    # Build file paths for outputs
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_vazulka["output"]["csv_subdir"]
    plot_dir = batch_dir / config_vazulka["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Aggregate duplicate genes by averaging log2FoldChange values
    df = df.groupby("Gene", as_index=False).agg({
        "log2FoldChange": "mean"
    })
    
    # Remove any remaining missing values and reset index
    df = df.dropna(subset=["log2FoldChange"]).reset_index(drop=True)

    # Threshold-based detection:
    # Upregulated if log2FoldChange > 0; Downregulated if log2FoldChange < 0.
    up_df = df[df["log2FoldChange"] > 0]
    down_df = df[df["log2FoldChange"] < 0]

    # Generate distribution plot with threshold-based outlier coloring
    vazulka_distribution_plot(
        df=df,
        log2_fc_col="log2FoldChange",
        up_idx=up_df.index,
        down_idx=down_df.index,
        output_path=plot_dir / "vazulka_fold_distribution_threshold.png"
    )

    # Build final data for CSV export
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

    # Save CSV outputs
    up_csv = csv_dir / "vazulka_upregulated_degs.csv"
    down_csv = csv_dir / "vazulka_downregulated_degs.csv"
    up_clean.to_csv(up_csv, index=False)
    down_clean.to_csv(down_csv, index=False)

    print(f"[Vazulka et al. Pipeline] Completed. Identified DEGs across 1 condition pair with threshold-based detection: {len(up_clean)} up, {len(down_clean)} down.")
    

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_vazulka_pipeline(full_config)
