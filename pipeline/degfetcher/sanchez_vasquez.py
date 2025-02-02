"""
--------------------------------------------------------------------------------
<deg2tfbs project>
sanchez_vasquez.py

Module for reading Sanchez‐Vazquez et al., which performed RNA‐seq of E. coli 
with and without ppGpp‐binding sites on RNAP, characterizing gene expression 
changes at 5–10 min.

The module identifies extreme up- and down-regulated genes using IQR-based outlier 
detection applied to the '1+2+ 5 min' column from the RNA-deg dataset.

"Genome-wide effects on Escherichia coli transcription from ppGpp binding to 
its two sites on RNA polymerase"
DOI: 10.1073/pnas.1819682116

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
from pathlib import Path

import yaml
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from deg2tfbs.pipeline.utils import load_dataset


def read_sanchez_vasquez_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Sanchez Vasquez dataset using keys from the YAML config.
    
    Expected columns: 'Gene', '1+2+ 5 min Category', '1+2+ 10 min Category', etc.
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def sanchez_distribution_plot(
    df: pd.DataFrame,
    log2_fc_col: str,
    up_mask,  # Change from up_idx to up_mask
    down_mask,  # Change from down_idx to down_mask
    output_path: Path
):
    sns.set_style("ticks")
    plt.figure(figsize=(5,7))

    # Create color map using direct boolean masks
    color_map = np.select(
        [up_mask, down_mask],
        ["red", "green"],
        default="gray"
    )

    xvals = np.random.uniform(-0.2, 0.2, size=len(df))
    plt.scatter(xvals, df[log2_fc_col], c=color_map, alpha=0.5, linewidth=0.5)

    plt.title("Sanchez-Vazquez et al.\nLog2 Fold Change (5 min) with IQR Filtering")
    plt.ylabel("log2(Fold Change)")
    plt.xticks([], [])
    
    sns.despine(top=True, right=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()



def run_sanchez_vasquez_pipeline(full_config: dict):
    """
    Reads the data, 
    Up => '1+2+ 5 min Category' == 'B'
    Down => '1+2+ 5 min Category' == 'A'

    comparison => "relA_overexpression_versus_control"
    source => "sanchez_vasquez"
    """
    config_sv = full_config['datasets'].get("sanchez_vasquez", None)
    if config_sv is None:
        print("[Sanchez-Vazquez Pipeline] No 'sanchez_vasquez' config found. Skipping.")
        return

    # Load
    df = read_sanchez_vasquez_data(config_sv["data"])

    assert "Gene" in df.columns, "Expected 'Gene' in Sanchez-Vazquez data"
    cat_col = "1+2+ 5 min Category"
    assert cat_col in df.columns, f"Expected '{cat_col}' in Sanchez-Vazquez data"

    # Build output path
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_sv["output"]["csv_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = batch_dir / config_sv["output"]["plot_subdir"]
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Plot distribution
    log2_fc_col = "1+2+ 5 min"  # Verify column name matches your data
    df[log2_fc_col] = pd.to_numeric(df[log2_fc_col], errors='coerce')
    df = df.dropna(subset=[log2_fc_col, cat_col])
    
    # Calculate IQR bounds for outlier detection
    q1, q3 = df[log2_fc_col].quantile([0.25, 0.75])
    iqr = q3 - q1
    upper_bound = q3 + 1.5*iqr
    lower_bound = q1 - 1.5*iqr
    
    up_filter = (df[log2_fc_col] > upper_bound)
    down_filter = (df[log2_fc_col] < lower_bound)
    df_up = df[up_filter].copy()
    df_down = df[down_filter].copy()

    # Plot
    plot_path = plot_dir / "sanchez_vasquez_iqr.png"
    sanchez_distribution_plot(
        df=df,
        log2_fc_col=log2_fc_col,
        up_mask=up_filter,  # Pass the boolean mask directly
        down_mask=down_filter,  # Pass the boolean mask directly
        output_path=plot_path
    )
    
    # Validation asserts
    assert up_filter.sum() > 0, "No upregulated genes passing both IQR and category filters"
    assert down_filter.sum() > 0, "No downregulated genes passing both filters"
    
    up_clean = pd.DataFrame({
        "gene": df_up["Gene"],
        "source": "sanchez_vasquez",
        "comparison": "relA_overexpression_versus_control"
    }).drop_duplicates()

    down_clean = pd.DataFrame({
        "gene": df_down["Gene"],
        "source": "sanchez_vasquez",
        "comparison": "relA_overexpression_versus_control"
    }).drop_duplicates()

    # Save
    up_csv = csv_dir / "sanchez_vasquez_upregulated_degs.csv"
    down_csv = csv_dir / "sanchez_vasquez_downregulated_degs.csv"
    up_clean.to_csv(up_csv, index=False)
    down_clean.to_csv(down_csv, index=False)

    print(f"[Sanchez-Vazquez et al. Pipeline] Completed. Identified DEGs across 1 condition pair with IQR-based outlier detection: {len(up_clean)} up, {len(down_clean)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_sanchez_vasquez_pipeline(full_config)