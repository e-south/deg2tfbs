"""
--------------------------------------------------------------------------------
<deg2tfbs project>
mori.py

Module for loading and analyzing data described in Mori et al., which quantified 
absolute protein abundances in E. coli for > 2,000 proteins across > 60 growth 
conditions, including nutrient limitations, non-metabolic stresses, and non-planktonic 
states. All sample descriptions are listed in Datasets EV2 and EV3, while the absolute 
protein mass fractions are reported in Datasets EV8 and EV9.

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and saves an MA plot (average expression vs. log2 Fold Change).

"From coarse to fine: the absolute Escherichia coli proteome under diverse growth conditions"
DOI: 10.15252/msb.20209536

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
from pathlib import Path

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from deg2tfbs.pipeline.utils import load_dataset


def read_mori_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Mori dataset using keys from the YAML config.

      E.g. config:
      {
        "dataset_key": "mori_ev8",
        "sheet_name": "EV8-AbsoluteMassFractions-1",
        "usecols": ["Gene name","Lib-02","Lib-06","Lib-14","Lib-17","Lib-28","Lib-19","Lib-20"]
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


def compare_and_plot(
    df: pd.DataFrame,
    condition: str,
    reference: str,
    threshold: float,
    descriptors: dict = None,
    save_plots: bool = False,
    plot_dir: Path = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compare protein abundances for strains grown in a given 'condition' or
    'reference' media by computing log2 fold-change. Return two DataFrames 
    of upregulated and downregulated genes. Optionally save an MA-like plot.

    Args:
        df (pd.DataFrame): Data with columns matching 'condition' and 'reference'.
        condition (str): The column (condition) of interest.
        reference (str): The column used as reference.
        threshold (float): log2 fold-change cutoff.
        descriptors (dict): Optional mapping from column name to descriptive label.
        save_plots (bool): If True, saves a scatter plot.
        plot_dir (Path): Where to save the plot.

    Returns:
        up (pd.DataFrame), down (pd.DataFrame)
    """

    # Safely drop rows with zero or NaN in either column
    df_sub = df[["Gene name", condition, reference]].replace(0, np.nan).dropna().copy()

    # Compute log2 fold change and average
    df_sub["log2_fc"] = np.log2(df_sub[condition]) - np.log2(df_sub[reference])
    df_sub["avg_expr"] = (df_sub[condition] + df_sub[reference]) / 2

    # Up- and down-regulated genes are colored red and green, respectively
    df_sub["color"] = np.where(
        df_sub["log2_fc"] >= threshold, "red",
        np.where(df_sub["log2_fc"] <= -threshold, "green", "gray")
    )

    # Partition up- and downregulated
    up = df_sub[df_sub["log2_fc"] >= threshold].copy()
    down = df_sub[df_sub["log2_fc"] <= -threshold].copy()

    # Optionally plot
    if save_plots and plot_dir is not None:
        # Construct descriptive text if possible
        cond_label = descriptors.get(condition, condition) if descriptors else condition
        ref_label = descriptors.get(reference, reference) if descriptors else reference
        subtitle = f"log2({cond_label} / {ref_label})"

        plt.figure(figsize=(6,5))
        sns.set_style("ticks")

        plt.scatter(
            df_sub["avg_expr"],
            df_sub["log2_fc"],
            c=df_sub["color"],
            alpha=0.25,
            edgecolors="none"
        )
        plt.xscale("log")
        plt.axhline(threshold, color="gray", linestyle="--")
        plt.axhline(-threshold, color="gray", linestyle="--")
        
        plt.title(f"Mori et al.\n{subtitle}")
        plt.xlabel("Average Protein Abundance")
        plt.ylabel("Log2 Fold Change")

        sns.despine()
        plot_fname = f"mori_{condition}_vs_{reference}.png"
        plt.savefig(plot_dir / plot_fname, dpi=300)
        plt.close()

    return up, down


def run_mori_pipeline(full_config: dict):
    """
    Loads Mori dataset(s) and performs thresholding
    for each user-defined comparison.

    Example from config:
      mori:
        data:
          dataset_key: "mori_ev8"
          sheet_name: "EV8-AbsoluteMassFractions-1"
          usecols:
            - "Gene name"
            - "Lib-02"
            - "Lib-06"
            ...
        comparisons:
          - { condition: "Lib-19", reference: "Lib-20" }
          - { condition: "Lib-02", reference: "Lib-00-A1" }
          ...
        descriptors:
          Lib-00-A1: "Glucose minimal medium"
          Lib-02: "High osmolarity"
          ...
        thresholds:
          log2_fc_threshold: 2.0
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """
    config_mori = full_config['datasets']["mori"]

    # Load data
    df = read_mori_data(config_mori["data"])

    # Build output paths
    project_root = Path(__file__).parent.parent.parent
    
    # Output folder, e.g. deg2tfbs/pipeline/deg_fetcher
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_mori["output"]["csv_subdir"]
    plot_dir = batch_dir / config_mori["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Retrieve media conditions and thresholds
    descriptors = config_mori.get("descriptors", {})
    threshold = config_mori["thresholds"]["log2_fc_threshold"]

    # For each comparison, do the analysis
    all_up, all_down = [], []
    for comp in config_mori["comparisons"]:
        cond = comp["condition"]
        ref = comp["reference"]
        up, down = compare_and_plot(
            df=df, 
            condition=cond, 
            reference=ref, 
            threshold=threshold,
            descriptors=descriptors,
            save_plots=True,
            plot_dir=plot_dir
        )
        
        required_columns = ["gene", "source", "thresholds", "comparison"]
        
        # Tidy up the DataFrames of up- and down-regulated genes
        up_clean = up.copy()
        up_clean["gene"] = up_clean["Gene name"]
        up_clean["source"] = "mori"
        up_clean["thresholds"] = threshold
        up_clean["comparison"] = f"{comp["condition"]}_versus_{comp["reference"]}"
        up_clean = up_clean[required_columns]
        all_up.append(up_clean)

        down_clean = down.copy()
        down_clean["gene"] = down_clean["Gene name"] 
        down_clean["source"] = "mori"
        down_clean["thresholds"] = threshold
        down_clean["comparison"] = f"{comp["condition"]}_versus_{comp["reference"]}"
        down_clean = down_clean[required_columns]
        all_down.append(down_clean)

    # Combine all up-/down-regulated
    df_up = pd.concat(all_up, ignore_index=True).drop_duplicates()
    df_down = pd.concat(all_down, ignore_index=True).drop_duplicates()

    # Save CSV outputs
    df_up.to_csv(csv_dir / "mori_upregulated_degs.csv", index=False)
    df_down.to_csv(csv_dir / "mori_downregulated_degs.csv", index=False)

    print(f"[Mori et al. Pipeline] Completed. Identified DEGs across {len(config_mori['comparisons'])} condition pairs at log2 ≥ {threshold}: {len(df_up)} up, {len(df_down)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_mori_pipeline(full_config)

