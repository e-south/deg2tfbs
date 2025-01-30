"""
--------------------------------------------------------------------------------
<deg2tfbs project>
zhu.py

Module for loading and analyzing data from Zhu et al., which compared abundances
of >2500 proteins in MG1655 during nutrient downshift. The proteomes of wild type 
cells at 0 min, 20 min, 40 min, and 80 min of postshift were analyzed. The proteomes 
of relA-deficient strain at 0 min, 1.5 h, 3 h, and 4.5 h of postshift were also 
analyzed. Supplementary Figure 5 contains iBAQ mass values.

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and saves an MA plot (average expression vs. log2 Fold Change).

"Stringent response ensures the timely adaptation of bacterial growth to nutrient 
downshift"
DOI: 10.1038/s41467-023-36254-0

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_zhu_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Zhu dataset using keys from the YAML config.
    
    E.g. config:
      {
        "dataset_key": "zhu",
        "sheet_name": "supplementary data 2",
        "usecols": "B,C,I,K",
        "header": 0
      }
    Rename columns for clarity: [Protein names, Gene names, iBAQ mass E1, iBAQ mass E3].
    Then we drop zeros and compute log2_fc.
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    # Rename columns if needed:
    df.columns = ["Protein names", "Gene names", "iBAQ mass E1", "iBAQ mass E3"]
    df = df.replace(0, np.nan).dropna(subset=["iBAQ mass E1","iBAQ mass E3"])
    return df

def zhu_threshold_analysis(df, threshold=3.0, save_plots=False, plot_dir=None):
    """
    Calculate log2 fold change = log2(E3 / E1).
    Classify up/down based on threshold, produce MA plot.
    Returns up & down DataFrames.
    """
    df["average_iBAQ"] = (df["iBAQ mass E1"] + df["iBAQ mass E3"]) / 2
    df["log2_fc"] = np.log2(df["iBAQ mass E3"] / df["iBAQ mass E1"])
    
    # Up- and down-regulated genes are colored red and green, respectively
    df["color"] = np.where(df["log2_fc"] >= threshold, "red",
                  np.where(df["log2_fc"] <= -threshold, "green", "gray"))

    up = df[df["log2_fc"] >= threshold]
    down = df[df["log2_fc"] <= -threshold]

    if save_plots and plot_dir:
        sns.set_style("ticks")
        plt.figure(figsize=(6,5))
        plt.scatter(df["average_iBAQ"], df["log2_fc"], c=df["color"], alpha=0.25, edgecolors='none')
        plt.xscale('log')
        plt.axhline(threshold, color='gray', linestyle='--')
        plt.axhline(-threshold, color='gray', linestyle='--')
        plt.title("Zhu et al.\n log2 (RelA Overproduction / WT)")
        plt.xlabel("Average iBAQ")
        plt.ylabel("Log2 Fold Change")
        plt.savefig(plot_dir / "zhu_relA_overexpression_versus_WT.png", dpi=150)
        plt.close()

    return up, down

def run_zhu_pipeline(full_config: dict):
    """
    Example config in YAML:
      zhu:
        data:
          dataset_key: "zhu"
          sheet_name: "supplementary data 2"
          usecols: "B,C,I,K"
          header: 0
        thresholds:
          log2_fc_threshold: 3
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """
    config_zhu = full_config["zhu"]
    df = read_zhu_data(config_zhu["data"])

    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_zhu["output"]["csv_subdir"]
    plot_dir = batch_dir / config_zhu["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    threshold = config_zhu["thresholds"]["log2_fc_threshold"]

    # Perform threshold analysis 
    up, down = zhu_threshold_analysis(df, threshold=threshold, save_plots=True, plot_dir=plot_dir)
    
    required_columns = ["gene", "source", "thresholds", "comparison"]
    
    # Define comparison
    target_condition = "RelA_Overproduction"
    reference_condition = "WT"
    comparison_str = f"{target_condition}_vs_{reference_condition}"
    
    # Tidy up the DataFrames of up- and down-regulated genes
    up_clean = up.copy()
    up_clean["gene"] = up_clean["Gene names"]
    up_clean["source"] = "zhu"
    up_clean["thresholds"] = threshold
    up_clean["comparison"] = comparison_str
    up_clean = up_clean[required_columns]

    down_clean = down.copy()
    down_clean["gene"] = down_clean["Gene names"] 
    down_clean["source"] = "zhu"
    down_clean["thresholds"] = threshold
    down_clean["comparison"] = comparison_str
    down_clean = down_clean[required_columns]
    
    up_clean.to_csv(csv_dir / "zhu_upregulated_degs.csv", index=False)
    down_clean.to_csv(csv_dir / "zhu_downregulated_degs.csv", index=False)

    print(f"[Zhu et al. Pipeline] Completed. Identified DEGs across 1 condition pair at log2 ≥ {threshold}: {len(up)} up, {len(down)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_zhu_pipeline(full_config)