"""
--------------------------------------------------------------------------------
<deg2tfbs project>
zhu.py

Module for reading and analyzing data from Zhu et al., where they compared abundances
of >2500 proteins in MG1655 during nutrient downshift (i.e., (p)ppGpp accumulation).
The proteomes of wild type cells at 0 min, 20 min, 40 min, and 80 min of postshift
were analyzed. The proteomes of relA-deficient strain at 0 min, 1.5 h, 3 h, and 4.5 h
of postshift were also analyzed. Supplementary Figure 5 contains iBAQ mass values.

"Stringent response ensures the timely adaptation of bacterial growth to nutrient 
downshift"
DOI: 10.1038/s41467-023-36254-0

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_zhu_data(config_data: dict) -> pd.DataFrame:
    """
    Reads Zhu dataset with keys from config.
    Typically:
      {
        "dataset_key": "zhu",
        "sheet_name": "supplementary data 2",
        "usecols": "B,C,I,K",
        "header": 0
      }
    We'll rename columns for clarity: [Protein names, Gene names, iBAQ mass E1, iBAQ mass E3].
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
    Classify up/down based on threshold, produce MA plot if requested.
    Returns up & down DataFrames.
    """
    df["average_iBAQ"] = (df["iBAQ mass E1"] + df["iBAQ mass E3"]) / 2
    df["log2_fc"] = np.log2(df["iBAQ mass E3"] / df["iBAQ mass E1"])
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
        plt.title("Zhu: Nutrient Downshift (RelA Overproduction vs WT)")
        plt.xlabel("Average iBAQ")
        plt.ylabel("Log2 FC")
        plt.savefig(plot_dir / "zhu_relA_overexpression_vs_WT.png", dpi=150)
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

    up, down = zhu_threshold_analysis(df, threshold=threshold, save_plots=True, plot_dir=plot_dir)
    up.to_csv(csv_dir / "DEGs_upregulated_zhu.csv", index=False)
    down.to_csv(csv_dir / "DEGs_downregulated_zhu.csv", index=False)

    print(f"[Zhu Pipeline] Completed. Found: {len(up)} up, {len(down)} down total.")
