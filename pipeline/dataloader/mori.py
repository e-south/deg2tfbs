"""
--------------------------------------------------------------------------------
<deg2tfbs project>
mori.py

Module for reading and analyzing data described in Mori et al., which accurately 
quantified absolute protein abundances in Escherichia coli for > 2,000 proteins
over > 60 growth conditions, including nutrient limitations, non-metabolic stresses,
and non-planktonic states. All sample descriptions are listed in Datasets
EV2 and EV3, while the absolute protein mass fractions are reported in Datasets EV8 and EV9.

"From coarse to fine: the absolute Escherichia coli proteome under diverse 
growth conditions"
DOI: 10.15252/msb.20209536

 EV8 dataset; columns:
         Lib-02, High osmolarity
         Lib-06, poor carbon (acetate)
         Lib-14, Oxidative stress
         Lib-17, Low pH (acetate)
         Lib-28, Xylose + oxaloacetate
         Lib-00-A1, Glucose minimal medium <- reference condition
         
         Lib-19, LB stationary
         Lib-20, LB log-phase <- reference condition
    
EV9 dataset; columns:
    F8, C-limitation, 0.34 (growth rate 1/h)
    C2, C-limitation, 0.91 (growth rate 1/h)


Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_mori_data(config_data: dict) -> pd.DataFrame:
    """
    Reads Mori dataset(s) using keys from the YAML config.

    Example config_data for EV8:
      {
        "dataset_key": "mori_ev8",
        "sheet_name": "EV8-AbsoluteMassFractions-1",
        "usecols": ["Gene name","Lib-02","Lib-06","Lib-14","Lib-17","Lib-28","Lib-19","Lib-20"]
      }

    If you have a second dataset (e.g., EV9), you might define a second config block 
    or read them in the same function. For simplicity, we handle just one dataset_key 
    or parse logic as needed.
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df

def mori_threshold_analysis(df, conditions, reference_col, threshold=2.0, save_plots=False, plot_dir=None):
    """
    For each column in 'conditions', compare to 'reference_col', compute log2 FC,
    and produce an MA-like plot. Returns an aggregated up/down list of genes.
    """

    all_up = []
    all_down = []

    # Avoid zero in reference
    df = df.replace(0, np.nan)
    for cond in conditions:
        if cond == reference_col:
            continue

        df_sub = df[[reference_col, cond]].dropna().copy()
        df_sub["log2_fc"] = np.log2(df_sub[cond]) - np.log2(df_sub[reference_col])
        df_sub["avg_expr"] = (df_sub[cond] + df_sub[reference_col]) / 2
        df_sub["color"] = np.where(df_sub["log2_fc"] >= threshold, "red",
                          np.where(df_sub["log2_fc"] <= -threshold, "green", "gray"))

        up = df_sub[df_sub["log2_fc"] >= threshold]
        down = df_sub[df_sub["log2_fc"] <= -threshold]
        all_up.append(up)
        all_down.append(down)

        # Optional plot
        if save_plots and plot_dir is not None:
            sns.set_style("ticks")
            plt.figure(figsize=(6,5))
            plt.scatter(df_sub["avg_expr"], df_sub["log2_fc"],
                        c=df_sub["color"], alpha=0.25, edgecolors='none')
            plt.xscale('log')
            plt.axhline(threshold, color='gray', linestyle='--')
            plt.axhline(-threshold, color='gray', linestyle='--')
            plt.title(f"Mori comparison: {cond} vs. {reference_col}")
            plt.xlabel("Average Abundance")
            plt.ylabel("Log2 FC")
            sns.despine()
            plot_path = plot_dir / f"mori_{cond}_vs_{reference_col}.png"
            plt.savefig(plot_path, dpi=150)
            plt.close()

    # Combine
    df_up = pd.concat(all_up, ignore_index=True).drop_duplicates()
    df_down = pd.concat(all_down, ignore_index=True).drop_duplicates()

    return df_up, df_down

def run_mori_pipeline(full_config: dict):
    """
    Orchestrates loading the Mori dataset(s) and 
    performing a threshold-based analysis.

    Example usage from config:
      mori:
        data:
          dataset_key: "mori_ev8"
          sheet_name: "EV8-AbsoluteMassFractions-1"
          usecols:
            - "Gene name"
            - "Lib-02"
            - "Lib-06"
            ...
        thresholds:
          reference_col: "Lib-20"
          conditions:
            - "Lib-02"
            - "Lib-06"
            ...
          log2_fc_threshold: 2.0
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """

    config_mori = full_config["mori"]
    df = read_mori_data(config_mori["data"])

    # Build output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_mori["output"]["csv_subdir"]
    plot_dir = batch_dir / config_mori["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Get thresholds
    threshold = config_mori["thresholds"]["log2_fc_threshold"]
    conditions = config_mori["thresholds"]["conditions"]
    reference_col = config_mori["thresholds"]["reference_col"]

    up_all, down_all = mori_threshold_analysis(
        df, 
        conditions=conditions,
        reference_col=reference_col,
        threshold=threshold,
        save_plots=True,
        plot_dir=plot_dir
    )

    up_all.to_csv(csv_dir / "DEGs_upregulated_all_mori.csv", index=False)
    down_all.to_csv(csv_dir / "DEGs_downregulated_all_mori.csv", index=False)

    print(f"[Mori Pipeline] Completed. Found: {len(up_all)} up, {len(down_all)} down total.")
