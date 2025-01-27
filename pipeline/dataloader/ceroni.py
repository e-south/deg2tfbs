"""
--------------------------------------------------------------------------------
<deg2tfbs project>
ceroni.py

Module for reading and analyzing data described in Ceroni et al., where they 
combined RNA-seq with an in vivo assay to identify transcriptional changes that
occur in Escherichia coli when burdensome synthetic constructs are expressed.

"Burden-driven feedback control of gene expression"
DOI: 10.1038/nmeth.4635

Author(s): Eric J. South
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

from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_ceroni_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the Ceroni dataset using keys from the YAML config.
    
    default config_data looks like:
      {
        "dataset_key": "ceroni_burden",
        "sheet_name": "Mean FPKM",
        "usecols": "A:P",
        "header": 0,
        "skiprows": 3
      }
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data["sheet_name"],
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    
    # Tidy up DataFrame
    plasmid_names = ['pSB1C3-H3','pLys-M1','pSB1C3','pLys','pSB1C3-Lux','pPD864-LacZ','pPD864']
    time_points = ['15','60']
    headers = ['gene_id','gene_name'] + [f"{p}_{t}" for p in plasmid_names for t in time_points]
    df.columns = headers
    
    df_melted = pd.melt(
        df,
        id_vars=['gene_id', 'gene_name'],
        var_name='Plasmid_Time',
        value_name='DE_of_Mean_FPKM'
    )
    df_melted[['Plasmid','Time']] = df_melted['Plasmid_Time'].str.rsplit('_', n=1, expand=True)
    df_melted['Time'] = df_melted['Time'].astype(int)

    return df_melted


def ceroni_ma_plot(df, plasmid1, plasmid2, threshold=2.0, epsilon=1e-5, plot_path=None):
    """
    Creates and saves an MA-plot comparing plasmid1 and plasmid2 from the Ceroni dataset. 
    Here, plasmid1 and plasmid2 represent strain genotypes during an RNA-seq experiment, 
    enabling differential gene expression analysis relative to plasmid burden.
    """
    sns.set(style="ticks")
    df['log2_ratio'] = np.log2((df[plasmid1] + epsilon) / (df[plasmid2] + epsilon))
    colors = np.where(
        df['log2_ratio'] >= threshold, 'red',
        np.where(df['log2_ratio'] <= -threshold, 'green','lightgray')
    )

    plt.figure(figsize=(6,5))
    plt.scatter(df['average_expression'], df['log2_ratio'], c=colors, alpha=0.25, edgecolors='none')
    plt.xscale('log')
    plt.axhline(threshold, color='gray', linestyle='--')
    plt.axhline(-threshold, color='gray', linestyle='--')
    plt.title(f"MA Plot: {plasmid1} vs. {plasmid2}")
    plt.xlabel("Average Expression")
    plt.ylabel("Log2 Fold Change")
    sns.despine()

    if plot_path:
        plt.savefig(plot_path, dpi=150)
    plt.close()


def ceroni_thresholding(df, plasmid1, plasmid2, threshold=2.0, time_point=15, drop_zeros=False):
    """
    For a given plasmid pair, filter df for time_point,
    pivot, compute log2 ratio, do up/down selection.
    Returns upregulated, downregulated DataFrames.
    """
    df_sub = df[(df['Time'] == time_point) & (df['Plasmid'].isin([plasmid1, plasmid2]))].copy()
    pivoted = df_sub.pivot_table(
        index=['gene_id','gene_name'], 
        columns='Plasmid',
        values='DE_of_Mean_FPKM'
    )

    if drop_zeros:
        pivoted = pivoted[(pivoted[plasmid1]!=0) & (pivoted[plasmid2]!=0)]

    pivoted['log2_ratio'] = np.log2((pivoted[plasmid1]+1e-5)/(pivoted[plasmid2]+1e-5))
    pivoted['average_expression'] = (pivoted[plasmid1]+pivoted[plasmid2])/2
    pivoted.reset_index(inplace=True)

    # Figure out upregulated, downregulated
    up = pivoted[pivoted['log2_ratio'] >= threshold]
    down = pivoted[pivoted['log2_ratio'] <= -threshold]

    return up, down, pivoted


def run_ceroni_pipeline(full_config: dict):
    """
    full_config is the entire YAML content, containing:
      - output: { root_dir: ..., batch_identifier: ... }
      - ceroni: { data, thresholds, output subdirs, etc. }
    """

    # 1) Extract the "ceroni" config subset
    config_ceroni = full_config["ceroni"]

    # 2) Load data
    df_melted = read_ceroni_data(config_ceroni["data"])

    # 3) Build the final output paths
    project_root = Path(__file__).parent.parent.parent
    # top-level output folder, e.g. deg2tfbs/output
    output_root = project_root / full_config["output"]["root_dir"]

    # Append the batch identifier
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    # For ceroni, we define subfolders like "csv" and "plots" from config
    csv_dir = batch_dir / config_ceroni["output"]["csv_subdir"]
    plot_dir = batch_dir / config_ceroni["output"]["plot_subdir"]

    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # 4) Retrieve thresholds
    plasmid_pairs = config_ceroni["thresholds"]["plasmid_pairs"]
    threshold = config_ceroni["thresholds"]["log2_fc_threshold"]
    time_point = config_ceroni["thresholds"]["time_point"]
    drop_zeros = config_ceroni["thresholds"]["drop_zeros"]

    all_up = []
    all_down = []

    for (p1, p2) in plasmid_pairs:
        up, down, pivoted = ceroni_thresholding(
            df_melted,
            p1, p2,
            threshold=threshold,
            time_point=time_point,
            drop_zeros=drop_zeros
        )

        # Make a plot
        plot_path = plot_dir / f"ceroni_{p1}_vs_{p2}.png"
        ceroni_ma_plot(pivoted, p1, p2, threshold=threshold, plot_path=plot_path)

        # Save CSV
        up.to_csv(csv_dir / f"up_{p1}_vs_{p2}.csv", index=False)
        down.to_csv(csv_dir / f"down_{p1}_vs_{p2}.csv", index=False)

        all_up.append(up)
        all_down.append(down)

    up_all = pd.concat(all_up, ignore_index=True).drop_duplicates()
    down_all = pd.concat(all_down, ignore_index=True).drop_duplicates()

    up_all.to_csv(csv_dir / "DEGs_upregulated_all.csv", index=False)
    down_all.to_csv(csv_dir / "DEGs_downregulated_all.csv", index=False)

    print(f"[Ceroni Pipeline] Completed. Found: {len(up_all)} up, {len(down_all)} down total.")


if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_ceroni_pipeline(full_config)

    
    