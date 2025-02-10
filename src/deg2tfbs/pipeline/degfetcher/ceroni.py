"""
--------------------------------------------------------------------------------
<deg2tfbs project>
ceroni.py

Module for loading and analyzing data described in Ceroni et al., which 
combined RNA-seq with in vivo assays to identify transcriptional changes that
occur in E. coli when burdensome synthetic constructs are expressed.

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and saves an MA plot (average expression vs. log2 Fold Change).

"Burden-driven feedback control of gene expression"
DOI: 10.1038/nmeth.4635

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


def read_ceroni_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Ceroni dataset using keys from the YAML config.
    
    E.g. config:
      {
        "dataset_key": "ceroni",
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
    
    # Specify which plasmids (i.e., strains) and time points we're interested in
    plasmid_names = ['pSB1C3-H3','pLys-M1','pSB1C3','pLys','pSB1C3-Lux','pPD864-LacZ','pPD864']
    time_points = ['15','60'] # minutes
    headers = ['gene_id','gene_name'] + [f"{p}_{t}" for p in plasmid_names for t in time_points]
    df.columns = headers
    
    # Tidy up the DataFrame
    df_melted = pd.melt(
        df,
        id_vars=['gene_id', 'gene_name'],
        var_name='Plasmid_Time',
        value_name='DE_of_Mean_FPKM'
    )
    df_melted[['Plasmid','Time']] = df_melted['Plasmid_Time'].str.rsplit('_', n=1, expand=True)
    df_melted['Time'] = df_melted['Time'].astype(int)

    return df_melted


def ceroni_ma_plot(
    df,
    plasmid1,
    plasmid2,
    threshold=2.0,
    epsilon=1e-5,
    plot_path=None,
    time_point=None
):
    """
    Creates and saves an MA-plot comparing strains carrying plasmid1 
    and plasmid2 from the Ceroni dataset.
    """
    # Set up the plot
    sns.set_theme(style="ticks")
    df['log2_ratio'] = np.log2((df[plasmid1] + epsilon) / (df[plasmid2] + epsilon))
    df['average_expression'] = (df[plasmid1] + df[plasmid2]) / 2

    # Up- and down-regulated genes are colored red and green, respectively
    colors = np.where(
        df['log2_ratio'] >= threshold, 'red',
        np.where(df['log2_ratio'] <= -threshold, 'green', 'lightgray')
    )
    # Plot 
    plt.figure(figsize=(6,5))
    plt.scatter(df['average_expression'], df['log2_ratio'], c=colors, alpha=0.25, edgecolors='none')
    plt.xscale('log')
    plt.axhline(threshold, color='gray', linestyle='--')
    plt.axhline(-threshold, color='gray', linestyle='--')

    # Title and subtitle
    if time_point is not None:
        plt.title(f"Ceroni et al.\nlog2({plasmid1} / {plasmid2}); {time_point} mins")
    else:
        plt.title(f"Ceroni et al.\nlog2({plasmid1} / {plasmid2})")

    plt.xlabel("Average Expression")
    plt.ylabel("Log2 Fold Change")
    sns.despine()

    if plot_path:
        plt.savefig(plot_path, dpi=300)
    plt.close()


def ceroni_thresholding(df, plasmid1, plasmid2, threshold=2.0, time_point=15, drop_zeros=False):
    """
    For a given pair of strains carrying different plasmids, filter df by time_point, 
    pivot the data, compute the log₂ ratio between average FPKM values, and classify 
    genes as upregulated or downregulated.

    Returns two DataFrames: upregulated and downregulated.
    """
    df_sub = df[(df['Time'] == time_point) & (df['Plasmid'].isin([plasmid1, plasmid2]))].copy()
    pivoted = df_sub.pivot_table(
        index=['gene_id','gene_name'], 
        columns='Plasmid',
        values='DE_of_Mean_FPKM'
    )

    # Drop rows with zeros in either plasmid
    if drop_zeros:
        pivoted = pivoted[(pivoted[plasmid1]!=0) & (pivoted[plasmid2]!=0)]

    # Compute log2 ratio and average expression
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

    # Extract the "ceroni" config subset
    config_ceroni = full_config["datasets"]["ceroni"]

    # Load data
    df_melted = read_ceroni_data(config_ceroni["data"])

    # Build the final output paths
    project_root = Path(__file__).parent.parent.parent
    
    # Output folder, e.g. deg2tfbs/pipeline/deg_fetcher
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']

    # Append the batch identifier
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    # Define subfolders for "csv" and "plots" values from config
    csv_dir = batch_dir / config_ceroni["output"]["csv_subdir"]
    plot_dir = batch_dir / config_ceroni["output"]["plot_subdir"]

    # Create the directories if they don't exist
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Retrieve thresholds
    plasmid_pairs = config_ceroni["thresholds"]["plasmid_pairs"]
    threshold = config_ceroni["thresholds"]["log2_fc_threshold"]
    time_point = config_ceroni["thresholds"]["time_point"]
    drop_zeros = config_ceroni["thresholds"]["drop_zeros"]
    time_point = config_ceroni["thresholds"]["time_point"]  # e.g. 15, 60

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

        required_columns = ["gene", "source", "thresholds", "comparison"]
        
        # Tidy up the DataFrames of up- and down-regulated genes
        up_clean = up.copy()
        up_clean["gene"] = up_clean["gene_name"]
        up_clean["source"] = "ceroni"
        up_clean["thresholds"] = threshold
        up_clean["comparison"] = f"{p1}_versus_{p2}"
        up_clean = up_clean[required_columns]
        all_up.append(up_clean)

        down_clean = down.copy()
        down_clean["gene"] = down_clean["gene_name"] 
        down_clean["source"] = "ceroni"
        down_clean["thresholds"] = threshold
        down_clean["comparison"] = f"{p1}_versus_{p2}"
        down_clean = down_clean[required_columns]
        all_down.append(down_clean)

        # Save MA plot: fold changes (M values) against average expression (A values)
        plot_path = plot_dir / f"ceroni_{p1}_versus_{p2}.png"
        ceroni_ma_plot(pivoted, p1, p2, threshold=threshold, plot_path=plot_path)

    # Combine all up and down DEGs
    up_all = pd.concat(all_up, ignore_index=True).drop_duplicates()
    down_all = pd.concat(all_down, ignore_index=True).drop_duplicates()

    up_all.to_csv(csv_dir / "ceroni_upregulated_degs.csv", index=False)
    down_all.to_csv(csv_dir / "ceroni_downregulated_degs.csv", index=False)

    print(f"[Ceroni et al. Pipeline] Completed. Identified DEGs across {len(plasmid_pairs)} condition pairs at log2 ≥ {threshold}: {len(up_all)} up, {len(down_all)} down.")  

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_ceroni_pipeline(full_config)
