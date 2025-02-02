"""
--------------------------------------------------------------------------------
<deg2tfbs project>
bie.py

Module for loading and analyzing data described in Bie et al. which performed 
transcriptomic analysis on how E. coli responds to nine representative 
classes of antibiotics (tetracycline, mitomycin C, imipenem, ceftazidime, kanamycin, 
ciprofloxacin, polymyxin E, erythromycin, and chloramphenicol).

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and saves an MA plot (average expression vs. log2 Fold Change).

"Comparative Analysis of Transcriptomic Response of Escherichia coli 
K-12 MG1655 to Nine Representative Classes of Antibiotics"
DOI: 10.1128/spectrum.00317-23

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


def read_bie_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the Bie dataset using keys from the YAML config.

    E.g. config:
      {
        "dataset_key": "bie",
        "sheet_name": "Sheet1",
        "usecols": [...],
        "header": 1
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


def bie_ma_plot(
    df_cond: pd.DataFrame,
    readcount_abx_col: str,
    readcount_h2o_col: str,
    log2_fc_col: str,
    threshold: float,
    antibiotic_name: str,
    output_path: Path
):
    """
    Creates and saves an MA-like plot for a given antibiotic vs. H2O.
    x-axis: average expression = (abx_readcount + h2o_readcount)/2
    y-axis: log2_fold_change from the dataset.
    """
    sns.set_style("ticks")

    # Compute average expression for the MA plot
    df_cond["avg_expr"] = (df_cond[readcount_abx_col] + df_cond[readcount_h2o_col]) / 2.0

    # Up- and down-regulated genes are colored red and green, respectively
    df_cond["color"] = np.where(
        df_cond[log2_fc_col] >= threshold, "red",
        np.where(df_cond[log2_fc_col] <= -threshold, "green", "gray")
    )

    # Plot the MA-like figure
    plt.figure(figsize=(6,5))
    plt.scatter(
        df_cond["avg_expr"],
        df_cond[log2_fc_col],
        c=df_cond["color"],
        alpha=0.25,
        edgecolors="none"
    )
    plt.xscale("log")
    plt.axhline(y=threshold, color="gray", linestyle="--")
    plt.axhline(y=-threshold, color="gray", linestyle="--")

    plt.title(f"Bie et al.\n log2 ({antibiotic_name} / Water)")
    plt.xlabel("Average Expression (read count)")
    plt.ylabel("log2 Fold Change")
    sns.despine()

    plt.savefig(output_path, dpi=150)
    plt.close()


def run_bie_pipeline(full_config: dict):
    """
      1) Reads the dataset from config
      2) For each antibiotic block in thresholds["conditions"], we:
         - Identify relevant columns (readcount_abx, readcount_h2o, log2_fc)
         - Generate an MA plot (avg_expr vs. log2_fc)
         - Determine up/down by threshold
         - Save partial CSVs for up/down
      3) Combine up/down across all antibiotics
      4) Save final aggregated CSVs
    """
    config_bie = full_config['datasets'].get("bie", None)
    if config_bie is None:
        print("[Bie Pipeline] No 'bie' config found. Skipping.")
        return

    # Load data
    df = read_bie_data(config_bie["data"])

    # Build output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_bie["output"]["csv_subdir"]
    plot_dir = batch_dir / config_bie["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Retrieve threshold + conditions
    threshold = config_bie["thresholds"].get("log2_fc_threshold", 2.0)
    conditions = config_bie["thresholds"].get("conditions", [])

    all_up = []
    all_down = []

    # For each antibiotic definition in the config
    for cond_info in conditions:
        antibiotic_name    = cond_info["antibiotic_name"]
        gene_col           = cond_info["gene_col"]
        abx_readcount_col  = cond_info["abx_readcount_col"]
        h2o_readcount_col  = cond_info["h2o_readcount_col"]
        log2_fc_col        = cond_info["log2_fc_col"]
        sig_col            = cond_info.get("sig_col", None)  # optional significance

        # Subset DataFrame to relevant columns
        req_cols = [gene_col, abx_readcount_col, h2o_readcount_col, log2_fc_col]
        if sig_col:
            req_cols.append(sig_col)
        
        # Check for missing columns
        missing = [c for c in req_cols if c not in df.columns]
        if missing:
            print(f"[Bie Pipeline] Missing columns {missing} for antibiotic {antibiotic_name}, skipping.")
            continue

        # Create a small sub-DataFrame
        df_cond = df[req_cols].dropna().copy()

        # Create an MA-like plot
        plot_fname = f"bie_ma_{antibiotic_name}_vs_H2O.png"
        plot_path  = plot_dir / plot_fname
        bie_ma_plot(
            df_cond,
            readcount_abx_col=abx_readcount_col,
            readcount_h2o_col=h2o_readcount_col,
            log2_fc_col=log2_fc_col,
            threshold=threshold,
            antibiotic_name=antibiotic_name,
            output_path=plot_path
        )

        # Identify up/down
        up = df_cond[df_cond[log2_fc_col] >= threshold].copy()
        down = df_cond[df_cond[log2_fc_col] <= -threshold].copy()

        # Generate comparison string (KANvsH2O, CIPvsH2O, etc)
        comparison_str = f"{antibiotic_name}vsH2O"

        required_columns = ["gene", "source", "thresholds", "comparison"]
        
        # Tidy up the DataFrames of up- and down-regulated genes
        up_clean = up[[gene_col]].copy()  # Use gene_col from config
        up_clean = up_clean.rename(columns={gene_col: "gene"})
        up_clean["source"] = "bie"
        up_clean["thresholds"] = threshold
        up_clean["comparison"] = comparison_str
        up_clean = up_clean[required_columns]

        down_clean = down[[gene_col]].copy()
        down_clean = down_clean.rename(columns={gene_col: "gene"})
        down_clean["source"] = "bie"
        down_clean["thresholds"] = threshold
        down_clean["comparison"] = comparison_str
        down_clean = down_clean[required_columns]  

        all_up.append(up_clean)
        all_down.append(down_clean)

    # Combine up-/down-regulated across all conditions
    df_up = pd.concat(all_up, ignore_index=True).drop_duplicates()
    df_down = pd.concat(all_down, ignore_index=True).drop_duplicates()

    df_up.to_csv(csv_dir / "bie_upregulated_degs.csv", index=False)
    df_down.to_csv(csv_dir / "bie_downregulated_degs.csv", index=False)

    print(f"[Bie et al. Pipeline] Completed. Identified DEGs across {len(config_bie['thresholds']['conditions'])} condition pairs at log2 ≥ {threshold}: {len(df_up)} up, {len(df_down)} down.")  

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_bie_pipeline(full_config)
