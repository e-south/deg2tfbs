"""
--------------------------------------------------------------------------------
<deg2tfbs project>
bie.py

Module for reading and analyzing data described in Bie et al.,
an RNA‐seq survey of how E. coli K‐12 MG1655 responds to multiple antibiotics.

We compute average FPKM across replicates, then compute log2 fold-change 
for antibiotic vs. water, returning up- and down-regulated genes based 
on a user-defined threshold (and optionally significance columns).
--------------------------------------------------------------------------------
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from deg2tfbs.pipeline.dataloader.utils import load_dataset


def read_bie_data(config_data: dict) -> pd.DataFrame:
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def compute_average_fpkm(df, prefix, replicate_suffixes):
    """
    Compute the average FPKM across replicates for a given prefix.
    e.g., prefix="KAN", replicate_suffixes=["_1_fpkm","_2_fpkm","_3_fpkm"]
    """
    cols = [f"{prefix}{suffix}" for suffix in replicate_suffixes]
    valid_cols = [c for c in cols if c in df.columns]
    if not valid_cols:
        raise KeyError(f"No valid columns found for prefix '{prefix}' with replicates {replicate_suffixes}")
    return df[valid_cols].mean(axis=1)


def add_average_fpkm_columns(df, config_bie):
    """
    For each antibiotic + the water prefix, compute average FPKM 
    using the replicate_dict from config. 
    E.g. replicate_dict might be:
        { "KAN": ["_1_fpkm","_2_fpkm","_3_fpkm"],
          "CIP": ["_2_fpkm","_3_fpkm","_4_fpkm"],
          "H2O": ["_2_fpkm","_3_fpkm","_4_fpkm"] }
    """
    replicate_dict = config_bie["thresholds"].get("replicate_dict", {})
    water_prefix = config_bie["thresholds"].get("water_prefix", "H2O")
    antibiotics = config_bie["thresholds"].get("antibiotics", [])

    # Gather all relevant prefixes (antibiotics plus water)
    all_prefixes = antibiotics + [water_prefix]

    for abx in all_prefixes:
        # If there's a custom replicate list for this prefix, use it.
        # Otherwise, fallback to something like ["_2_fpkm","_3_fpkm","_4_fpkm"].
        replicate_suffixes = replicate_dict.get(abx, ["_2_fpkm","_3_fpkm","_4_fpkm"])
        
        avg_col = f"{abx}_avg_fpkm"
        df[avg_col] = compute_average_fpkm(df, abx, replicate_suffixes)

    return df


def filter_and_rank_genes(
    df,
    abx_prefix,
    water_prefix,
    log2_fc_col,
    sig_col,
    threshold=2.0
):
    """
    Filter and rank up/down regulated genes for antibiotic vs. water 
    based on log2 fold change and significance columns.

    This function:
    1) Subsets rows that meet a condition (log2_fc >= threshold or <= -threshold, if you want).
    2) (Optional) Uses 'sig_col' to confirm significance if that is part of your logic.
    3) Calculates an effect size or rank based on control FPKM + log2 FC.

    Returns: A tuple (up_df, down_df) with relevant columns.
    """
    # e.g. "KAN_avg_fpkm"
    abx_col = f"{abx_prefix}_avg_fpkm"
    h2o_col = f"{water_prefix}_avg_fpkm"

    # Optionally, we might require sig_col to say 'UP' or 'DOWN' 
    # or use numerical p-value columns. 
    # For now, let's assume "UP" means up, "DOWN" means down 
    # or skip if you only want fold-change.
    sub = df.copy()

    # Up = log2_fc >= threshold
    up = sub[
        (sub[log2_fc_col] >= threshold) &
        (sub[sig_col] == "UP")     # <--- or (True) if you don't need significance
    ].copy()

    # Down = log2_fc <= -threshold
    down = sub[
        (sub[log2_fc_col] <= -threshold) &
        (sub[sig_col] == "DOWN")   # <--- or (True) if you don't need significance
    ].copy()

    # Optionally rank them, e.g. "big movers" with high water expression + large fold change
    # We can add some "Effect_Size" metric
    if abx_col in up.columns and h2o_col in up.columns:
        up["FPKM_Difference"] = up[abx_col] - up[h2o_col]
        up["Rank_FPKM"] = up[h2o_col].rank(ascending=False)
        up["Rank_Log2FC"] = up[log2_fc_col].rank(ascending=False)  # up means bigger is better
        up["Effect_Size"] = up["Rank_FPKM"] + up["Rank_Log2FC"]

    if abx_col in down.columns and h2o_col in down.columns:
        down["FPKM_Difference"] = down[abx_col] - down[h2o_col]
        down["Rank_FPKM"] = down[h2o_col].rank(ascending=False)
        down["Rank_Log2FC"] = down[log2_fc_col].rank(ascending=True)  # down means smaller is better
        down["Effect_Size"] = down["Rank_FPKM"] + down["Rank_Log2FC"]

    # Sort by effect size
    if "Effect_Size" in up.columns:
        up.sort_values("Effect_Size", inplace=True, ascending=True)
    if "Effect_Size" in down.columns:
        down.sort_values("Effect_Size", inplace=True, ascending=True)

    return up, down


def bie_ma_plot(
    df: pd.DataFrame,
    condition_col: str,
    reference_col: str,
    threshold: float,
    output_path: Path = None
):
    """
    Creates and saves an MA-like plot comparing two columns: antibiotic vs. water,
    coloring points that exceed log2_fc >= threshold (red) or <= -threshold (green).
    Assumes df has 'log2_fc' and 'avg_expr' columns from antibiotic vs. water.

    Args:
        df (pd.DataFrame): DataFrame with columns [condition_col, reference_col, log2_fc, avg_expr].
        condition_col (str): e.g. "KAN_avg_fpkm"
        reference_col (str): e.g. "H2O_avg_fpkm"
        threshold (float): log2 fold-change cutoff.
        output_path (Path, optional): Save the figure if provided.
    """
    sns.set_style("ticks")

    colors = np.where(
        df["log2_fc"] >= threshold, "red",
        np.where(df["log2_fc"] <= -threshold, "green", "gray")
    )

    plt.figure(figsize=(6,5))
    plt.scatter(df["avg_expr"], df["log2_fc"], c=colors, alpha=0.25, edgecolors="none")
    plt.xscale("log")
    plt.axhline(threshold, color="gray", linestyle="--")
    plt.axhline(-threshold, color="gray", linestyle="--")

    abx_name = condition_col.replace("_avg_fpkm","")
    water_name = reference_col.replace("_avg_fpkm","")

    plt.title(f"BIE et al.\n{abx_name} vs. {water_name}")
    plt.xlabel("Average Expression")
    plt.ylabel("Log2 Fold Change")
    sns.despine()

    if output_path is not None:
        plt.savefig(output_path, dpi=150)
    plt.close()


def run_bie_pipeline(full_config: dict):
    """
    Orchestrates:

      1) Reading the Bie dataset from config
      2) Computing average FPKM for antibiotic vs. water
      3) For each antibiotic in config, compute log2(FC) vs. water
      4) Identify up/down regulated genes
      5) Generate MA plots
      6) Save results to CSV
    """
    config_bie = full_config.get("bie", None)
    if config_bie is None:
        print("[BIE Pipeline] No 'bie' config found. Skipping.")
        return

    # 1) Load data
    df = read_bie_data(config_bie["data"])

    # 2) Build output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_bie["output"]["csv_subdir"]
    plot_dir = batch_dir / config_bie["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # 3) Setup thresholds
    threshold = config_bie["thresholds"].get("log2_fc_threshold", 2.0)
    antibiotic_list = config_bie["thresholds"].get("antibiotics", ["KAN","CIP"])
    water_prefix = config_bie["thresholds"].get("water_prefix","H2O")

    # 4) Compute average FPKM columns
    df = add_average_fpkm_columns(df, config_bie)

    all_up = []
    all_down = []

    # 5) For each antibiotic, do the comparison vs. water
    sig_cols = config_bie["thresholds"].get("significance_cols", {})
    fc_cols = config_bie["thresholds"].get("log2_fc_cols", {})

    for abx in antibiotic_list:
        log2_fc_col = fc_cols.get(abx, None)
        sig_col = sig_cols.get(abx, None)

        # If we have precomputed log2 fold-change + significance:
        if log2_fc_col and log2_fc_col in df.columns and sig_col in df.columns:
            up, down = filter_and_rank_genes(
                df=df,
                abx_prefix=abx,
                water_prefix=water_prefix,
                log2_fc_col=log2_fc_col,
                sig_col=sig_col,
                threshold=threshold
            )
            # For MA plot, let's rename columns so it matches [log2_fc, avg_expr].
            plot_df = df[[log2_fc_col]].copy()
            plot_df["avg_expr"] = (
                df[f"{abx}_avg_fpkm"] + df[f"{water_prefix}_avg_fpkm"]
            ) / 2.0
            plot_df.rename(columns={log2_fc_col: "log2_fc"}, inplace=True)

        else:
            # Otherwise compute on the fly
            df_comp = df[[f"{abx}_avg_fpkm", f"{water_prefix}_avg_fpkm"]].dropna().copy()
            df_comp["log2_fc"] = (
                np.log2(df_comp[f"{abx}_avg_fpkm"]) - 
                np.log2(df_comp[f"{water_prefix}_avg_fpkm"])
            )
            df_comp["avg_expr"] = (
                df_comp[f"{abx}_avg_fpkm"] + df_comp[f"{water_prefix}_avg_fpkm"]
            ) / 2.0

            up = df_comp[df_comp["log2_fc"] >= threshold].copy()
            down = df_comp[df_comp["log2_fc"] <= -threshold].copy()
            plot_df = df_comp

        # Save partial CSVs
        up.to_csv(csv_dir / f"bie_up_{abx}_vs_{water_prefix}.csv", index=False)
        down.to_csv(csv_dir / f"bie_down_{abx}_vs_{water_prefix}.csv", index=False)
        all_up.append(up)
        all_down.append(down)

        # 6) Create MA plot
        plot_path = plot_dir / f"bie_ma_{abx}_vs_{water_prefix}.png"
        bie_ma_plot(
            df=plot_df,
            condition_col=f"{abx}_avg_fpkm",
            reference_col=f"{water_prefix}_avg_fpkm",
            threshold=threshold,
            output_path=plot_path
        )

    # 7) Combine all up/down across conditions
    df_up = pd.concat(all_up, ignore_index=True).drop_duplicates()
    df_down = pd.concat(all_down, ignore_index=True).drop_duplicates()

    df_up.to_csv(csv_dir / "bie_DEGs_upregulated_all.csv", index=False)
    df_down.to_csv(csv_dir / "bie_DEGs_downregulated_all.csv", index=False)

    print(f"[Bie Pipeline] Completed. Found {len(df_up)} total up and {len(df_down)} total down.")