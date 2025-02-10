"""
--------------------------------------------------------------------------------
<deg2tfbs project>
wu.py

Module for loading and analyzing data from Wu et al., which investigated the 
allocation of the proteome of E. coli grown in rich or minimal media during steady-state
growth and transitions between them. Proteomic fractions for many genes are available in
Supplementary Table 9, specifically samples M1 (Carbon limited), O8 (LB), and N5
(nutrient-rich WT). See Supplemental Table 2 for details about environmental conditions
for each proteomic sample.

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and saves an MA plot (average expression vs. log2 Fold Change).

"Enzyme expression kinetics by Escherichia coli during transition from rich to 
minimal media depends on proteome reserves"
DOI: 10.1038/s41564-022-01310-w

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""
from pathlib import Path

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from deg2tfbs.pipeline.utils import load_dataset


def read_wu_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Wu dataset using keys from the YAML config.

    E.g. config:
      {
        "dataset_key": "wu",
        "sheet_name": "Group #1 protein mass fractions"
      }
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    # Drop rows with any missing values
    df = df.dropna()
    return df

def wu_pairwise_comparison(df, comparisons, threshold=2.0, plot_dir=None):
    """
    comparisons is a dict like:
      {
        "Clim_MOPS_Over_Rich_MOPS": ("M1","N5"),
        "Clim_MOPS_Over_WT_Rich_LB": ("M1","O8")
      }

    For each comparison key, we do log2( numerator / denominator ).
    Return aggregated up and down sets.
    """
    all_up = []
    all_down = []

    for label, (num, den) in comparisons.items():
        sub = df[[num, den]].copy()
        # Handle zeros before calculations
        sub = df[["Gene name", num, den]].copy()
        sub = sub.replace(0, np.nan).dropna()
        # Add pseudocount to avoid log(0)
        epsilon = 1e-5
        sub["log2_fc"] = np.log2((sub[num]+epsilon)/(sub[den]+epsilon))
        sub["avg_expr"] = (sub[num] + sub[den])/2
        sub["comparison"] = label
        
        # Up- and down-regulated genes are colored red and green, respectively
        sub["color"] = np.where(sub["log2_fc"] >= threshold, "red",
                      np.where(sub["log2_fc"] <= -threshold, "green", "gray"))

        up = sub[sub["log2_fc"] >= threshold]
        down = sub[sub["log2_fc"] <= -threshold]
        all_up.append(up)
        all_down.append(down)
         
        if plot_dir:
            sns.set_style("ticks")
            plt.figure(figsize=(6,5))
            plt.scatter(sub["avg_expr"], sub["log2_fc"],
                        c=sub["color"], alpha=0.25, edgecolors='none')
            plt.xscale('log')
            plt.axhline(threshold, color='gray', linestyle='--')
            plt.axhline(-threshold, color='gray', linestyle='--')
            plt.title(f"Wu et al.\n log2 ({num} / {den})\n{label}")
            plt.xlabel("Average Protein Mass Fraction")
            plt.ylabel("Log2 Fold Change")
            sns.despine(top=True, right=True)
            plt.savefig(plot_dir / f"wu_{label}.png", dpi=300)
            plt.close()

    df_up = pd.concat(all_up, ignore_index=True).drop_duplicates()
    df_down = pd.concat(all_down, ignore_index=True).drop_duplicates()
    return df_up, df_down

def run_wu_pipeline(full_config: dict):
    """
    Example config in YAML:
      wu:
        data:
          dataset_key: "wu"
          sheet_name: "Group #1 protein mass fractions"
        thresholds:
          log2_fc_threshold: 2.0
          comparisons:
            Clim_MOPS_Over_Rich_MOPS: ["M1", "N5"]
            Clim_MOPS_Over_WT_Rich_LB: ["M1", "O8"]
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """
    config_wu = full_config['datasets']["wu"]
    df = read_wu_data(config_wu["data"])

    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_wu["output"]["csv_subdir"]
    plot_dir = batch_dir / config_wu["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    threshold = config_wu["thresholds"]["log2_fc_threshold"]
    comparisons = config_wu["thresholds"]["comparisons"]

    up_all, down_all = wu_pairwise_comparison(df, comparisons, threshold=threshold, plot_dir=plot_dir)
    
    required_columns = ["gene", "source", "thresholds", "comparison"]
    
    # Tidy up the DataFrames of up- and down-regulated genes
    up_clean = up_all[["Gene name", "comparison"]].copy()
    up_clean.rename(columns={"Gene name": "gene"}, inplace=True)
    up_clean["source"] = "wu"
    up_clean["thresholds"] = threshold
    up_clean = up_clean[required_columns]

    down_clean = down_all[["Gene name", "comparison"]].copy()
    down_clean.rename(columns={"Gene name": "gene"}, inplace=True)
    down_clean["source"] = "wu"
    down_clean["thresholds"] = threshold
    down_clean = down_clean[required_columns]

    up_clean.to_csv(csv_dir / "wu_upregulated_degs.csv", index=False)
    down_clean.to_csv(csv_dir / "wu_downregulated_degs.csv", index=False)

    print(f"[Wu et al. Pipeline] Completed. Identified DEGs across {len(comparisons)} condition pairs at log2 ≥ {threshold}: {len(up_all)} up, {len(down_all)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_wu_pipeline(full_config)
