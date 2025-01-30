"""
--------------------------------------------------------------------------------
<deg2tfbs project>
emani.py

Module for loading and analyzing data described in Emani et al., which investigated
the relationship between growth and protein secretion in E. coli Nissle 1917 by
secreting super-folder green fluorescent protein (sfGFP) through the native curli 
secretion machinery. RNA sequencing data from batch culture experiments comparing
the all_tags and non-tagged GFP production strains provide omic snapshots of 
periplasmic stress.

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and plots a volcano figure (MinusLog10PAd vs. log2 Fold Change).

"Periplasmic stress contributes to a trade-off between protein secretion and cell 
growth in Escherichia coli Nissle 1917"
DOI: 10.1093/synbio/ysad013

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import yaml
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from deg2tfbs.pipeline.degfetcher.utils import load_dataset

def read_emani_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Emani dataset using keys from the YAML config.
    
    E.g. config:
      {
        "dataset_key": "emani_et_al",
        "sheet_name": "Figure 5_alltags v control",
        "usecols": [
            "Gene",
            "log2FoldChange (Green = Upregulated in 10uM IPTG All Tags)",
            "pvalue",
            "padj"
        ]
      }
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    # Rename columns for clarity
    df.rename(columns={
        "log2FoldChange (Green = Upregulated in 10uM IPTG All Tags)": "Log2FoldChange",
        "pvalue": "PValue",
        "padj": "PAdj"
    }, inplace=True)

    # Convert PAdj to -log10 scale for plotting
    df["MinusLog10PAdj"] = -np.log10(df["PAdj"])
    return df

def run_emani_pipeline(full_config: dict):
    """
    Loads Emani dataset, applies thresholds, and stores results.

    Expected structure in the YAML:
      emani:
        data:
          dataset_key: "emani"
          sheet_name: "Figure 5_alltags v control"
          usecols:
            - "Gene"
            - "log2FoldChange (Green = Upregulated in 10uM IPTG All Tags)"
            - "pvalue"
            - "padj"
        thresholds:
          log2_fc_threshold: 2.0
        output:
          csv_subdir: "csv"
          plot_subdir: "plots"
    """
    config_emani = full_config["emani"]
    df = read_emani_data(config_emani["data"])

    # Build output paths for CSV and plots
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_emani["output"]["csv_subdir"]
    plot_dir = batch_dir / config_emani["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Threshold for up/down regulation
    threshold = config_emani["thresholds"]["log2_fc_threshold"]
    
    # Up- and down-regulated genes are colored red and green, respectively
    df["color"] = np.where(
        df["Log2FoldChange"] >= threshold, "red",
        np.where(df["Log2FoldChange"] <= -threshold, "green", "lightgray")
    )

    # Plot a "volcano" style figure with -log10(PAdj) versus Log2FoldChange
    plt.figure(figsize=(8, 6))
    plt.scatter(
        df["Log2FoldChange"], df["MinusLog10PAdj"],
        c=df["color"], alpha=0.25, edgecolors="none"
    )
    plt.axvline(x=threshold, color='gray', linestyle='--')
    plt.axvline(x=-threshold, color='gray', linestyle='--')
    plt.title("Emani et al.\n -log10(PAdj) versus (N22_GFP / WT)")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-log10(PAdj)")
    sns.despine()

    plot_path = plot_dir / "emani_N22_GFP_versus_WT.png"
    plt.savefig(plot_path, dpi=150)
    plt.close()

    # Filter up/down
    upregulated = df[df["Log2FoldChange"] >= threshold].copy()
    downregulated = df[df["Log2FoldChange"] <= -threshold].copy()

    # Clean up gene column
    upregulated.dropna(subset=["Gene"], inplace=True)
    upregulated = upregulated[upregulated["Gene"].astype(str).str.strip().ne('?')]
    downregulated.dropna(subset=["Gene"], inplace=True)
    downregulated = downregulated[downregulated["Gene"].astype(str).str.strip().ne('?')]

    required_columns = ["gene", "source", "thresholds", "comparison"]
    
    # Define comparison
    target_condition = "all_tags"
    reference_condition = "control"
    comparison_str = f"{target_condition}_vs_{reference_condition}"
    
    # Tidy up the DataFrames of up- and down-regulated genes
    up_clean = upregulated.copy()
    up_clean["gene"] = up_clean["Gene"]
    up_clean["source"] = "emani"
    up_clean["thresholds"] = threshold
    up_clean["comparison"] = comparison_str
    up_clean = up_clean[required_columns]

    down_clean = downregulated.copy()
    down_clean["gene"] = down_clean["Gene"] 
    down_clean["source"] = "emani"
    down_clean["thresholds"] = threshold
    down_clean["comparison"] = comparison_str
    down_clean = down_clean[required_columns]

    # Save filtered data
    up_clean.to_csv(csv_dir / "emani_upregulated_degs.csv", index=False)
    down_clean.to_csv(csv_dir / "emani_downregulated_degs.csv", index=False)

    print(f"[Emani et al. Pipeline] Completed. Identified DEGs across 1 condition pair at log2 ≥ {threshold}: {len(upregulated)} up, {len(downregulated)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_emani_pipeline(full_config)