"""
--------------------------------------------------------------------------------
<deg2tfbs project>
emani.py

Module for reading and analyzing data described in Emani et al., which investigated
the relationship between growth and protein secretion in E. coli Nissle 1917 by
secreting super-folder green fluorescent protein (sfGFP) through the native curli 
secretion machinery. RNA sequencing data from batch culture experiments comparing
the all_tags and N22_GFP variants provide further evidence of periplasmic stress as
a mechanism for growth rate depression. Their results from RNA sequencing point to
specific pathways that maybe useful targets for mitigating maladaptive stress response 
for secretion systems.

"Periplasmic stress contributes to a trade-off between protein secretion and cell 
growth in Escherichia coli Nissle 1917"
DOI: 10.1093/synbio/ysad013

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_emani_data(config_data: dict) -> pd.DataFrame:
    """
    Loads Emani data from the YAML config, e.g.:

        "dataset_key": "emani_et_al",
        "sheet_name": "Figure 5_alltags v control",
        "usecols": [
            "Gene",
            "log2FoldChange (Green = Upregulated in 10uM IPTG All Tags)",
            "pvalue",
            "padj"
        ]

    Returns a DataFrame with columns [Gene, Log2FoldChange, PValue, PAdj, MinusLog10PAdj].
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
    Orchestrates loading Emani dataset, applying threshold, and storing results.

    Expected structure in the YAML:
      emani:
        data:
          dataset_key: "emani_et_al"
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

    # Build output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_emani["output"]["csv_subdir"]
    plot_dir = batch_dir / config_emani["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Threshold
    threshold = config_emani["thresholds"]["log2_fc_threshold"]
    df["color"] = np.where(
        df["Log2FoldChange"] >= threshold, "red",
        np.where(df["Log2FoldChange"] <= -threshold, "green", "lightgray")
    )

    # Plot a "volcano" style figure
    plt.figure(figsize=(8, 6))
    plt.scatter(
        df["Log2FoldChange"], df["MinusLog10PAdj"],
        c=df["color"], alpha=0.25, edgecolors="none"
    )
    plt.axvline(x=threshold, color='gray', linestyle='--')
    plt.axvline(x=-threshold, color='gray', linestyle='--')
    plt.title("Periplasmic Stress Genes: log2FC vs. -log10(PAdj)")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-log10(PAdj)")
    sns.despine()

    plot_path = plot_dir / "emani_DEGs_et_al.png"
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

    # Save
    upregulated.to_csv(csv_dir / "DEGs_upregulated_Emani_et_al.csv", index=False)
    downregulated.to_csv(csv_dir / "DEGs_downregulated_Emani_et_al.csv", index=False)

    print(f"[Emani Pipeline] Completed. Found: {len(upregulated)} up, {len(downregulated)} down total.")
