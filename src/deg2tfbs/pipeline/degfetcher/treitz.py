"""
--------------------------------------------------------------------------------
<deg2tfbs project>
treitz.py

Module for processing and analyzing data from Treitz et al, which quantified 
relative protein abundances of MG1655 growing exponentially on minimal medium 
with acetate or glucose as the sole carbon source.

The module loads an Excel file, isolates the necessary columns, plots a volcano plot 
(using log2 fold change versus (-log) p-value), and saves tidy CSV outputs of 
upregulated, downregulated, and fold change data.

"Differential quantitative proteome analysis of Escherichia coli grown on acetate
versus glucose"
DOI: 10.1002/pmic.201600303

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import yaml

from deg2tfbs.pipeline.utils import load_dataset


def read_treitz_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the Treitz dataset using keys from the YAML config.

    Expected config structure:
      {
         "dataset_key": "treitz",
         "sheet_name": "TableSI;Data",
         "usecols": [
            "Gene name",
            "# Proteins",
            "(A/G)-ratio",
            "log2 (A/G)-ratio",
            "Ratio Count",
            "(-log) p-value"
         ],
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
    # Rename columns for clarity and consistency
    df.rename(columns={
        "Gene name": "gene",
        "# Proteins": "proteins",
        "(A/G)-ratio": "ratio",
        "log2 (A/G)-ratio": "log2_fc",
        "Ratio Count": "ratio_count",
        "(-log) p-value": "neg_log_pval"
    }, inplace=True)
    return df


def run_treitz_pipeline(full_config: dict):
    """
    Loads the Treitz dataset, applies thresholds, plots a volcano plot, and saves results as CSVs.

    Expected YAML config snippet for treitz:
      treitz:
        data:
          dataset_key: "treitz"
          sheet_name: "TableSI;Data"
          usecols:
            - "Gene name"
            - "# Proteins"
            - "(A/G)-ratio"
            - "log2 (A/G)-ratio"
            - "Ratio Count"
            - "(-log) p-value"
          header: 1
        thresholds:
          log2_fc_threshold: 1.5
        output:
          csv_subdir: "csvs"
          plot_subdir: "plots"
    """
    config_treitz = full_config['datasets']["treitz"]
    df = read_treitz_data(config_treitz["data"])

    # Build output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_treitz["output"]["csv_subdir"]
    plot_dir = batch_dir / config_treitz["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Retrieve threshold for log2 fold change (asserting numeric type for robustness)
    threshold = config_treitz["thresholds"]["log2_fc_threshold"]
    assert isinstance(threshold, (int, float)), "Threshold must be numeric."

    # Assign colors for the volcano plot: red for up, green for down, lightgray otherwise
    df["color"] = np.where(
        df["log2_fc"] >= threshold, "red",
        np.where(df["log2_fc"] <= -threshold, "green", "lightgray")
    )

    # Set plot font to Arial
    plt.rcParams['font.family'] = 'Arial'

    # Create a volcano plot: x-axis is log2 fold change and y-axis is (-log) p-value
    plt.figure(figsize=(6, 5))
    plt.scatter(
        df["log2_fc"], df["neg_log_pval"],
        c=df["color"], alpha=0.25, edgecolors="none"
    )
    plt.axvline(x=threshold, color='gray', linestyle='--')
    plt.axvline(x=-threshold, color='gray', linestyle='--')
    plt.title("Protein abundances (M9-acetate versus M9-glucose)")
    plt.xlabel("log2(Acetate/Glucose)")
    plt.ylabel("-log p-value")
    sns.despine()

    plot_path = plot_dir / "treitz_volcano.png"
    plt.savefig(plot_path, dpi=300)
    plt.close()

    # Filter up- and down-regulated genes based on the log2_fc threshold
    upregulated = df[df["log2_fc"] >= threshold].copy()
    downregulated = df[df["log2_fc"] <= -threshold].copy()

    # Clean up the gene column: remove missing or ambiguous entries
    upregulated.dropna(subset=["gene"], inplace=True)
    upregulated = upregulated[upregulated["gene"].astype(str).str.strip().ne('?')]
    downregulated.dropna(subset=["gene"], inplace=True)
    downregulated = downregulated[downregulated["gene"].astype(str).str.strip().ne('?')]

    # Create tidy DataFrames for up- and down-regulated genes with metadata
    required_columns = ["gene", "source", "thresholds", "comparison"]
    comparison_str = "Acetate_over_Glucose"

    up_clean = upregulated.copy()
    up_clean["source"] = "treitz"
    up_clean["thresholds"] = threshold
    up_clean["comparison"] = comparison_str
    up_clean = up_clean[required_columns]

    down_clean = downregulated.copy()
    down_clean["source"] = "treitz"
    down_clean["thresholds"] = threshold
    down_clean["comparison"] = comparison_str
    down_clean = down_clean[required_columns]

    # Save tidy lists of up- and down-regulated genes as CSV files
    up_csv_path = csv_dir / "treitz_upregulated_degs.csv"
    down_csv_path = csv_dir / "treitz_downregulated_degs.csv"
    up_clean.to_csv(up_csv_path, index=False)
    down_clean.to_csv(down_csv_path, index=False)

    # Save a CSV with fold change and log fold change values for all genes
    fold_change_data = df[["gene", "ratio", "log2_fc"]].dropna(subset=["gene"]).copy()
    fold_change_csv_path = csv_dir / "treitz_fold_change_data.csv"
    fold_change_data.to_csv(fold_change_csv_path, index=False)

    print(f"[Treitz Pipeline] Completed. Processed {len(df)} entries with threshold log2_fc ≥ {threshold}.")


if __name__ == "__main__":
    # Load configuration from the example YAML file
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_treitz_pipeline(full_config)
