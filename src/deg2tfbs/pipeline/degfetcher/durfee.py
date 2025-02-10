"""
--------------------------------------------------------------------------------
<deg2tfbs project>
durfee.py

Module for loading data from Durfee et al., which performed
a transcriptomic study by inducing ppGpp accumulation; curated a list of 
differentially expressed genes within 5 minutes of stringent response onset.

"Transcription profiling of the stringent response in Escherichia coli"
DOI: 10.1128/JB.01092-07

List of up- or down-regulated genes (column 'regulation') is already provided. 
No threshold logic or plot needed.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
from pathlib import Path

import yaml
import numpy as np
import pandas as pd

from deg2tfbs.pipeline.utils import load_dataset


def read_durfee_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Durfee dataset using keys from the YAML config.
    
    E.g. config:
      {
        "dataset_key": "durfee",
        "sheet_name": null,   # or "Sheet1" if needed
        "usecols": [...],
        "header": 0
      }
    The dataset includes columns:
      b_number, gene, gene_product_description,
      Log2_Ratio, Log2_Signal_t0, Log2_Signal_t5,
      regulation ('upregulated' or 'downregulated')
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def run_durfee_pipeline(full_config: dict):
    """
    Reads Durfee data, extracts up/down calls from 'regulation',
    and saves minimal CSVs with columns [gene, source, comparison].

    comparison = "ppGpp_versus_control"
    source = "durfee"

    Output:
      durfee_upregulated_degs.csv
      durfee_downregulated_degs.csv
    """
    # Get config
    config_durfee = full_config['datasets'].get("durfee", None)
    if config_durfee is None:
        print("[Durfee Pipeline] No 'durfee' config found. Skipping.")
        return

    # Load data
    df = read_durfee_data(config_durfee["data"])
    assert "gene" in df.columns, "Expected 'gene' column in Durfee dataset"
    assert "regulation" in df.columns, "Expected 'regulation' column in Durfee dataset"

    # Build output paths
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_durfee["output"]["csv_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)

    # Parse up-/down-regulated
    df_up = df[df["regulation"].str.lower() == "upregulated"].copy()
    df_down = df[df["regulation"].str.lower() == "downregulated"].copy()

    # Tidy final output
    # required_columns = ["gene", "source", "comparison"]
    up_clean = pd.DataFrame({
        "gene": df_up["gene"],
        "source": "durfee",
        "comparison": "ppGpp_versus_control"
    }).drop_duplicates()

    down_clean = pd.DataFrame({
        "gene": df_down["gene"],
        "source": "durfee",
        "comparison": "ppGpp_versus_control"
    }).drop_duplicates()

    # Save
    up_csv = csv_dir / "durfee_upregulated_degs.csv"
    down_csv = csv_dir / "durfee_downregulated_degs.csv"
    up_clean.to_csv(up_csv, index=False)
    down_clean.to_csv(down_csv, index=False)

    print(f"[Durfee et al. Pipeline] Completed. Gathered list of reported DEGs across 1 condition pair: {len(up_clean)} up, {len(down_clean)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_durfee_pipeline(full_config)