"""
--------------------------------------------------------------------------------
<deg2tfbs project>
houser.py

Module for reading Houser et al., which performed a time‐course characterization 
(two weeks) of E. coli growth and starvation.

List of up- or down-regulated genes is already provided. No threshold logic
or plot needed.

"Controlled Measurement and Comparative Analysis of Cellular Components in E. coli"
DOI: 10.1371/journal.pcbi.1004400

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
import yaml
import pandas as pd
from pathlib import Path

from deg2tfbs.pipeline.degfetcher.utils import load_dataset

def read_houser_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Houser dataset using keys from the YAML config.
    
    Source columns: 'Common Name','Entrez ID','Gene Name','Prediction','Sorted Category'
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def run_houser_pipeline(full_config: dict):
    """
    Reads Houser data, identifies up/down from 'Sorted Category',
    and saves minimal CSV with columns [gene, source, comparison].

    comparison => "two_week_starvation_versus_control"
    source => "houser"
    """
    config_houser = full_config.get("houser", None)
    if config_houser is None:
        print("[Houser Pipeline] No 'houser' config found. Skipping.")
        return

    # Load
    df = read_houser_data(config_houser["data"])

    # Rename 'Common Name' => 'gene'
    if "Common Name" in df.columns:
        df.rename(columns={"Common Name": "gene"}, inplace=True)
    assert "gene" in df.columns, "Expect 'gene' or 'Common Name' in Houser data"
    assert "Sorted Category" in df.columns, "Expect 'Sorted Category' in Houser data"

    # Build output path
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_houser["output"]["csv_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)

    # Parse up/down
    df_up = df[df["Sorted Category"] == "Up Regulated"].copy()
    df_down = df[df["Sorted Category"] == "Down Regulated"].copy()

    up_clean = pd.DataFrame({
        "gene": df_up["gene"],
        "source": "houser",
        "comparison": "two_week_starvation_versus_control"
    }).drop_duplicates()

    down_clean = pd.DataFrame({
        "gene": df_down["gene"],
        "source": "houser",
        "comparison": "two_week_starvation_versus_control"
    }).drop_duplicates()

    # Save
    up_out = csv_dir / "houser_upregulated_degs.csv"
    down_out = csv_dir / "houser_downregulated_degs.csv"
    up_clean.to_csv(up_out, index=False)
    down_clean.to_csv(down_out, index=False)

    print(f"[Houser et al. Pipeline] Completed. Gathered list of reported DEGs across 1 condition pair: {len(up_clean)} up, {len(down_clean)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_houser_pipeline(full_config)