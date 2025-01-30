"""
--------------------------------------------------------------------------------
<deg2tfbs project>
zhang.py

Module for reading Zhang et al. dataset, which performed RNA-seq on E. coli expressing
σ32-I54N, a mutant form of the heat-shock sigma factor σ32, and WT σ32 at 42°C.


List of up-regulated genes (no downregulated) is already provided. No threshold logic
or plot needed.

"Heat-Shock Response Transcriptional Program Enables High-Yield and 
High-Quality Recombinant Protein Production in Escherichia coli"
DOI: 10.1021/cb5004477

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
import yaml
import pandas as pd
from pathlib import Path

from deg2tfbs.pipeline.degfetcher.utils import load_dataset

def read_zhang_data(config_data: dict) -> pd.DataFrame:
    """
    Reads from 'zhang' key in the YAML, 
    referencing TableS1.xlsx, sheet 'Sheet1', header=0
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df

def run_zhang_pipeline(full_config: dict):
    """
    Reads Zhang data, 
    which are all reported to be upregulated (no down).
    """
    config_zhang = full_config.get("zhang", None)
    if config_zhang is None:
        print("[Zhang Pipeline] No 'zhang' config found. Skipping.")
        return

    # Load
    df = read_zhang_data(config_zhang["data"])
    if "Gene" in df.columns:
        df.rename(columns={"Gene": "gene"}, inplace=True)
    assert "gene" in df.columns, "Expect 'Gene' or 'gene' in Zhang data"

    # Build output path
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_zhang["output"]["csv_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)

    # All are upregulated
    up_clean = pd.DataFrame({
        "gene": df["gene"],
        "source": "zhang",
        "comparison": "sigma32-I54N_expression_versus_control"
    }).drop_duplicates()

    # Save only up. No down CSV needed
    up_csv = csv_dir / "zhang_upregulated_degs.csv"
    up_clean.to_csv(up_csv, index=False)

    print(f"[Zhang et al. Pipeline] Completed. Gathered list of reported DEGs across 1 condition pair: {len(up_clean)} up, 0 down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_zhang_pipeline(full_config)