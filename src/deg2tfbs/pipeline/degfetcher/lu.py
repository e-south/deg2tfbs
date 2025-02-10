"""
--------------------------------------------------------------------------------
<deg2tfbs project>
lu.py

Module for reading Lu et al., which performed transcriptomic analysis of E. coli 
exposed to 200 mM glyphosate revealed differential expression of 1,040 genes 
(~23% of the genome).

List of up- or down-regulated genes is already provided. No threshold logic
or plot needed.

"Genome-wide transcriptional responses of Escherichia coli to glyphosate, a 
potent inhibitor of the shikimate pathway enzyme 5-enolpyruvylshikimate-3-phosphate
synthase"
DOI: 10.1039/C2MB25374G

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
import yaml
import pandas as pd
from pathlib import Path

from deg2tfbs.pipeline.utils import load_dataset


def read_lu_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Lu dataset using keys from the YAML config.
    
    Source columns: 'Gene name','Gene ID','Functional description','Fold change','Direction'
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def run_lu_pipeline(full_config: dict):
    """
    Reads Lu data, identifies up/down from 'Direction',
    and saves minimal CSV with columns [gene, source, comparison].

    comparison => "glyphosate_shock_versus_control"
    source => "lu"
    """
    config_lu = full_config['datasets'].get("lu", None)
    if config_lu is None:
        print("[Lu Pipeline] No 'lu' config found. Skipping.")
        return

    # Load
    df = read_lu_data(config_lu["data"])

    # Rename 'Gene name' => 'gene'
    if "Gene name" in df.columns:
        df.rename(columns={"Gene name": "gene"}, inplace=True)
    assert "gene" in df.columns, "Expected 'Gene name' in Lu dataset"
    assert "Direction" in df.columns, "Expected 'Direction' in Lu dataset"

    # Build output path
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_lu["output"]["csv_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)

    # Parse up/down
    df_up = df[df["Direction"].str.lower() == "upregulated"].copy()
    df_down = df[df["Direction"].str.lower() == "downregulated"].copy()

    up_clean = pd.DataFrame({
        "gene": df_up["gene"],
        "source": "lu",
        "comparison": "glyphosate_shock_versus_control"
    }).drop_duplicates()

    down_clean = pd.DataFrame({
        "gene": df_down["gene"],
        "source": "lu",
        "comparison": "glyphosate_shock_versus_control"
    }).drop_duplicates()

    # Save
    up_out = csv_dir / "lu_upregulated_degs.csv"
    down_out = csv_dir / "lu_downregulated_degs.csv"
    up_clean.to_csv(up_out, index=False)
    down_clean.to_csv(down_out, index=False)

    print(f"[Lu et al. Pipeline] Completed. Gathered list of reported DEGs across 1 condition pair: {len(up_clean)} up, {len(down_clean)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_lu_pipeline(full_config)