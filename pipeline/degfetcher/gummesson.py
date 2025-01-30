"""
--------------------------------------------------------------------------------
<deg2tfbs project>
gummesson.py

Module for reading Gummesson et al. data, which performed RNA-sequencing following 
isoleucine starvation in E. coli. The dataset includes up- and down-regulated genes
after 80 minutes.

List of up- or down-regulated genes is already provided. No threshold logic
or plot needed.

"Valine‐Induced Isoleucine Starvation in E. coli Studied by Spike‐In 
Normalized RNA Sequencing"
DOI: 10.3389/fgene.2020.00144

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
import yaml
import pandas as pd
from pathlib import Path

from deg2tfbs.pipeline.degfetcher.utils import load_dataset

def read_gummesson_data(config_data: dict, sheet_name: str) -> pd.DataFrame:
    """
    Reads the sourced Gummesson dataset using keys from the YAML config.
    
    E.g. config:
      {
        "dataset_key": "gummesson",
        "header": 3,  # columns start at row 4
        ...
      }
    We rename 'Name' => 'gene'.
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=sheet_name,
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    # rename 'Name' => 'gene' 
    if "Name" in df.columns:
        df.rename(columns={"Name": "gene"}, inplace=True)
    else:
        assert "gene" in df.columns, "Expect 'Name' or 'gene' in Gummesson data"

    return df


def run_gummesson_pipeline(full_config: dict):
    """
    - Reads from 'TableS8' => upregulated
    - Reads from 'TableS9' => downregulated
    - Keeps only [gene, source, comparison]
    - Writes gummesson_upregulated_degs.csv, gummesson_downregulated_degs.csv

    comparison => "isoleucine_starvation_versus_control"
    source => "gummesson"
    """
    config_gum = full_config.get("gummesson", None)
    if config_gum is None:
        print("[Gummesson Pipeline] No 'gummesson' config found. Skipping.")
        return

    # Build output path
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_gum["output"]["csv_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)

    # Read TableS8 => up
    df_up = read_gummesson_data(config_gum["data"], sheet_name="Table_S8")
    up_clean = pd.DataFrame({
        "gene": df_up["gene"],
        "source": "gummesson",
        "comparison": "isoleucine_starvation_versus_control"
    }).drop_duplicates()

    # Read TableS9 => down
    df_down = read_gummesson_data(config_gum["data"], sheet_name="Table_S9")
    down_clean = pd.DataFrame({
        "gene": df_down["gene"],
        "source": "gummesson",
        "comparison": "isoleucine_starvation_versus_control"
    }).drop_duplicates()

    # Save
    up_csv = csv_dir / "gummesson_upregulated_degs.csv"
    down_csv = csv_dir / "gummesson_downregulated_degs.csv"
    up_clean.to_csv(up_csv, index=False)
    down_clean.to_csv(down_csv, index=False)

    print(f"[Gummesson et al. Pipeline] Completed. Gathered list of reported DEGs across 1 condition pair: {len(up_clean)} up, {len(down_clean)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_gummesson_pipeline(full_config)