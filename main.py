"""
--------------------------------------------------------------------------------
<deg2tfbs project>
deg2tfbs/main.py

Module Author(s): Eric J. South
Dunlop Lab

Coordinates the overall pipelines.
--------------------------------------------------------------------------------
"""

import yaml
from pathlib import Path

# Import the dynamic runner from degfetcher
from deg2tfbs.pipeline.degfetcher import run_pipeline_for
from deg2tfbs.pipeline.tffetcher.tffetcher import run_tffetcher_stage
from deg2tfbs.pipeline.tfbsfetcher.tfbsfetcher import run_tfbsfetcher_stage

def main(config_file: Path):
    with open(config_file, "r") as f:
        full_config = yaml.safe_load(f)

    # Step 1. Isolate DEGs
    degfetcher_conf = full_config["pipeline"]["stages"].get("degfetcher", {})
    modules_to_run = degfetcher_conf.get("modules", [])

    for mod_name in modules_to_run:
        run_pipeline_for(mod_name, full_config)

    # Step 2. Map DEGs to TFs
    tffetcher_conf = full_config["pipeline"]["stages"].get("tffetcher", {})
    if tffetcher_conf:
        run_tffetcher_stage(tffetcher_conf)

    # Step 2. Map TFs to TFBSs
    tfbs_conf = full_config["pipeline"]["stages"].get("tfbsfetcher", {})
    if tfbs_conf:
        run_tfbsfetcher_stage(tfbs_conf)
 
 
if __name__ == "__main__":
    config_path = Path(__file__).parent / "configs" / "example.yaml"
    main(config_path)
