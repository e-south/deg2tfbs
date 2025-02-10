#!/usr/bin/env python
"""
--------------------------------------------------------------------------------
<deg2tfbs project>
deg2tfbs/main.py

Coordinates the overall deg2tfbs pipeline.
Supports grouping of DEG CSV files so that you can process:
  - all generated files in a degbatch directory (global),
  - subsets based on a filter (e.g. "up" or "down"),
  - explicit combinations (with customized comparison overrides).

Module Author(s): Eric J. South (refined)
--------------------------------------------------------------------------------
"""

import yaml
import copy
from pathlib import Path

# Import the stage runners
from deg2tfbs.pipeline.degfetcher import run_pipeline_for
from deg2tfbs.pipeline.tffetcher.tffetcher import run_tffetcher_stage
from deg2tfbs.pipeline.tfbsfetcher.tfbsfetcher import run_tfbsfetcher_stage


def load_config(config_file: Path) -> dict:
    with open(config_file, "r") as f:
        return yaml.safe_load(f)


def run_degfetcher(full_config: dict):
    """Run the DEG isolation stage (degfetcher) once for all modules."""
    degfetcher_conf = full_config["pipeline"]["stages"].get("degfetcher", {})
    modules_to_run = degfetcher_conf.get("modules", [])
    for mod_name in modules_to_run:
        # Ensure that mod_name is a non-empty string.
        assert isinstance(mod_name, str) and mod_name.strip(), "Module name must be a non-empty string"
        run_pipeline_for(mod_name, full_config)


def run_downstream(full_config: dict):
    """
    For each group defined in tffetcher.input.deg_csv_groups, update the tffetcher
    and tfbsfetcher configurations and run the downstream stages.
    
    Groups can be defined in one of three ways:
      - A group with an empty definition (or not defined) means use all CSVs in a degbatch subdirectory.
      - A group with a "filter" key to automatically select CSVs by substring.
      - A group with a "files" key listing one or more dictionaries of the form:
           { file: "<csv_filename>", comparison: "<desired_comparison>" }
    """
    base_tffetcher_conf = full_config["pipeline"]["stages"].get("tffetcher", {})
    base_tfbs_conf = full_config["pipeline"]["stages"].get("tfbsfetcher", {})

    input_conf = base_tffetcher_conf.get("input", {})

    # Look for grouping key in config file
    groups = input_conf.get("deg_csv_groups", None)    
    if not isinstance(groups, dict):
        raise ValueError("`deg_csv_groups` must be a mapping (dictionary) of group names to definitions.")

    total = len(groups)
    for idx, (group_name, group_def) in enumerate(groups.items(), start=1):
        # Each group definition can be:
        #   {} or None: process all CSVs.
        #   { filter: "<substring>" }: filter filenames.
        #   { files: [ {file: <name>, comparison: <override>}, ... ] }: explicit list.
        if group_def is None:
            group_def = {}

        # (At this point, group_def is a dict, and its contents will be handled by the stage code.)
        print(f"Processing group {idx} of {total}: '{group_name}'")

        # Deep-copy the base configurations so each group is processed independently.
        group_tffetcher_conf = copy.deepcopy(base_tffetcher_conf)
        group_tfbs_conf = copy.deepcopy(base_tfbs_conf)

        # Inject the group definition into tffetcher configuration.
        # (The tffetcher stage should check for the key 'deg_csv_group' and act accordingly.)
        group_tffetcher_conf.setdefault("input", {})["deg_csv_group"] = group_def

        # Update batch IDs to include the group name so output directories (and files) don’t collide.
        base_tf_batch_id = group_tffetcher_conf.get("batch_id", "tfbatch_default")
        group_tffetcher_conf["batch_id"] = f"{base_tf_batch_id}_{group_name}"
        group_tfbs_conf.setdefault("input", {})["tf_batch_id"] = group_tffetcher_conf["batch_id"]
        base_tfbs_batch_id = group_tfbs_conf.get("batch_id", "tfbsbatch_default")
        group_tfbs_conf["batch_id"] = f"{base_tfbs_batch_id}_{group_name}"

        # Run stage 2: Map DEGs to TFs.
        run_tffetcher_stage(group_tffetcher_conf)
        # Run stage 3: Map TFs to TFBSs.
        run_tfbsfetcher_stage(group_tfbs_conf)
        print(f"Completed group '{group_name}'\n")


def main(config_file: Path):
    full_config = load_config(config_file)

    print("Starting degfetcher stage...")
    run_degfetcher(full_config)
    print("Completed degfetcher stage.\n")

    print("Starting downstream (tffetcher & tfbsfetcher) stages for each group...")
    run_downstream(full_config)
    print("All groups processed.")


if __name__ == "__main__":
    # The config file is assumed to be in deg2tfbs/configs/example.yaml
    config_path = Path(__file__).parent / "configs" / "example.yaml"
    main(config_path)
