"""
--------------------------------------------------------------------------------
<deg2tfbs project>
deg2tfbs/main.py

Coordinates the overall deg2tfbs pipeline.
Supports grouping of DEG CSV files so that you can process:
  - all generated files in a degbatch directory (global),
  - subsets based on a filter (e.g. "up" or "down"),
  - explicit combinations (with customized comparison overrides).

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import yaml
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
    """
    base_tffetcher_conf = full_config["pipeline"]["stages"].get("tffetcher", {})
    base_tfbs_conf = full_config["pipeline"]["stages"].get("tfbsfetcher", {})

    print("Processing consolidated deg_csv_groups:")
    run_tffetcher_stage(base_tffetcher_conf)
    run_tfbsfetcher_stage(base_tfbs_conf)
    print("Completed downstream processing.\n")


def main(config_file: Path):
    full_config = load_config(config_file)

    print("Starting degfetcher stage...")
    run_degfetcher(full_config)
    print("Completed degfetcher stage.\n")

    print("Starting downstream (tffetcher & tfbsfetcher) stages for each group...")
    run_downstream(full_config)
    print("All groups processed.")


if __name__ == "__main__":
    config_path = Path(__file__).parent / "configs" / "example.yaml"
    main(config_path)
