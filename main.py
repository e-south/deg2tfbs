"""
--------------------------------------------------------------------------------
<deg2tfbs project>

main.py

Module Author(s): Eric J. South
Dunlop Lab

Coordinates the overall pipelines for ceroni, mori, wu, zhu, emani, schmidt, radzikowski, etc.
--------------------------------------------------------------------------------
"""
import yaml
from pathlib import Path

from deg2tfbs.pipeline.dataloader.ceroni import run_ceroni_pipeline
from deg2tfbs.pipeline.dataloader.mori import run_mori_pipeline
from deg2tfbs.pipeline.dataloader.wu import run_wu_pipeline
from deg2tfbs.pipeline.dataloader.zhu import run_zhu_pipeline
from deg2tfbs.pipeline.dataloader.emani import run_emani_pipeline
from deg2tfbs.pipeline.dataloader.schmidt import run_schmidt_pipeline
from deg2tfbs.pipeline.dataloader.radzikowski import run_radzikowski_pipeline
from deg2tfbs.pipeline.dataloader.bie import run_bie_pipeline
from deg2tfbs.pipeline.dataloader.deter import run_deter_pipeline
from deg2tfbs.pipeline.dataloader.jovanovic import run_jovanovic_pipeline
from deg2tfbs.pipeline.dataloader.rajacharya import run_rajacharya_pipeline

def main(config_file: Path):
    with open(config_file, "r") as f:
        full_config = yaml.safe_load(f)

    # Run whichever pipelines are specified or we want to do in sequence
    # run_ceroni_pipeline(full_config)
    # run_mori_pipeline(full_config)
    # run_wu_pipeline(full_config)
    # run_zhu_pipeline(full_config)
    # run_emani_pipeline(full_config)
    # run_schmidt_pipeline(full_config)
    # run_radzikowski_pipeline(full_config)
    # run_bie_pipeline(full_config)
    # run_deter_pipeline(full_config)
    # run_jovanovic_pipeline(full_config)
    run_rajacharya_pipeline(full_config)

    print("All pipelines completed.")

if __name__ == "__main__":
    # Example usage:
    config_path = Path(__file__).parent / "configs" / "example.yaml"
    main(config_path)
