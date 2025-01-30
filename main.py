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

from deg2tfbs.pipeline.degfetcher.ceroni import run_ceroni_pipeline
from deg2tfbs.pipeline.degfetcher.mori import run_mori_pipeline
from deg2tfbs.pipeline.degfetcher.wu import run_wu_pipeline
from deg2tfbs.pipeline.degfetcher.zhu import run_zhu_pipeline
from deg2tfbs.pipeline.degfetcher.emani import run_emani_pipeline
from deg2tfbs.pipeline.degfetcher.schmidt import run_schmidt_pipeline
from deg2tfbs.pipeline.degfetcher.radzikowski import run_radzikowski_pipeline
from deg2tfbs.pipeline.degfetcher.bie import run_bie_pipeline
from deg2tfbs.pipeline.degfetcher.deter import run_deter_pipeline
from deg2tfbs.pipeline.degfetcher.jovanovic import run_jovanovic_pipeline
from deg2tfbs.pipeline.degfetcher.rajacharya import run_rajacharya_pipeline
from deg2tfbs.pipeline.degfetcher.durfee import run_durfee_pipeline
from deg2tfbs.pipeline.degfetcher.gummesson import run_gummesson_pipeline
from deg2tfbs.pipeline.degfetcher.houser import run_houser_pipeline
from deg2tfbs.pipeline.degfetcher.lu import run_lu_pipeline
from deg2tfbs.pipeline.degfetcher.sanchez_vasquez import run_sanchez_vasquez_pipeline
from deg2tfbs.pipeline.degfetcher.vazulka import run_vazulka_pipeline
from deg2tfbs.pipeline.degfetcher.kim import run_kim_pipeline
from deg2tfbs.pipeline.degfetcher.zhang import run_zhang_pipeline


def main(config_file: Path):
    with open(config_file, "r") as f:
        full_config = yaml.safe_load(f)

    # Run whichever pipelines are specified or we want to do in sequence
    run_ceroni_pipeline(full_config)
    run_mori_pipeline(full_config)
    run_wu_pipeline(full_config)
    run_zhu_pipeline(full_config)
    run_emani_pipeline(full_config)
    run_schmidt_pipeline(full_config)
    run_radzikowski_pipeline(full_config)
    run_bie_pipeline(full_config)
    run_deter_pipeline(full_config)
    run_jovanovic_pipeline(full_config)
    run_rajacharya_pipeline(full_config)
    run_durfee_pipeline(full_config)
    run_gummesson_pipeline(full_config)
    run_houser_pipeline(full_config)
    run_lu_pipeline(full_config)
    run_sanchez_vasquez_pipeline(full_config)
    run_vazulka_pipeline(full_config)
    run_kim_pipeline(full_config)
    run_zhang_pipeline(full_config)
    
    print("All pipelines completed.")

if __name__ == "__main__":
    # Example usage:
    config_path = Path(__file__).parent / "configs" / "example.yaml"
    main(config_path)
