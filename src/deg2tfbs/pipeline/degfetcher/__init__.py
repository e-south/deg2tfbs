"""
--------------------------------------------------------------------------------
<deg2tfbs project>
deg2tfbs/pipeline/degfetcher/__init__.py

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

from typing import Dict, Callable
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
from deg2tfbs.pipeline.degfetcher.treitz import run_treitz_pipeline
from deg2tfbs.pipeline.degfetcher.treitz_schmidt_concordant import run_treitz_schmidt_concordant_pipeline

# A dictionary mapping each key in your config to the actual function
MODULE_MAP: Dict[str, Callable] = {
    "ceroni": run_ceroni_pipeline,
    "mori": run_mori_pipeline,
    "wu": run_wu_pipeline,
    "zhu": run_zhu_pipeline,
    "emani": run_emani_pipeline,
    "schmidt": run_schmidt_pipeline,
    "radzikowski": run_radzikowski_pipeline,
    "bie": run_bie_pipeline,
    "deter": run_deter_pipeline,
    "jovanovic": run_jovanovic_pipeline,
    "rajacharya": run_rajacharya_pipeline,
    "durfee": run_durfee_pipeline,
    "gummesson": run_gummesson_pipeline,
    "houser": run_houser_pipeline,
    "lu": run_lu_pipeline,
    "sanchez_vasquez": run_sanchez_vasquez_pipeline,
    "vazulka": run_vazulka_pipeline,
    "kim": run_kim_pipeline,
    "zhang": run_zhang_pipeline,
    "treitz": run_treitz_pipeline,
    "treitz_schmidt_concordant": run_treitz_schmidt_concordant_pipeline
}

def run_pipeline_for(module_name: str, full_config: dict) -> None:
    """
    Dynamically run the pipeline for the specified module_name
    using the appropriate function from MODULE_MAP.
    """
    if module_name not in MODULE_MAP:
        raise ValueError(f"No pipeline function found for module: {module_name}")
    # Call the mapped function, passing in the config
    MODULE_MAP[module_name](full_config)
