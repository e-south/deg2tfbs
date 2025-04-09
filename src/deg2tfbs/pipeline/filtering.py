"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/filtering.py

Used by tffetcher.py, accessing downstream parsers from the tfbsfetcher directory,
to check if any transcription factors do not have any binding site data, in which
case they are flagged and not carried forward to the next stage.

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import logging
import yaml
from pathlib import Path

from deg2tfbs.pipeline.tfbsfetcher.parsers import get_tfbs_parser

logger = logging.getLogger(__name__)

def filter_mapping_rows_by_tfbs(rows):
    """
    Given a list of regulator mapping rows, load the TFBS configuration from YAML,
    load all enabled TFBS sources, and filter out mapping rows for regulators
    that do not have any binding site sequences.
    
    Parameters:
      rows (List[Dict]): List of mapping rows (each a dict with at least a 'regulator' key).
      
    Returns:
      List[Dict]: The filtered mapping rows.
    """
    project_root = Path(__file__).parent.parent
    full_config_path = project_root / "configs" / "example.yaml"
    # Load the YAML file into a dictionary
    with open(full_config_path, "r") as f:
        full_config = yaml.safe_load(f)

    tfbsfetcher_config = full_config["pipeline"]["stages"].get("tfbsfetcher", {})
    bs_sources = tfbsfetcher_config.get("sources", {}).get("binding_sites", {})
    if not bs_sources:
        logger.warning("[filtering] No binding_sites sources configured; returning original rows.")
        return rows

    binding_sites = {}
    
    # Iterate over each configured binding site source.
    for source_key, source_conf in bs_sources.items():
        # Accept both boolean or dict type configurations.
        enabled = source_conf if isinstance(source_conf, bool) else source_conf.get("enabled", True)
        if not enabled:
            logger.info(f"[filtering] Skipping disabled TFBS source: {source_key}")
            continue
        try:
            parser = get_tfbs_parser(source_key)
            bs_data = parser.parse()  # Expected: {tf: set(binding_site_seq)}
            # For each TF, update the binding site set filtering out empty strings.
            for tf, sites in bs_data.items():
                tf_norm = tf.lower().strip()
                # Filter out entries that are just empty strings (after stripping).
                non_empty_sites = {site for site in sites if site.strip()}
                if non_empty_sites:  # Only update if there is at least one non-empty site.
                    if tf_norm not in binding_sites:
                        binding_sites[tf_norm] = set()
                    binding_sites[tf_norm].update(non_empty_sites)
        except Exception as e:
            logger.error(f"[filtering] Error parsing TFBS source '{source_key}': {e}")
            continue

    # Build the set of regulators having at least one binding site.
    regulators_with_bs = {tf for tf, sites in binding_sites.items() if sites}
    logger.info(f"[filtering] Regulators with binding sites: {', '.join(sorted(regulators_with_bs))}")

    # Filter out rows whose regulator (normalized) is not in the binding site data.
    filtered_rows = [row for row in rows if row.get("regulator", "").lower().strip() in regulators_with_bs]
    
    num_removed = len(rows) - len(filtered_rows)
    if num_removed:
        logger.info(f"[filtering] Filtered out {num_removed} mapping row(s) without any binding site evidence.")
    else:
        logger.info("[filtering] All mapping rows have at least one binding site.")
    
    return filtered_rows