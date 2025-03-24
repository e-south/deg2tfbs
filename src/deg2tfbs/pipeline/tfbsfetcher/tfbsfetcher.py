"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tfbsfetcher.py

Maps TFs to TFBSs by processing one or more binding site resources (e.g., EcoCyc,
RegulonDB) as configured in the YAML file. The module loads a TF mapping CSV (from 
tffetcher) that contains at least the columns:
    [gene, regulator, deg_source]
Optionally it also contains:
    [polarity, is_global_regulator, is_sigma_factor]

For each TF (i.e. regulator), the module uses one or more parsers to retrieve the 
corresponding binding site sequences. It then joins the TF mapping with the TFBS 
information. The output CSV contains (for example):
    tf, tfbs, tfbs_source, gene, deg_source, polarity, is_sigma_factor, is_global_regulator

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import logging
from pathlib import Path
from collections import defaultdict
from typing import Dict, Set, List

import pandas as pd

from deg2tfbs.pipeline.tfbsfetcher.parsers import get_tfbs_parser

# Configure the module-level logger with a simple formatter.
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.propagate = False  # Avoid duplicate logs.
if not logger.handlers:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(handler)
else:
    for h in logger.handlers:
        h.setFormatter(logging.Formatter("%(message)s"))

def aggregate_deg_source(series: pd.Series) -> str:
    seen = []
    for s in series.dropna():
        s = s.strip()
        if not s:
            continue
        tokens = s.split("-")
        for token in tokens:
            token = token.strip()
            if token and token not in seen:
                seen.append(token)
    return "-".join(seen)

def run_tfbsfetcher_stage(config: dict):
    """
    Main entry point for the TFBS-fetcher stage.
    
    Steps:
      1. Validate configuration and determine output directories.
      2. Load the TF mapping CSV from the tffetcher stage.
      3. If the configuration option 'use_only_enriched_tfs' is True,
          load the TF enrichment summary and filter the mapping to include only
          enriched TFs that:
             - Appear in the enrichment 'regulator' column,
             - Do NOT have "yes" in the 'is_sigma_factor' column,
             - Do NOT have "yes" in the 'is_nucleoid_regulator' column.
          Then, sort these rows by the 'fdr' column in ascending order and select
          the top N rows based on the 'top_n' parameter in config.
      4. Load and parse each enabled TFBS dataset.
      5. Merge the TF mapping with the parsed TFBS data.
      6. Deduplicate and write the final TFBS mapping CSV.
    """
    # Validate configuration keys.
    if "root_dir" not in config or "batch_id" not in config:
        raise ValueError("tfbsfetcher config must have 'root_dir' and 'batch_id' keys")

    root_dir = config["root_dir"]
    batch_id = config["batch_id"]
    project_root = Path(__file__).parent.parent.parent
    tfbs_root_dir = (project_root / root_dir).resolve()
    tfbs_out_dir = tfbs_root_dir / batch_id
    tfbs_out_dir.mkdir(parents=True, exist_ok=True)

    # Load input TF mapping CSV from the tffetcher stage.
    input_conf = config.get("input", {})
    tf_batch_id = input_conf.get("tf_batch_id", None)
    if not tf_batch_id:
        raise ValueError("[TFBSFetcher] No 'tf_batch_id' specified in config under 'input'.")

    # Look for the mapping CSV inside the "csvs" subdirectory.
    tf_mapping_csv = (
        tfbs_out_dir.parent.parent / "tffetcher" / tf_batch_id / "csvs" / "deg2tf_mapping.csv"
    )
    if not tf_mapping_csv.exists():
        raise FileNotFoundError(f"[TFBSFetcher] Cannot find TF mapping CSV: {tf_mapping_csv.relative_to(project_root)}")

    logger.info(f"[TFBSFetcher] Reading TF mapping from {tf_mapping_csv.relative_to(project_root)}")
    df_tfmap = pd.read_csv(tf_mapping_csv)

    # Verify required columns.
    required_cols = {"gene", "regulator", "deg_source"}
    missing = required_cols - set(df_tfmap.columns)
    if missing:
        raise ValueError(f"[TFBSFetcher] deg2tf mapping is missing columns: {missing}")

    # Standardize the TF names and deg_source.
    df_tfmap["regulator"] = df_tfmap["regulator"].str.lower().str.strip()
    df_tfmap["deg_source"] = df_tfmap["deg_source"].astype(str).str.lower().str.strip()

    # Check for user-defined option "use_only_enriched_tfs"
    use_only_enriched = config.get("use_only_enriched_tfs", False)
    enriched_tfs = set()
    if use_only_enriched:
        # Assume the enrichment summary is located in the same tfbatch as tffetcher output.
        tf_enrichment_csv = (
            tfbs_out_dir.parent.parent / "tffetcher" / tf_batch_id / "csvs" / "tf_enrichment_summary.csv"
        )
        if not tf_enrichment_csv.exists():
            raise FileNotFoundError(f"[TFBSFetcher] Cannot find TF enrichment summary CSV: {tf_enrichment_csv.relative_to(project_root)}")
        df_enrich = pd.read_csv(tf_enrichment_csv)
        
        # Standardize relevant columns.
        df_enrich["regulator"] = df_enrich["regulator"].str.lower().str.strip()
        df_enrich["is_sigma_factor"] = df_enrich["is_sigma_factor"].str.lower().str.strip()
        df_enrich["is_nucleoid_regulator"] = df_enrich["is_nucleoid_regulator"].str.lower().str.strip()
        
        # Apply filtering: retain rows where neither 'is_sigma_factor' nor 'is_nucleoid_regulator' is "yes".
        df_enrich_filtered = df_enrich[
            (df_enrich["is_sigma_factor"] != "yes") &
            (df_enrich["is_nucleoid_regulator"] != "yes")
        ]
        
        # Sort the filtered enrichment data by the 'fdr' column (ascending).
        df_enrich_filtered = df_enrich_filtered.sort_values(by="fdr", ascending=True)
        
        # Get the top N rows if specified in config.
        top_n = config.get("params", {}).get("top_n", None)
        if top_n is not None:
            df_enrich_filtered = df_enrich_filtered.head(top_n)
        
        # Create a set of enriched TFs from the filtered enrichment summary.
        enriched_tfs = set(df_enrich_filtered["regulator"])
        
        original_count = len(df_tfmap)
        df_tfmap = df_tfmap[df_tfmap["regulator"].isin(enriched_tfs)]
        logger.info(f"[TFBSFetcher] use_only_enriched_tfs enabled: Filtered mapping from {original_count} to {len(df_tfmap)} rows based on enriched TFs (top {top_n} by fdr).")

    # Load each TFBS dataset as configured.
    sources_conf = config.get("sources", {}).get("binding_sites", {})
    if not sources_conf:
        logger.info("[TFBSFetcher] No TFBS binding_sites sources configured, skipping.")
        return

    tfbs_sources = []
    for key, conf in sources_conf.items():
        enabled = conf.get("enabled", True) if isinstance(conf, dict) else conf
        if not enabled:
            logger.info(f"[TFBSFetcher] Skipping TFBS dataset: {key} (not enabled).")
            continue
        parser_label = conf.get("parser", key)
        try:
            parser = get_tfbs_parser(key)
            tfbs_data = parser.parse()
            tfbs_sources.append((parser_label, tfbs_data))
            logger.info(f"[TFBSFetcher] Parser {parser_label} => found {len(tfbs_data)} unique TFs")
        except Exception as e:
            logger.error(f"[TFBSFetcher] Failed to parse {parser_label}: {str(e)}")
            raise

    # Merge TFBS data and track which source provided each binding site.
    source_tracker = defaultdict(lambda: defaultdict(set))
    empty_tfbs_count = 0
    empty_tfbs_by_source = defaultdict(int)
    for source_label, tfbs_data in tfbs_sources:
        for tf, sites in tfbs_data.items():
            for site in sites:
                if not site:
                    empty_tfbs_count += 1
                    empty_tfbs_by_source[source_label] += 1
                    continue
                source_tracker[tf][site].add(source_label)
    if empty_tfbs_count:
        breakdown = ", ".join([f"{src}: {cnt}" for src, cnt in empty_tfbs_by_source.items()])
        logger.info(f"[TFBSFetcher] Skipped {empty_tfbs_count} empty TFBS entries during merge (breakdown by source: {breakdown}).")
    
    # If filtering by enriched TFs, warn if any enriched TF does not have TFBS information.
    if use_only_enriched:
        missing_tfs = enriched_tfs - set(source_tracker.keys())
        if missing_tfs:
            logger.warning("[TFBSFetcher] Warning: The following enriched TFs do not have any TFBS information:")
            for tf in sorted(missing_tfs):
                logger.warning(f" - {tf}")

    # Join the TF mapping with TFBS data.
    rows = []
    for _, row in df_tfmap.iterrows():
        tf = row["regulator"]
        if not tf:
            raise AssertionError("[TFBSFetcher] Encountered an empty TF value in the mapping CSV.")
        if tf not in source_tracker:
            continue
        for site_seq, sources in source_tracker[tf].items():
            if not site_seq:
                logger.error(f"[TFBSFetcher] Encountered an empty TFBS for TF '{tf}'. This should be caught upstream.")
                continue
            row_data = {
                "tf": tf,
                "tfbs": site_seq,
                "gene": row["gene"],
                "deg_source": row["deg_source"],
                "polarity": row.get("polarity", "NA"),
                "tfbs_source": "_".join(sorted(sources))
            }
            # Propagate optional columns if present.
            if "is_sigma_factor" in df_tfmap.columns:
                row_data["is_sigma_factor"] = row["is_sigma_factor"]
            if "is_global_regulator" in df_tfmap.columns:
                row_data["is_global_regulator"] = row["is_global_regulator"]
            rows.append(row_data)

    if not rows:
        logger.warning("[TFBSFetcher] No matching TFBS entries found for any TF in mapping file.")
        return

    df_out = pd.DataFrame(rows)
    expected_columns = {"tf", "tfbs", "gene", "deg_source", "tfbs_source"}
    missing_output = expected_columns - set(df_out.columns)
    if missing_output:
        raise AssertionError(f"[TFBSFetcher] Missing columns in output: {missing_output}")

    # Final assertions: make sure no empty values exist in the tf or tfbs columns.
    if df_out["tf"].str.strip().eq("").any():
        raise AssertionError("[TFBSFetcher] Found empty values in 'tf' column.")
    if df_out["tfbs"].str.strip().eq("").any():
        raise AssertionError("[TFBSFetcher] Found empty values in 'tfbs' column.")

    # Deduplicate rows.
    # Uniqueness is defined by the tuple: (tf, tfbs)
    group_cols = ["tf", "tfbs"]
    duplicate_mask = df_out.duplicated(subset=group_cols, keep=False)
    num_duplicates = duplicate_mask.sum()
    if num_duplicates:
        logger.info(f"[TFBSFetcher] Found {num_duplicates} duplicate rows; aggregating duplicates based on {group_cols}")
        logger.warning("[TFBSFetcher] Note: The duplication check is stringent and does not account for cases where binding sites differ by ±1 nucleotide at either end.")
        aggregation_rules = {
            "gene": lambda x: "_".join(sorted(set(x))),
            "deg_source": aggregate_deg_source,
            "polarity": lambda x: x.mode().iloc[0] if not x.mode().empty else "NA",
            "tfbs_source": lambda x: "_".join(sorted(set(x))),
            "is_sigma_factor": "first",
            "is_global_regulator": "first"
        }
        agg_rules = {col: rule for col, rule in aggregation_rules.items() if col not in group_cols and col in df_out.columns}
        df_clean = df_out.groupby(group_cols, as_index=False).agg(agg_rules)
    else:
        df_clean = df_out.drop_duplicates(subset=group_cols)

    final_duplicates = df_clean.duplicated(subset=group_cols).sum()
    logger.info(f"[TFBSFetcher] Final dataset contains {len(df_clean)} rows with {final_duplicates} duplicates")
    if final_duplicates:
        raise AssertionError(f"[TFBSFetcher] Post-processing duplicates detected: {final_duplicates}")

    # Log summary statistics.
    num_unique_tfs = df_clean["tf"].nunique()
    logger.info(f"Unique TFs found: {num_unique_tfs}")
    
    # Format TFBS counts as a grid with 5 columns.
    tfbs_counts = df_clean.groupby("tf")["tfbs"].nunique()
    max_tf_len = max(len(tf) for tf in tfbs_counts.index)
    max_count_len = max(len(str(count)) for count in tfbs_counts.values)
    formatted_entries = [f"{tf:<{max_tf_len}} : {count:>{max_count_len}}" for tf, count in tfbs_counts.items()]
    # Group the entries in chunks of 5, preserving a grid-like structure.
    for i in range(0, len(formatted_entries), 5):
        chunk = "    ".join(formatted_entries[i:i+5])
        logger.info(chunk)

    multi_source = df_clean[df_clean["tfbs_source"].str.contains("_")]
    logger.info(f"Number of (tf, tfbs) pairs found in multiple sources: {len(multi_source)}")
    if not multi_source.empty:
        multi_source_breakdown = multi_source["tfbs_source"].value_counts()
        breakdown_str = " | ".join([f"{k}: {v}" for k, v in multi_source_breakdown.items()])
        logger.info(f"Multi-source pairs breakdown: {breakdown_str}")
        logger.info("Sample of multi-source pairs:\n" + multi_source.head().to_string())

    out_csv = tfbs_out_dir / "tf2tfbs_mapping.csv"
    # Log relative path for output.
    logger.info(f"[TFBSFetcher] Saved {len(df_clean)} unique TFBS mappings to {out_csv.relative_to(project_root)}")
    df_clean.to_csv(out_csv, index=False)

def merge_tfbs_dicts(tfbs_list: List[Dict[str, Set[str]]]) -> Dict[str, Set[str]]:
    """
    Union-based merge: If a TF appears in multiple dictionaries, return the union
    of its binding site sets.
    """
    merged = {}
    for tdict in tfbs_list:
        for tf, site_set in tdict.items():
            if tf not in merged:
                merged[tf] = set()
            merged[tf].update(site_set)
    return merged
