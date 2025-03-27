"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tfbsfetcher.py

Maps TFs to TFBSs by processing one or more binding site resources as configured 
in the YAML file.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import yaml
import logging
from pathlib import Path
from collections import defaultdict
import copy
import pandas as pd
import torch

from deg2tfbs.pipeline.tfbsfetcher.parsers import get_tfbs_parser
from deg2tfbs.pipeline.tfbsfetcher import tfbs_dedup

# Configure logger.
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.propagate = False
if not logger.handlers:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(handler)
else:
    for h in logger.handlers:
        h.setFormatter(logging.Formatter("%(message)s"))

def print_diagnostics(df: pd.DataFrame):
    """
    For each transcription factor in the final mapping, print counts of binding site pairs 
    (where one is a substring) that differ by exactly 1, exactly 2, and 3 or more nucleotides in length.
    """
    for tf, group in df.groupby("tf"):
        rows = group.to_dict(orient="records")
        diff1 = diff2 = diff3plus = 0
        for i in range(len(rows)):
            seq_i = rows[i]["tfbs"]
            for j in range(i+1, len(rows)):
                seq_j = rows[j]["tfbs"]
                if seq_i in seq_j or seq_j in seq_i:
                    diff = abs(len(seq_i) - len(seq_j))
                    if diff == 1:
                        diff1 += 1
                    elif diff == 2:
                        diff2 += 1
                    elif diff >= 3:
                        diff3plus += 1
        total = diff1 + diff2 + diff3plus
        logger.info(f"TF '{tf}': diff==1: {diff1}, diff==2: {diff2}, diff>=3: {diff3plus} (total pairs: {total})")

def run_tfbsfetcher_stage(config: dict):
    """
    Main entry point for the TFBS-fetcher stage.
    Steps:
      1. Validate configuration and create output directories.
      2. Load TF mapping CSV.
      3. Filter mapping by enriched TFs if specified.
      4. Load and parse each enabled TFBS dataset.
      5. Merge TF mapping with parsed TFBS data.
      6. Deduplicate and filter TFBS using Jaccard-similarity (using k-mer size from config).
      7. Write final mapping CSV and filtered-out CSV.
      8. Generate a three-panel Jaccard summary plot.
      9. Export cluster strips (alignment-like text summary of clusters).
      10. Print diagnostic metrics regarding length differences.
    """
    if "root_dir" not in config or "batch_id" not in config:
        raise ValueError("tfbsfetcher config must have 'root_dir' and 'batch_id' keys")

    root_dir = config["root_dir"]
    batch_id = config["batch_id"]
    project_root = Path(__file__).parent.parent.parent
    tfbs_root_dir = (project_root / root_dir).resolve()
    tfbs_out_dir = tfbs_root_dir / batch_id
    tfbs_out_dir.mkdir(parents=True, exist_ok=True)
    
    out_csv_dir = tfbs_out_dir / "csvs"
    out_plots_dir = tfbs_out_dir / "plots"
    out_csv_dir.mkdir(parents=True, exist_ok=True)
    out_plots_dir.mkdir(parents=True, exist_ok=True)

    # Load TF mapping CSV.
    input_conf = config.get("input", {})
    tf_batch_id = input_conf.get("tf_batch_id", None)
    if not tf_batch_id:
        raise ValueError("[TFBSFetcher] No 'tf_batch_id' specified in config under 'input'.")
    tf_mapping_csv = (tfbs_out_dir.parent.parent / "tffetcher" / tf_batch_id / "csvs" / "deg2tf_mapping.csv")
    if not tf_mapping_csv.exists():
        raise FileNotFoundError(f"[TFBSFetcher] Cannot find TF mapping CSV: {tf_mapping_csv.relative_to(project_root)}")
    logger.info(f"[TFBSFetcher] Reading TF mapping from {tf_mapping_csv.relative_to(project_root)}")
    df_tfmap = pd.read_csv(tf_mapping_csv)

    required_cols = {"gene", "regulator", "deg_source"}
    missing = required_cols - set(df_tfmap.columns)
    if missing:
        raise ValueError(f"[TFBSFetcher] deg2tf mapping is missing columns: {missing}")

    df_tfmap["regulator"] = df_tfmap["regulator"].str.lower().str.strip()
    df_tfmap["deg_source"] = df_tfmap["deg_source"].astype(str).str.lower().str.strip()

    use_only_enriched = config.get("use_only_enriched_tfs", False)
    enriched_tfs = set()
    if use_only_enriched:
        tf_enrichment_csv = (tfbs_out_dir.parent.parent / "tffetcher" / tf_batch_id / "csvs" / "tf_enrichment_summary.csv")
        if not tf_enrichment_csv.exists():
            raise FileNotFoundError(f"[TFBSFetcher] Cannot find TF enrichment summary CSV: {tf_enrichment_csv.relative_to(project_root)}")
        df_enrich = pd.read_csv(tf_enrichment_csv)
        df_enrich["regulator"] = df_enrich["regulator"].str.lower().str.strip()
        df_enrich["is_sigma_factor"] = df_enrich["is_sigma_factor"].str.lower().str.strip()
        df_enrich["is_nucleoid_regulator"] = df_enrich["is_nucleoid_regulator"].str.lower().str.strip()
        df_enrich_filtered = df_enrich[(df_enrich["is_sigma_factor"] != "yes") & (df_enrich["is_nucleoid_regulator"] != "yes")]
        df_enrich_filtered = df_enrich_filtered.sort_values(by="fdr", ascending=True)
        top_n = config.get("params", {}).get("top_n", None)
        if top_n is not None:
            df_enrich_filtered = df_enrich_filtered.head(top_n)
        enriched_tfs = set(df_enrich_filtered["regulator"])
        original_count = len(df_tfmap)
        df_tfmap = df_tfmap[df_tfmap["regulator"].isin(enriched_tfs)]
        logger.info(f"[TFBSFetcher] use_only_enriched_tfs enabled: Filtered mapping from {original_count} to {len(df_tfmap)} rows based on enriched TFs (top {top_n} by fdr).")

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

    if use_only_enriched:
        missing_tfs = enriched_tfs - set(source_tracker.keys())
        if missing_tfs:
            logger.warning("[TFBSFetcher] Warning: The following enriched TFs do not have any TFBS information:")
            for tf in sorted(missing_tfs):
                logger.warning(f" - {tf}")

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
    if df_out["tf"].str.strip().eq("").any():
        raise AssertionError("[TFBSFetcher] Found empty values in 'tf' column.")
    if df_out["tfbs"].str.strip().eq("").any():
        raise AssertionError("[TFBSFetcher] Found empty values in 'tfbs' column.")

    group_cols = ["tf", "tfbs"]
    duplicate_mask = df_out.duplicated(subset=group_cols, keep=False)
    num_duplicates = duplicate_mask.sum()
    if num_duplicates:
        logger.info(f"[TFBSFetcher] Found {num_duplicates} duplicate rows; aggregating duplicates based on {group_cols}")
        aggregation_rules = {
            "gene": lambda x: '_'.join(sorted(set(x))),
            "deg_source": lambda x: '-'.join(sorted(set(x))),
            "polarity": lambda x: x.mode().iloc[0] if not x.mode().empty else "NA",
            "tfbs_source": lambda x: '_'.join(sorted(set(x))),
            "is_sigma_factor": "first",
            "is_global_regulator": "first"
        }
        agg_rules = {col: rule for col, rule in aggregation_rules.items() if col not in group_cols and col in df_out.columns}
        df_out = df_out.groupby(group_cols, as_index=False).agg(agg_rules)
    else:
        df_out = df_out.drop_duplicates(subset=group_cols)

    if config.get("apply_jaccard_filter", False):
        logger.info("[TFBSFetcher] Applying Jaccard similarity filtering to reduce redundant TFBS entries...")
        df_before = df_out.copy()
        df_clean, df_removed = tfbs_dedup.deduplicate_tfbs(df_out, config)
        logger.info(f"[TFBSFetcher] Jaccard filtering reduced TFBS entries from {len(df_before)} to {len(df_clean)}.")
        tfbs_counts_before = df_before.groupby("tf")["tfbs"].nunique()
        tfbs_counts_after = df_clean.groupby("tf")["tfbs"].nunique()
        tfs = sorted(set(tfbs_counts_before.index) & set(tfbs_counts_after.index))
        logger.info("[TFBSFetcher] TF-level pruning summary (TF: before -> after [# filtered])")
        for tf in tfs:
            before = tfbs_counts_before.get(tf, 0)
            after = tfbs_counts_after.get(tf, 0)
            filtered = before - after
            perc = (filtered / before * 100) if before > 0 else 0
            logger.info(f"  {tf:<10} : {before} -> {after}  ({filtered} filtered, {perc:.1f}%)")
        filtered_out_csv = out_csv_dir / "tf2tfbs_filtered_out.csv"
        df_removed.to_csv(filtered_out_csv, index=False)
        logger.info(f"[TFBSFetcher] Saved filtered-out TFBS entries to {filtered_out_csv.relative_to(project_root)}")
        # Define two groups:
        clustered = pd.concat([df_clean[df_clean["cluster_id"].notnull()], df_removed[df_removed["cluster_id"].notnull()]])
        non_clustered = df_clean[df_clean["cluster_id"].isnull()]
    else:
        logger.info("[TFBSFetcher] Jaccard filtering not enabled; skipping redundancy pruning.")
        df_clean = df_out
        clustered = pd.DataFrame()
        non_clustered = df_out

    out_csv = out_csv_dir / "tf2tfbs_mapping.csv"
    df_clean.to_csv(out_csv, index=False)
    logger.info(f"[TFBSFetcher] Saved {len(df_clean)} unique TFBS mappings to {out_csv.relative_to(project_root)}")

    if config.get("apply_jaccard_filter", False):
        summary_plot_out = out_plots_dir / "tfbs_jaccard_summary.png"
        tfbs_dedup.plot_jaccard_summary(non_clustered, clustered, config, summary_plot_out)
        logger.info(f"[TFBSFetcher] Saved Jaccard summary plot to {summary_plot_out.relative_to(project_root)}")
    else:
        logger.info("[TFBSFetcher] Jaccard filtering not enabled; skipping Jaccard summary plot.")

    if config.get("apply_jaccard_filter", False):
        strips_out = out_csv_dir / "tfbs_cluster_strips.txt"
        tfbs_dedup.export_cluster_strips(clustered, strips_out)
        logger.info(f"[TFBSFetcher] Saved cluster strips to {strips_out.relative_to(project_root)}")

    print_diagnostics(df_clean)
    num_unique_tfs = df_clean["tf"].nunique()
    logger.info(f"[TFBSFetcher] Unique TFs in final dataset: {num_unique_tfs}")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with config_path.open("r") as f:
        full_config = yaml.safe_load(f)
    tfbs_config = full_config.get("pipeline", {}).get("stages", {}).get("tfbsfetcher", {})
    run_tfbsfetcher_stage(tfbs_config)
