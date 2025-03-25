"""
--------------------------------------------------------------------------------
<deg2tfbs project>
src/deg2tfbs/analysis/main.py

Loads TFBS mapping files and, based on the configuration mode, performs either:
  • Multi-comparison analysis:
      - Builds TF rosters
      - Computes pairings and UMAP clustering
      - Generates a UMAP plot, a TF roster heatmap, and TFBS counts plots
      - Merges tf2tfbs_mapping.csv files by cluster
  • Single-batch analysis:
      - Uses a single reference set to generate a TFBS counts plot
      - Generates a TFBS length density plot and writes a summary CSV

Outputs are organized into subdirectories based on the mode (e.g. "all_regs", "intersects_removed", or "single_batch").

Module Author(s): Eric J. South, Dunlop Lab
--------------------------------------------------------------------------------
"""

import sys
import logging
from pathlib import Path
import pandas as pd

from deg2tfbs.pipeline.utils import load_config, resolve_path
from deg2tfbs.analysis import data, roster, compare, tf_cluster_merger
from deg2tfbs.analysis import plot_tf_roster_clustering, plot_tf_roster_heatmap, plot_tfbs_counts
from deg2tfbs.analysis import plot_tfbs_length_density

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def run_multi_batch_analysis(rosters, regs_reference, total_regs, output_root: Path, run_label: str, 
                             reference_set: str, title_suffix: str, config):
    """
    Run multi-comparison analyses (pairing comparisons, clustering, heatmap, and TFBS counts)
    using the provided rosters. Output is saved under output_root/run_label.
    """
    out_dir = output_root / run_label
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = out_dir / "data"
    csv_dir.mkdir(exist_ok=True)
    
    comparisons = config["analysis"].get("comparisons", [])
    # Process pairings.
    pairing_summary_text = ""
    for pair in comparisons:
        if not (isinstance(pair, (list, tuple)) and len(pair) == 2):
            logger.warning(f"Invalid pairing {pair}; skipping.")
            continue
        group1, group2 = pair
        if group1 not in rosters or group2 not in rosters:
            logger.warning(f"Pairing {pair} not found in the rosters; skipping.")
            continue
        try:
            vec1 = rosters[group1]
            vec2 = rosters[group2]
            details = roster.compute_pairing_details(vec1, vec2)
            df_pair = compare.generate_pairing_dataframe(regs_reference, details)
            pair_name = f"{group1}_vs_{group2}"
            csv_file = csv_dir / f"{pair_name}_detailed.csv"
            compare.save_pairing_dataframe(df_pair, csv_file)
            summary = compare.print_pairing_summary(pair_name, details, total_regs)
            pairing_summary_text += summary + "\n"
        except Exception as e:
            logger.error(f"Error processing pairing {pair}: {e}")
    
    (csv_dir / "pairing_summaries.txt").write_text(pairing_summary_text)
    
    # Run UMAP clustering.
    selected_keys = set()
    for pair in comparisons:
        if isinstance(pair, (list, tuple)) and len(pair) == 2:
            selected_keys.update(pair)
    selected_keys.add(reference_set)
    missing_keys = selected_keys - set(rosters.keys())
    if missing_keys:
        logger.warning(f"The following keys from comparisons are missing: {missing_keys}")
    selected_keys = selected_keys.intersection(rosters.keys())
    if not selected_keys:
        logger.error("No valid keys found from comparisons to include in the UMAP.")
        sys.exit(1)
    
    selected_rosters = {key: rosters[key] for key in selected_keys}
    matrix_raw, sample_names_raw = roster.build_matrix_from_rosters(selected_rosters)
    
    # UMAP parameters (set in __main__)
    adata_raw = plot_tf_roster_clustering.run_scanpy_clustering(matrix_raw, sample_names_raw,
                                                                 n_neighbors, min_dist)
    # Add additional info.
    reg_count_raw = matrix_raw.sum(axis=1)
    adata_raw.obs["reg_count"] = reg_count_raw
    adata_raw.obs["sample"] = list(adata_raw.obs["sample"])  # ensure strings
    
    # Generate UMAP plot.
    umap_title = f"UMAP Clustering of TF Rosters {title_suffix}"
    umap_output_file = out_dir / "umap_clustering.png"
    plot_tf_roster_clustering.plot_umap(adata_raw, umap_title,
                                          str(umap_output_file),
                                          figsize=umap_figsize, alpha=umap_alpha, dpi=dpi,
                                          reference_set=reference_set)
    
    # Generate heatmap.
    cluster_mapping = dict(zip(adata_raw.obs["sample"], adata_raw.obs["leiden"]))
    heatmap_output_file = out_dir / "tf_roster_heatmap.png"
    
    group_defs = config["pipeline"]["stages"]["tffetcher"]["input"]["deg_csv_groups"]
    group_labels = { key: group_def.get("plot_name", key)
                     for key, group_def in group_defs.items() if isinstance(group_def, dict) }
    
    plot_tf_roster_heatmap.plot_heatmap(
        rosters=rosters,
        regs_reference=regs_reference,
        cluster_mapping=cluster_mapping,
        output_file=str(heatmap_output_file),
        reference_set=reference_set,
        group_labels=group_labels
    )
    
    # TFBS Counts plots.
    mapping_files = data.get_tfbs_mapping_files(config["pipeline"]["stages"]["tfbsfetcher"]["root_dir"])
    for key, roster_csv in mapping_files.items():
        try:
            mapping_file = roster_csv.parent / "tf2tfbs_mapping.csv"
            if not mapping_file.exists():
                logger.error(f"Mapping file {mapping_file} not found for {key}. Skipping TFBS counts plot.")
                continue
            out_file = out_dir / f"{key}_tfbs_counts.png"
            # If intersections are excluded, compute the set of regulators removed.
            exclude_set = None
            if exclude_intersections:
                orig = all_rosters[key]
                mod = rosters_to_use[key]
                exclude_set = {reg for reg, o, m in zip(regs_reference, orig, mod) if o == 1 and m == 0}
            plot_tfbs_counts.plot_tfbs_counts(
                source_csv=roster_csv,
                mapping_csv=mapping_file,
                output_file=str(out_file),
                exclude_regulators=exclude_set,
                group_labels=group_labels
            )
        except Exception as e:
            logger.error(f"Error generating TFBS counts plot for {key}: {e}")
    
    return cluster_mapping


def run_single_batch_analysis(config, repo_root):
    """
    Run single-batch analysis:
      - Loads mapping for the reference_set.
      - Errors if the reference_set is missing.
      - Creates output directories under "single_batch".
      - Generates the TFBS counts plot and the TFBS length density plot with summary CSV.
    """
    analysis_config = config.get("analysis", {})
    tfbs_dir = resolve_path(analysis_config.get("tfbsfetcher_dir", "pipeline/tfbsfetcher"), repo_root)
    mapping = data.get_tfbs_mapping_files(tfbs_dir)
    reference_set = analysis_config.get("reference_set")
    
    if not reference_set or reference_set not in mapping:
        logger.error("reference_set must be defined and present in tfbsfetcher output when multi_comparison is false.")
        sys.exit(1)
    
    plot_base_dir = resolve_path(analysis_config.get("plot_output_dir", "analysis/plots"), repo_root)
    csv_base_dir = resolve_path(analysis_config.get("csv_output_dir", "analysis/data"), repo_root)
    output_plots_dir = plot_base_dir / "single_batch"
    output_csv_dir = csv_base_dir / "single_batch"
    output_plots_dir.mkdir(parents=True, exist_ok=True)
    output_csv_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate TFBS counts plot.
    tfbs_counts_csv = mapping[reference_set]
    counts_output_file = output_plots_dir / "tfbs_counts.png"
    try:
        mapping_file = tfbs_counts_csv.parent / "tf2tfbs_mapping.csv"
        if not mapping_file.exists():
            logger.error(f"Mapping file {mapping_file} not found for {reference_set}.")
            sys.exit(1)
        group_defs = config["pipeline"]["stages"]["tffetcher"]["input"]["deg_csv_groups"]
        group_labels = { key: group_def.get("plot_name", key)
                         for key, group_def in group_defs.items() if isinstance(group_def, dict) }
        plot_tfbs_counts.plot_tfbs_counts(
            source_csv=tfbs_counts_csv,
            mapping_csv=mapping_file,
            output_file=str(counts_output_file),
            group_labels=group_labels
        )
    except Exception as e:
        logger.error(f"Error generating tfbs_counts plot: {e}")
        sys.exit(1)
    
    # Generate TFBS length density plot and summary CSV.
    length_plot_output = output_plots_dir / "tfbs_length_density.png"
    summary_csv_output = output_csv_dir / "tfbs_length_summary.csv"
    try:
        plot_tfbs_length_density.plot_and_save_tfbs_length_analysis(
            mapping_csv=mapping[reference_set],
            plot_output=length_plot_output,
            summary_output=summary_csv_output
        )
    except Exception as e:
        logger.error(f"Error generating tfbs_length_density plot and summary: {e}")
        sys.exit(1)
    
    logger.info("Single-batch analysis complete.")


# --------------------------
# Main Execution
# --------------------------
if __name__ == "__main__":
    # Determine repository root from this script's location.
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]
    logger.info(f"Repository root: {repo_root}")
    
    # Load configuration from the default config file.
    config_path = resolve_path("configs/example.yaml", repo_root)
    config = load_config(str(config_path))
    analysis_config = config.get("analysis", {})
    
    # Common preparatory steps: get reference_set and mapping.
    reference_set = analysis_config.get("reference_set")
    if not reference_set:
        logger.error("No reference_set defined in the configuration.")
        sys.exit(1)
    
    tfbs_dir = resolve_path(analysis_config.get("tfbsfetcher_dir", "pipeline/tfbsfetcher"), repo_root)
    mapping = data.get_tfbs_mapping_files(tfbs_dir)
    if not mapping:
        logger.error("No tfbs mapping CSV files found.")
        sys.exit(1)
    
    include_unassigned = analysis_config.get("include_unassigned", False)
    comparisons = analysis_config.get("comparisons", [])
    if not include_unassigned:
        valid_keys = set([reference_set])
        for pair in comparisons:
            if isinstance(pair, (list, tuple)) and len(pair) == 2:
                valid_keys.update(pair)
        mapping = {k: v for k, v in mapping.items() if k in valid_keys}
    
    if reference_set not in mapping:
        logger.error(f"Reference set {reference_set} not found in mapping keys.")
        sys.exit(1)
    
    # Branch early based on the mode.
    multi_comparison = analysis_config.get("multi_comparison", True)
    
    if multi_comparison:
        # Multi-comparison mode: perform UMAP and roster computations.
        n_neighbors = int(analysis_config.get("umap_n_neighbors", 15))
        min_dist = float(analysis_config.get("umap_min_dist", 0.1))
        raw_figsize = analysis_config.get("umap_figsize", [8, 6])
        if isinstance(raw_figsize, str):
            umap_figsize = tuple(float(x.strip()) for x in raw_figsize.split(","))
        else:
            umap_figsize = tuple(float(x) for x in raw_figsize)
        umap_alpha = float(analysis_config.get("umap_alpha", 0.5))
        dpi = int(analysis_config.get("umap_dpi", 150))
        
        # Build TF rosters.
        regs_reference, ref_vector = roster.build_reference_roster(mapping, reference_set)
        total_regs = len(regs_reference)
        logger.info(f"Reference set '{reference_set}' has {total_regs} unique regulators.")
        
        all_rosters = {}
        for key, csv_file in mapping.items():
            try:
                all_rosters[key] = roster.create_boolean_vector(regs_reference, csv_file)
            except Exception as e:
                logger.error(f"Error building roster for {key}: {e}")
                sys.exit(1)
        
        exclude_intersections = analysis_config.get("exclude_intersections", False)
        if exclude_intersections:
            rosters_to_use = roster.exclude_intersections(all_rosters)
            title_suffix = "(Intersections Removed)"
            run_label = "intersects_removed"
        else:
            rosters_to_use = all_rosters
            title_suffix = "(All Regulators)"
            run_label = "all_regs"
        
        final_tf_sets = { source: { reg for reg, val in zip(regs_reference, vec) if val == 1 }
                          for source, vec in rosters_to_use.items() }
        
        base_plot_output_dir = resolve_path(analysis_config.get("plot_output_dir", "analysis/plots"), repo_root)
        base_plot_output_dir.mkdir(parents=True, exist_ok=True)
        
        cluster_mapping = run_multi_batch_analysis(rosters_to_use, regs_reference, total_regs, base_plot_output_dir,
                                                    run_label, reference_set, title_suffix, config)
        
        merged_outputs_dir = resolve_path("analysis/outputs", repo_root)
        tf_cluster_merger.merge_tfbs_by_cluster(mapping, cluster_mapping, merged_outputs_dir, run_label, final_tf_sets)
        
        logger.info("TFBS batch analysis, plotting, and merged CSV creation complete.")
    else:
        # Single-batch mode: minimal processing.
        run_single_batch_analysis(config, repo_root)
        sys.exit(0)
