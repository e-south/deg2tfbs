"""
--------------------------------------------------------------------------------
<deg2tfbs project>
src/deg2tfbs/analysis/main.py

Loads TFBS mapping files, builds TF rosters, computes pairings,
performs clustering, and then generates three types of plots:
  - UMAP clustering plot (plot_tf_roster_clustering)
  - TF roster heatmap (plot_tf_roster_heatmap)
  - TFBS counts per source (plot_tfbs_counts)

If the configuration key "exclude_intersections" is true, the rosters are
computed after removing intersections, and all downstream analyses (including
Leiden clustering) use these modified rosters. The results are saved in a subdirectory
named "intersects_removed". Otherwise, results are saved in a subdirectory named "all_regs".
  
After the analysis and plotting, the script also merges the original
tf2tfbs_mapping.csv files by cluster (using the same clustering from UMAP/Leiden)
and writes merged, deduplicated CSVs into analysis/outputs.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import sys
import argparse
import logging
from pathlib import Path

import pandas as pd

from deg2tfbs.pipeline.utils import load_config, resolve_path
from deg2tfbs.analysis import data, roster, compare, tf_cluster_merger
from deg2tfbs.analysis import plot_tf_roster_clustering, plot_tf_roster_heatmap, plot_tfbs_counts

logger = logging.getLogger(__name__)

def run_analysis(rosters, regs_reference, total_regs, output_root: Path, run_label: str, 
                 reference_set: str, title_suffix: str):
    """
    Run downstream analyses (pairing comparisons, clustering, heatmap, and TFBS counts)
    using the provided rosters dictionary. The output is saved under output_root/run_label.
    
    The title of the UMAP plot is augmented with title_suffix to indicate whether
    intersections were removed.

    Returns:
      cluster_mapping: A dict mapping each source (roster) to its Leiden cluster label.
    """
    out_dir = output_root / run_label
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = out_dir / "data"
    csv_dir.mkdir(exist_ok=True)
    
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
    adata_raw = plot_tf_roster_clustering.run_scanpy_clustering(matrix_raw, sample_names_raw,
                                                                 n_neighbors=n_neighbors,
                                                                 min_dist=min_dist)
    # Add additional information.
    reg_count_raw = matrix_raw.sum(axis=1)
    adata_raw.obs["reg_count"] = reg_count_raw
    adata_raw.obs["sample"] = list(adata_raw.obs["sample"])  # ensure strings
    
    # Append a descriptive suffix to the UMAP title.
    umap_title = f"UMAP Clustering of TF Rosters {title_suffix}"
    umap_output_file = out_dir / "umap_clustering.png"
    plot_tf_roster_clustering.plot_umap(adata_raw, umap_title,
                                          str(umap_output_file),
                                          figsize=umap_figsize, alpha=umap_alpha, dpi=dpi,
                                          reference_set=reference_set)
    
    # Heatmap.
    cluster_mapping = dict(zip(adata_raw.obs["sample"], adata_raw.obs["leiden"]))
    heatmap_output_file = out_dir / "tf_roster_heatmap.png"
    plot_tf_roster_heatmap.plot_heatmap(
        rosters=rosters,
        regs_reference=regs_reference,
        cluster_mapping=cluster_mapping,
        output_file=str(heatmap_output_file),
        reference_set=reference_set
    )
    
    # TFBS Counts plots.
    for key, roster_csv in mapping.items():
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
                exclude_set = {reg for reg, o, m in zip(regs_reference, orig, mod) if o==1 and m==0}
            plot_tfbs_counts.plot_tfbs_counts(
                source_csv=roster_csv,
                mapping_csv=mapping_file,
                output_file=str(out_file),
                exclude_regulators=exclude_set
            )
        except Exception as e:
            logger.error(f"Error generating TFBS counts plot for {key}: {e}")
    
    return cluster_mapping

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TFBS Batch Analysis and Plotting")
    parser.add_argument("--config", type=str, default="configs/example.yaml",
                        help="Path to the global configuration YAML file (relative to repo root)")
    args = parser.parse_args()
    
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[1]
    logger.info(f"Repository root: {repo_root}")
    
    config_path = resolve_path(args.config, repo_root)
    config = load_config(str(config_path))
    analysis_config = config.get("analysis", {})
    
    # Directories.
    tfbs_dir = resolve_path(analysis_config.get("tfbsfetcher_dir", "pipeline/tfbsfetcher"), repo_root)
    base_plot_output_dir = resolve_path(analysis_config.get("plot_output_dir", "analysis/plots"), repo_root)
    base_plot_output_dir.mkdir(parents=True, exist_ok=True)
    
    # UMAP parameters.
    n_neighbors = int(analysis_config.get("umap_n_neighbors", 15))
    min_dist = float(analysis_config.get("umap_min_dist", 0.1))
    raw_figsize = analysis_config.get("umap_figsize", [8, 6])
    if isinstance(raw_figsize, str):
        umap_figsize = tuple(float(x.strip()) for x in raw_figsize.split(","))
    else:
        umap_figsize = tuple(float(x) for x in raw_figsize)
    umap_alpha = float(analysis_config.get("umap_alpha", 0.5))
    dpi = int(analysis_config.get("umap_dpi", 150))
    
    reference_set = analysis_config.get("reference_set")
    comparisons = analysis_config.get("comparisons", [])
    include_unassigned = analysis_config.get("include_unassigned", False)
    exclude_intersections = analysis_config.get("exclude_intersections", False)
    
    if not reference_set:
        logger.error("No reference_set defined in the configuration.")
        sys.exit(1)
    
    logger.info(f"TFBS directory: {tfbs_dir}")
    logger.info(f"Base plot output directory: {base_plot_output_dir}")
    
    # -----------------------
    # Get mapping.
    # -----------------------
    mapping = data.get_tfbs_mapping_files(tfbs_dir)
    if not mapping:
        logger.error("No tfbs mapping CSV files found.")
        sys.exit(1)
    if not include_unassigned:
        valid_keys = set([reference_set])
        for pair in comparisons:
            if isinstance(pair, (list, tuple)) and len(pair) == 2:
                valid_keys.update(pair)
        mapping = {k: v for k, v in mapping.items() if k in valid_keys}
    
    if reference_set not in mapping:
        logger.error(f"Reference set {reference_set} not found in mapping keys.")
        sys.exit(1)
    
    # -----------------------
    # Build TF rosters.
    # -----------------------
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
    
    # If exclude_intersections is true, perform the removal BEFORE clustering.
    if exclude_intersections:
        rosters_to_use = roster.exclude_intersections(all_rosters)
        title_suffix = "(Intersections Removed)"
        run_label = "intersects_removed"
    else:
        rosters_to_use = all_rosters
        title_suffix = "(All Regulators)"
        run_label = "all_regs"
    
    # Compute final TF sets based on the filtered binary vectors.
    final_tf_sets = { source: { reg for reg, val in zip(regs_reference, vec) if val == 1 }
                      for source, vec in rosters_to_use.items() }
    
    # Run the analysis (which returns the cluster mapping)
    cluster_mapping = run_analysis(rosters_to_use, regs_reference, total_regs, base_plot_output_dir, 
                                   run_label, reference_set, title_suffix)
    
    # Merge tf2tfbs_mapping.csv files by cluster using the merger module.
    # Pass the final_tf_sets so that only rows corresponding to a final '1' are retained.
    merged_outputs_dir = resolve_path("analysis/outputs", repo_root)
    tf_cluster_merger.merge_tfbs_by_cluster(mapping, cluster_mapping, merged_outputs_dir, run_label, final_tf_sets)
    
    logger.info("TFBS batch analysis, plotting, and merged CSV creation complete.")
