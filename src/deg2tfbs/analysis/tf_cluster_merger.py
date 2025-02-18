"""
--------------------------------------------------------------------------------
<deg2tfbs project>
src/deg2tfbs/analysis/tf_cluster_merger.py

This module merges tf2tfbs_mapping.csv files using Leidon cluster results from the 
main analysis script—which compares and groups TF rosters derived from all tfbsfetcher 
output files. 

Specifically, the module groups source mapping CSV files by cluster (as determined 
by the TF roster clustering), concatenates them, and removes duplicate rows based on 
the 'tf' and 'tfbs' columns. For each cluster, the merged CSV file is saved in a 
subdirectory along with a text file that lists the contributing sources.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import logging
from pathlib import Path
import pandas as pd

logger = logging.getLogger(__name__)

def merge_tfbs_by_cluster(mapping: dict, cluster_mapping: dict, outputs_dir: Path, 
                          run_label: str, final_tf_sets: dict = None) -> None:
    """
    For each cluster, group all tf2tfbs_mapping CSV files (from the provided mapping dict)
    corresponding to sources in that cluster, concatenate them, remove duplicate rows based
    on the 'tf' and 'tfbs' columns, and, if intersections were removed (run_label=="intersects_removed")
    and final_tf_sets is provided, filter out rows whose 'tf' is not in the final set for that source.
    Also, write a sources.txt file listing which source keys contributed.
    
    The merged outputs are saved under a subdirectory outputs_dir/<run_label>/cluster_<id>.
    
    Preconditions:
      - mapping: dict mapping source keys (str) to CSV file Paths.
      - cluster_mapping: dict mapping source keys (str) to cluster labels (str).
      - outputs_dir: a Path where outputs can be created.
      - run_label: a string indicating the run configuration ("intersects_removed" or "all_regs").
      - final_tf_sets: optional dict mapping source keys to a set of TFs that remain (if intersections removed).
    Postconditions:
      - For each cluster, a subdirectory outputs_dir/<run_label>/cluster_<id> is created with:
          * tf2tfbs_mapping.csv (merged, deduplicated, and filtered CSV)
          * sources.txt (list of contributing sources)
    """
    # Create a run-specific subdirectory.
    run_dir = outputs_dir / run_label
    run_dir.mkdir(parents=True, exist_ok=True)
    
    clusters = {}
    for source, cluster in cluster_mapping.items():
        clusters.setdefault(cluster, []).append(source)
    
    for cluster, sources in clusters.items():
        merged_df = pd.DataFrame()
        # Process each source.
        for src in sources:
            if src in mapping:
                try:
                    df = pd.read_csv(mapping[src])
                    # If final_tf_sets is provided, filter the rows to only those TFs in the final set.
                    if final_tf_sets and src in final_tf_sets:
                        before_rows = df.shape[0]
                        df = df[df['tf'].isin(final_tf_sets[src])]
                        after_rows = df.shape[0]
                        logger.info(f"Source '{src}' in cluster {cluster}: Filtered {before_rows - after_rows} rows not in final roster.")
                    merged_df = pd.concat([merged_df, df], ignore_index=True)
                except Exception as e:
                    logger.error(f"Error reading mapping file for source '{src}': {e}")
            else:
                logger.warning(f"Source '{src}' not found in mapping; skipping.")
        if merged_df.empty:
            logger.warning(f"No data merged for cluster {cluster}.")
            continue
        
        # If run_label indicates intersections were removed, optionally log the overall common TFs.
        if run_label == "intersects_removed" and final_tf_sets and len(sources) > 1:
            tf_sets = [final_tf_sets[src] for src in sources if src in final_tf_sets]
            if tf_sets:
                common_tfs = set.intersection(*tf_sets)
                logger.info(f"Cluster {cluster}: Common TFs across final rosters (if any) are: {common_tfs}")
        
        # Remove duplicates based on 'tf' and 'tfbs'
        merged_df = merged_df.drop_duplicates(subset=['tf', 'tfbs'])
        cluster_dir = run_dir / f"cluster_{cluster}"
        cluster_dir.mkdir(exist_ok=True)
        output_csv = cluster_dir / "tf2tfbs_mapping.csv"
        try:
            merged_df.to_csv(output_csv, index=False)
            logger.info(f"Merged tfbs mapping for cluster {cluster} saved to {output_csv}")
        except Exception as e:
            logger.error(f"Error saving merged CSV for cluster {cluster}: {e}")
        # Write contributing sources.
        sources_file = cluster_dir / "sources.txt"
        try:
            with sources_file.open("w") as f:
                for src in sources:
                    f.write(src + "\n")
            logger.info(f"Sources for cluster {cluster} saved to {sources_file}")
        except Exception as e:
            logger.error(f"Error saving sources for cluster {cluster}: {e}")
