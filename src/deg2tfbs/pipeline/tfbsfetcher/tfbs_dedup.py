"""
--------------------------------------------------------------------------------
<deg2tfbs project>
tfbs_dedup.py

This module contains functions to perform Jaccard-similarity based deduplication 
of TFBS entries (per transcription factor) and to generate a network-based visualization 
of the clusters in the second subplot. It is designed to be used as part of the 
TFBS Fetcher pipeline.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import itertools
from pathlib import Path
from typing import Set, Tuple, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

# Set a global random seed
np.random.seed(42)

def get_kmers(seq: str, k: int) -> Set[str]:
    """Return the set of all k-mers in the sequence."""
    seq = seq.upper()
    if len(seq) < k:
        return set()
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}

def deduplicate_tfbs(df: pd.DataFrame, config: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Apply Jaccard-similarity based deduplication on the TFBS mapping DataFrame.
    For each transcription factor (tf), binding sites with pairwise Jaccard similarity above
    the threshold (from config) are grouped into clusters.
    
    The k-mer size is taken from the configuration.
    For clusters of size >1, a unique cluster_id is assigned and all members are marked as clustered.
    The centroid (the one with the highest average similarity) is retained; non-centroid members
    are recorded in a separate DataFrame.

    Returns:
      - df_filtered: deduplicated mapping DataFrame (only centroids or singletons).
      - df_removed: DataFrame of filtered-out (non-centroid) TFBSs.
    """
    config_k = config.get("kmer_size")
    threshold = config.get("jaccard_similarity_threshold")
    if config_k is None or threshold is None:
        raise ValueError("Jaccard filtering enabled but 'kmer_size' and/or 'jaccard_similarity_threshold' is not specified in config.")

    filtered_rows = []
    removed_rows = []

    for tf, group in df.groupby("tf"):
        group = group.copy()
        k = config_k  
        group["kmers"] = group["tfbs"].apply(lambda s: get_kmers(s, k))
        
        cluster_counter = 1
        if len(group) <= 1:
            # Singletons
            group["is_centroid"] = True
            group["clustered"] = False
            group["cluster_id"] = None
            filtered_rows.append(group)
            continue

        # Build clusters via greedy set-lumping
        indices = group.index.tolist()
        clusters = []
        used = set()
        for idx in indices:
            if idx in used:
                continue
            cluster = {idx}
            for jdx in indices:
                if jdx == idx or jdx in used:
                    continue
                s1 = group.at[idx, "kmers"]
                s2 = group.at[jdx, "kmers"]
                union = s1 | s2
                sim = len(s1 & s2) / len(union) if union else 0
                if sim >= threshold:
                    cluster.add(jdx)
            clusters.append(cluster)
            used.update(cluster)

        # Mark centroids vs. non-centroids
        for cluster in clusters:
            if len(cluster) == 1:
                idx = list(cluster)[0]
                group.loc[idx, "is_centroid"] = True
                group.loc[idx, "clustered"] = False
                group.loc[idx, "cluster_id"] = None
                filtered_rows.append(group.loc[[idx]])
            else:
                cid = f"{tf}_cluster_{cluster_counter}"
                cluster_counter += 1
                cluster_list = list(cluster)
                best_idx = None
                best_avg_sim = -1
                # pick centroid: highest average similarity to others
                for idx in cluster_list:
                    sims = []
                    for jdx in cluster_list:
                        if idx == jdx:
                            continue
                        s1 = group.at[idx, "kmers"]
                        s2 = group.at[jdx, "kmers"]
                        union = s1 | s2
                        sim = len(s1 & s2) / len(union) if union else 0
                        sims.append(sim)
                    avg_sim = np.mean(sims) if sims else 0
                    if avg_sim > best_avg_sim:
                        best_avg_sim = avg_sim
                        best_idx = idx

                group.loc[best_idx, "is_centroid"] = True
                group.loc[best_idx, "clustered"] = True
                group.loc[best_idx, "cluster_id"] = cid
                filtered_rows.append(group.loc[[best_idx]])

                for idx in cluster_list:
                    if idx != best_idx:
                        row = group.loc[idx].copy()
                        row["is_centroid"] = False
                        row["clustered"] = True
                        row["cluster_id"] = cid
                        row["reason"] = "jaccard_similarity_clustered_non_centroid"
                        removed_rows.append(row)

    df_filtered = pd.concat(filtered_rows, ignore_index=True) if filtered_rows else pd.DataFrame()
    df_removed = pd.DataFrame(removed_rows) if removed_rows else pd.DataFrame()

    # Drop the 'kmers' column from both DataFrames
    for df_ in (df_filtered, df_removed):
        if "kmers" in df_.columns:
            df_.drop(columns=["kmers"], inplace=True)
    return df_filtered, df_removed

def plot_jaccard_summary(non_clustered: pd.DataFrame, clustered: pd.DataFrame, config: dict, out_path: Path):
    """
    Create a 3-panel summary figure:
    
      Panel 0 (Swarm):
        - Title: "TFBS Pairwise Similarity Scores"
        - x-axis: "Jaccard Similarity"
        - y-axis: "Transcription Factor" 
          with TFs ordered by increasing proportion of high-similarity pairs then reversed so that the TF with most high-similarity pairs appears at the top
        - Legend: "Lower Similarity (< threshold)" / "Higher Similarity (≥ threshold)"
      
      Panel 1 (Network Graph):
        - Title: "TFBS Similarity Network for Selected Transcription Factor"
        - x-axis label: "Nodes = TFBS sequences, Edges = Jaccard similarity (≥ threshold)"
          (this annotation appears right below the network plot)
        - Uses spring_layout with edge weights for positioning.
        - Node sizes: centroid = 80, non-centroid = 50.
        - Edge thicknesses are normalized modestly (light gray, alpha ≈ 0.5).
        - All axis splines (frames) are removed.
      
      Panel 2 (Cluster Snapshot):
        - Title: "Largest Cluster Snapshot for Selected Transcription Factor"
        - In the snapshot, the centroid sequence is annotated with " <-- RETAINED"
    
    Font sizes for titles, axis labels, ticks, and annotations have been increased.
    The transcription factor for Panels 1 and 2 is determined by the config key "selected_tf".
    If not provided or not found, a random TF is selected.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    import networkx as nx

    # Define increased font sizes
    title_fs = 14
    label_fs = 12
    tick_fs = 10
    annotation_fs = 12

    threshold = config.get("jaccard_similarity_threshold", 0.8)
    k = config.get("kmer_size", 6)
    fig_width = config.get("plot_width", 18)
    fig_height = config.get("plot_height", 5)

    # Combine data for pairwise comparisons
    df_raw = pd.concat([non_clustered, clustered], ignore_index=True)

    # -------------------------------
    # Panel 0: Swarm Plot with reversed y ticks
    # -------------------------------
    jaccard_data = []
    for tf, group in df_raw.groupby("tf"):
        seqs = group["tfbs"].tolist()
        if len(seqs) < 2:
            continue
        for i in range(len(seqs)):
            kmers_i = get_kmers(seqs[i], k)
            for j in range(i+1, len(seqs)):
                kmers_j = get_kmers(seqs[j], k)
                union = kmers_i | kmers_j
                sim = len(kmers_i & kmers_j) / len(union) if union else 0
                jaccard_data.append({"tf": tf, "similarity": sim})
    df_jaccard = pd.DataFrame(jaccard_data)
    if not df_jaccard.empty:
        df_prop = (
            df_jaccard
            .groupby("tf")["similarity"]
            .apply(lambda x: (x >= threshold).sum() / len(x))
            .reset_index(name="prop_above")
        )
        # Sort by increasing proportion then reverse the order so that TF with most high-similarity pairs is at the top
        df_prop = df_prop.sort_values("prop_above", ascending=True)
        tf_order = df_prop["tf"].tolist()
    else:
        df_prop = pd.DataFrame(columns=["tf", "prop_above"])
        tf_order = []

    def label_func(x):
        return "Higher Similarity (≥ threshold)" if x >= threshold else "Lower Similarity (< threshold)"
    df_jaccard["group"] = df_jaccard["similarity"].apply(label_func)
    color_map = {
        "Lower Similarity (< threshold)": "#D3D3D3",
        "Higher Similarity (≥ threshold)": "#F4A6A6"
    }

    fig, axes = plt.subplots(1, 3, figsize=(fig_width, fig_height))

    ax0 = axes[0]
    if df_prop.empty:
        ax0.text(0.5, 0.5, "No pairwise data", ha="center", va="center", transform=ax0.transAxes)
        ax0.set_axis_off()
    else:
        existing_tfs = df_jaccard["tf"].unique()
        order_filtered = [x for x in tf_order if x in existing_tfs]
        sns.stripplot(
            data=df_jaccard,
            x="similarity",
            y="tf",
            order=order_filtered,
            hue="group",
            palette=color_map,
            ax=ax0,
            size=4,
            jitter=True
        )
        ax0.axvline(threshold, color="lightgray", linestyle="--", linewidth=1.5)
        ax0.set_title("TFBS Pairwise Similarity Scores", fontsize=title_fs)
        ax0.set_xlabel("Jaccard Similarity", fontsize=label_fs)
        ax0.set_ylabel("Transcription Factor", fontsize=label_fs)
        ax0.tick_params(axis='both', labelsize=tick_fs)
        leg = ax0.legend(
            title="Pairwise Sequence Similarity",
            frameon=True,
            facecolor="white",
            edgecolor="black",
            loc="lower right",
            fontsize=tick_fs,
            title_fontsize=label_fs
        )
        ax0.spines["top"].set_visible(False)
        ax0.spines["right"].set_visible(False)

    # -------------------------------
    # Panel 1: Network Graph with User-Defined TF Selection
    # -------------------------------
    ax1 = axes[1]
    ax1.set_title("TFBS Similarity Network for Selected TF", fontsize=title_fs)
    # Set the x-axis label as the annotation (appears below the network plot)
    ax1.set_xlabel("Nodes = TFBS sequences, Edges = Jaccard similarity (≥ threshold)", fontsize=label_fs)
    
    if df_jaccard.empty:
        df_count = pd.DataFrame(columns=["tf", "count_above"])
    else:
        df_count = (
            df_jaccard
            .groupby("tf")["similarity"]
            .apply(lambda x: (x >= threshold).sum())
            .reset_index(name="count_above")
        )
    valid_tfs = clustered["tf"].unique() if not clustered.empty else []
    df_count_valid = df_count[df_count["tf"].isin(valid_tfs)]
    if df_count_valid.empty:
        ax1.text(0.5, 0.5, "No TF with clustered edges found", ha="center", va="center", transform=ax1.transAxes)
        ax1.set_axis_off()
        selected_tf = None
    else:
        selected_tf = config.get("selected_tf", None)
        if selected_tf is None or selected_tf not in df_count_valid["tf"].tolist():
            selected_tf = np.random.choice(df_count_valid["tf"].tolist())
        sub_tf = clustered[clustered["tf"] == selected_tf].copy()
        if sub_tf.empty:
            ax1.text(0.5, 0.5, f"No data for TF '{selected_tf}' after clustering", ha="center", va="center", transform=ax1.transAxes)
            ax1.set_axis_off()
        else:
            G = nx.Graph()
            for i, row in sub_tf.iterrows():
                G.add_node(row["tfbs"], is_centroid=row["is_centroid"])
            seqs = sub_tf["tfbs"].tolist()

            def scale_linewidth(sim):
                clipped = max(threshold, min(1.0, sim))
                frac = (clipped - threshold) / (1.0 - threshold)
                return 0.3 + 1.0 * frac

            kmers_cache = {s: get_kmers(s, k) for s in seqs}
            for i in range(len(seqs)):
                for j in range(i+1, len(seqs)):
                    s_i = seqs[i]
                    s_j = seqs[j]
                    kmers_i = kmers_cache[s_i]
                    kmers_j = kmers_cache[s_j]
                    union = kmers_i | kmers_j
                    sim = len(kmers_i & kmers_j) / len(union) if union else 0
                    if sim >= threshold:
                        G.add_edge(s_i, s_j, weight=sim)
            # Use spring_layout with edge weights for positioning
            pos = nx.spring_layout(G, seed=42, k=0.35, weight="weight")
            centroids = [n for n in G.nodes if G.nodes[n].get("is_centroid", False)]
            noncents = list(set(G.nodes) - set(centroids))
            for (u, v, d) in G.edges(data=True):
                lw = scale_linewidth(d.get('weight', 1.0))
                ax1.plot([pos[u][0], pos[v][0]],
                         [pos[u][1], pos[v][1]],
                         color="#CCCCCC", linewidth=lw, alpha=0.5)
            nx.draw_networkx_nodes(
                G, pos, nodelist=noncents,
                node_color="#F4A6A6", node_size=50, alpha=0.8,
                ax=ax1
            )
            nx.draw_networkx_nodes(
                G, pos, nodelist=centroids,
                node_color="#CD5C5C", node_size=80, alpha=0.9,
                ax=ax1
            )
            ax1.set_xticks([])
            ax1.set_yticks([])
            # Remove all spines (axis frames)
            for spine in ax1.spines.values():
                spine.set_visible(False)
            components = list(nx.connected_components(G))
            num_clusters = len(components)
            an_text = (f"TF: {selected_tf}\n"
                       f"Nodes: {len(G.nodes)}\n"
                       f"Clusters: {num_clusters}")
            ax1.text(
                0.95, 0.05, an_text, transform=ax1.transAxes,
                ha="right", va="bottom", fontsize=annotation_fs,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.7)
            )

    # -------------------------------
    # Panel 2: Cluster Snapshot from Selected TF
    # -------------------------------
    ax2 = axes[2]
    ax2.set_title("Representative Cluster for Selected TF", fontsize=title_fs)
    ax2.set_axis_off()
    if selected_tf is None:
        cluster_text = "No clusters found."
    else:
        sub_tf = clustered[clustered["tf"] == selected_tf].copy()
        sub_tf = sub_tf[sub_tf["cluster_id"].notnull()]
        if sub_tf.empty:
            cluster_text = "No clusters found for selected TF."
        else:
            c_sizes = sub_tf.groupby("cluster_id").size().reset_index(name="count")
            largest_cluster = c_sizes.sort_values("count", ascending=False).iloc[0]["cluster_id"]
            example_grp = sub_tf[sub_tf["cluster_id"] == largest_cluster]
            lines = [f"TF: {selected_tf} | Cluster: {largest_cluster}"]
            cent = example_grp[example_grp["is_centroid"] == True]
            noncent = example_grp[example_grp["is_centroid"] == False]
            if not cent.empty:
                c_seq = cent.iloc[0]["tfbs"]
                m_len = max(len(s) for s in example_grp["tfbs"])
                lines.append(f"{c_seq.ljust(m_len)}  <-- RETAINED")
                for _, rowx in noncent.iterrows():
                    seqxx = rowx["tfbs"]
                    lines.append(seqxx.ljust(m_len))
            else:
                lines.append("No centroid found?!")
            cluster_text = "\n".join(lines)
    ax2.text(
        0.05, 0.95, cluster_text,
        transform=ax2.transAxes,
        fontsize=annotation_fs, fontfamily="monospace",
        verticalalignment="top"
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

def export_cluster_strips(df_clusters: pd.DataFrame, out_path: Path):
    """
    Export a text-based summary of clusters.
    For each TF and each cluster, the centroid and its non-centroid sequences 
    are printed in an aligned, left-justified format.
    """
    if "type" not in df_clusters.columns:
        if "is_centroid" in df_clusters.columns:
            df_clusters = df_clusters.copy()
            df_clusters["type"] = df_clusters["is_centroid"].apply(lambda x: "centroid" if x else "non-centroid")
        else:
            df_clusters["type"] = "non-centroid"

    with open(out_path, "w") as f:
        for tf, group in df_clusters.groupby("tf"):
            f.write(f"TF: {tf}\n")
            for cid, subgrp in group.groupby("cluster_id"):
                f.write(f"  Cluster ID: {cid}\n")
                seqs = subgrp["tfbs"].tolist()
                max_len = max(len(s) for s in seqs)
                centroid = subgrp[subgrp["type"] == "centroid"]
                noncentroid = subgrp[subgrp["type"] == "non-centroid"]
                if not centroid.empty:
                    cent_seq = centroid.iloc[0]["tfbs"]
                    f.write(f"    {cent_seq.ljust(max_len)}  <-- RETAINED\n")
                for _, row in noncentroid.iterrows():
                    seq = row["tfbs"]
                    f.write(f"    {seq.ljust(max_len)}\n")
            f.write("\n")
