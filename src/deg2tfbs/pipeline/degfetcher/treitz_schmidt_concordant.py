"""
--------------------------------------------------------------------------------
<deg2tfbs project>
treitz_schmidt_concordant.py

Module to merge mRNA fold change data from Treitz et al. with protein fold change 
data from Schmidt et al. for the M9-acetate versus M9-glucose comparison. The module 
merges the two datasets based on gene, categorizes genes (concordant, RNA specific, 
Protein specific, or unchanged) based on a threshold, creates a concordant plot 
(with Pearson correlation annotation and top 5 annotations for both up and down),
and saves separate CSV files for concordant up-regulated and concordant down-regulated genes.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from scipy.stats import pearsonr
import seaborn as sns

def run_treitz_schmidt_concordant_pipeline(full_config: dict):
    """
    Loads the Treitz fold change data and the Schmidt Acetate_vs_Glucose tidy data,
    merges them on gene (with robust matching), creates a concordant plot with point 
    sizes scaled by the logarithm of 'avg_exp', annotates the top five points for both 
    up and down concordance, and outputs CSV files for concordant up- and down-regulated genes.
    
    Assumes that:
      - Treitz output CSV: "treitz_fold_change_data.csv" (columns: gene, ratio, log2_fc)
      - Schmidt tidy output CSV: "schmidt_acetate_vs_glucose_tidy.csv" (columns including 
        gene, log2FC, avg_exp, etc.)
    
    The comparative condition is "Acetate_vs_Glucose".
    """
    # Build paths from pipeline configuration.
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config['pipeline']['stages']['degfetcher']['root_dir']
    batch_id = full_config['pipeline']['stages']['degfetcher']['batch_id']
    batch_dir = output_root / batch_id
    # Assume both modules use the same csv subdirectory.
    csv_dir = batch_dir / full_config["datasets"]["treitz"]["output"]["csv_subdir"]
    
    # Load Treitz fold change data.
    treitz_file = csv_dir / "treitz_fold_change_data.csv"
    df_treitz = pd.read_csv(treitz_file)
    # Normalize gene names robustly and rename log2_fc column.
    df_treitz["gene"] = df_treitz["gene"].astype(str).str.lower().str.strip()
    df_treitz = df_treitz.rename(columns={"log2_fc": "log2FC_mrna"})
    # Drop any existing 'color' column.
    df_treitz = df_treitz.drop(columns=["color"], errors="ignore")
    
    # Load Schmidt tidy data for Acetate_vs_Glucose.
    schmidt_file = csv_dir / "schmidt_acetate_vs_glucose_tidy.csv"
    df_schmidt = pd.read_csv(schmidt_file)
    # Normalize gene names and rename log2FC column.
    df_schmidt["gene"] = df_schmidt["gene"].astype(str).str.lower().str.strip()
    df_schmidt = df_schmidt.rename(columns={"log2FC": "log2FC_protein"})
    # Drop any existing 'color' column.
    df_schmidt = df_schmidt.drop(columns=["color"], errors="ignore")
    
    # Merge on gene.
    merged = pd.merge(df_treitz, df_schmidt, on="gene", how="inner")
    
    # Save genes that did not match for debugging purposes.
    # Genes in Treitz not found in the merged data.
    non_matched_treitz = df_treitz[~df_treitz["gene"].isin(merged["gene"])].copy()
    non_matched_treitz["dataset"] = "Treitz"
    # Genes in Schmidt not found in the merged data.
    non_matched_schmidt = df_schmidt[~df_schmidt["gene"].isin(merged["gene"])].copy()
    non_matched_schmidt["dataset"] = "Schmidt"
    # Combine the non-matching genes.
    non_matches = pd.concat([non_matched_treitz, non_matched_schmidt])
    non_matches = non_matches.drop_duplicates(subset=["gene", "dataset"])
    # Save the non-matching genes to CSV.
    non_matches.to_csv(csv_dir / "treitz_schmidt_no_matches.csv", index=False)
    
    # Report matching statistics.
    genes_treitz = set(df_treitz["gene"])
    genes_schmidt = set(df_schmidt["gene"])
    union_genes = genes_treitz.union(genes_schmidt)
    merged_count = len(merged)
    
    # Use the Treitz threshold for concordance.
    threshold = full_config["datasets"]["treitz"]["thresholds"]["log2_fc_threshold"]
    
    # Categorize genes.
    def assign_category(row, t):
        r = row["log2FC_mrna"]
        p = row["log2FC_protein"]
        if (abs(r) >= t) and (abs(p) >= t) and (np.sign(r) == np.sign(p)):
            return "concordant"
        elif abs(r) >= t and abs(p) < t:
            return "RNA specific"
        elif abs(p) >= t and abs(r) < t:
            return "Protein specific"
        else:
            return "unchanged"
    
    merged["category"] = merged.apply(lambda row: assign_category(row, threshold), axis=1)
    
    # Compute Pearson correlation for concordant genes.
    concordant = merged[merged["category"] == "concordant"]
    if len(concordant) > 1:
        r_val, _ = pearsonr(concordant["log2FC_mrna"], concordant["log2FC_protein"])
    else:
        r_val = np.nan

    # Compute concordant up- and down-regulated subsets.
    concordant_up = merged[(merged["category"] == "concordant") &
                           (merged["log2FC_mrna"] > 0) &
                           (merged["log2FC_protein"] > 0)]
    concordant_down = merged[(merged["category"] == "concordant") &
                             (merged["log2FC_mrna"] < 0) &
                             (merged["log2FC_protein"] < 0)]
    
    # Create the concordant plot.
    plt.figure(figsize=(7,6))
    cat_colors = {
        "concordant": "purple",
        "RNA specific": "blue",
        "Protein specific": "red",
        "unchanged": "gray"
    }
    # For each category, scale marker sizes using the log of avg_exp (from Schmidt data).
    for cat, color in cat_colors.items():
        subset = merged[merged["category"] == cat]
        sizes = np.log10(subset["avg_exp"] + 1) * 20  # Adjust scaling factor as needed.
        plt.scatter(
            subset["log2FC_mrna"], subset["log2FC_protein"],
            color=color, alpha=0.35, edgecolors="none", s=sizes, label=cat
        )
    
    # Annotate top five up-regulated concordant points.
    top5_up = concordant_up.sort_values(by="log2FC_mrna", ascending=False).head(5)
    for _, row in top5_up.iterrows():
        plt.annotate(row["gene"],
                     xy=(row["log2FC_mrna"], row["log2FC_protein"]),
                     xytext=(row["log2FC_mrna"]+0.1, row["log2FC_protein"]+0.1),
                     fontsize=9, color="black")
    
    # Annotate top five down-regulated concordant points.
    top5_down = concordant_down.sort_values(by="log2FC_mrna", ascending=True).head(5)
    for _, row in top5_down.iterrows():
        plt.annotate(row["gene"],
                     xy=(row["log2FC_mrna"], row["log2FC_protein"]),
                     xytext=(row["log2FC_mrna"]-0.1, row["log2FC_protein"]-0.1),
                     fontsize=9, color="black")
    
    plt.axvline(0, color="gray", linestyle="--")
    plt.axhline(0, color="gray", linestyle="--")
    plt.xlabel("Protein Log2 Fold Change (Treitz et al.)")
    plt.ylabel("Protein Log2 Fold Change (Schmidt et al.)")
    plt.title("Concordance of proteomic datasets\nM9-Acetate versus M9-Glucose")
    plt.legend(frameon=False)
    plt.annotate(f"Pearson r = {r_val:.2f}", xy=(0.05, 0.95), xycoords="axes fraction",
                 fontsize=12, color="black")
    # Remove top and right spines.
    sns.despine()
    
    # Save the plot.
    plot_dir = batch_dir / full_config["datasets"]["treitz"]["output"]["plot_subdir"]
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_file = plot_dir / "treitz_schmidt_acetate_vs_glucose_concordant.png"
    plt.savefig(plot_file, dpi=300)
    plt.close()

    # Remove any color column from merged before saving CSVs.
    merged = merged.drop(columns=["color"], errors="ignore")
    
    # Sort the output CSVs.
    concordant_up = concordant_up.sort_values(by=["log2FC_mrna", "avg_exp"], ascending=[False, False])
    concordant_down = concordant_down.sort_values(by=["log2FC_mrna", "avg_exp"], ascending=[True, False])
    
    merged.to_csv(csv_dir / "treitz_schmidt_merged.csv", index=False)
    concordant_up.to_csv(csv_dir / "treitz_schmidt_concordant_up.csv", index=False)
    concordant_down.to_csv(csv_dir / "treitz_schmidt_concordant_down.csv", index=False)

    union_count = len(union_genes)
    print(f"[Treitz-Schmidt Concordant] Completed. Merged {merged_count} genes "
          f"({merged_count/union_count:.1%} of union).")


if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)
    run_treitz_schmidt_concordant_pipeline(full_config)
