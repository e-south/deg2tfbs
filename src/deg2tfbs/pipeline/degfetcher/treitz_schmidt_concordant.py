"""
--------------------------------------------------------------------------------
<deg2tfbs project>
treitz_schmidt_concordant.py

Module to merge proteomic fold change data from Treitz et al. with proteomic fold 
change data from Schmidt et al. for the M9-acetate versus M9-glucose comparison. 
The module merges the two datasets based on gene, categorizes genes (concordant, 
Treitz-specific, Schmidt-specific, or unchanged) based on a threshold, creates a 
concordant plot (with Pearson correlation annotation), and saves separate CSV 
files for concordant up-regulated and concordant down-regulated genes.

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
    Loads Treitz fold change data and the Schmidt Acetate_vs_Glucose tidy data,
    merges them on gene (with robust matching), creates a concordant plot with point 
    sizes scaled by the logarithm of 'avg_exp', annotates the top five points for both 
    up and down concordance, and outputs CSV files for concordant up- and down-regulated genes.
    
    Assumes that:
      - Treitz output CSV: "treitz_fold_change_data.csv" 
        (columns: gene, ratio, log2_fc)
      - Schmidt tidy output CSV: "schmidt_acetate_vs_glucose_tidy.csv" 
        (columns including gene, log2FC, avg_exp, etc.)

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
    # Normalize gene names robustly and rename the log2_fc column.
    df_treitz["gene"] = df_treitz["gene"].astype(str).str.lower().str.strip()
    df_treitz = df_treitz.rename(columns={"log2_fc": "log2FC_treitz"})
    # Drop any existing 'color' column.
    df_treitz = df_treitz.drop(columns=["color"], errors="ignore")
    
    # Load Schmidt tidy data for Acetate_vs_Glucose.
    schmidt_file = csv_dir / "schmidt_acetate_vs_glucose_tidy.csv"
    df_schmidt = pd.read_csv(schmidt_file)
    # Normalize gene names and rename the log2FC column.
    df_schmidt["gene"] = df_schmidt["gene"].astype(str).str.lower().str.strip()
    df_schmidt = df_schmidt.rename(columns={"log2FC": "log2FC_schmidt"})
    # Drop any existing 'color' column.
    df_schmidt = df_schmidt.drop(columns=["color"], errors="ignore")
    
    # Merge on gene.
    merged = pd.merge(df_treitz, df_schmidt, on="gene", how="inner")
    
    # Save genes that did not match for debugging purposes.
    non_matched_treitz = df_treitz[~df_treitz["gene"].isin(merged["gene"])].copy()
    non_matched_treitz["dataset"] = "Treitz"
    non_matched_schmidt = df_schmidt[~df_schmidt["gene"].isin(merged["gene"])].copy()
    non_matched_schmidt["dataset"] = "Schmidt"
    non_matches = pd.concat([non_matched_treitz, non_matched_schmidt])
    non_matches = non_matches.drop_duplicates(subset=["gene", "dataset"])
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
        r = row["log2FC_treitz"]
        p = row["log2FC_schmidt"]
        if (abs(r) >= t) and (abs(p) >= t) and (np.sign(r) == np.sign(p)):
            return "Concordant"
        elif abs(r) >= t and abs(p) < t:
            return "Treitz et al. specific"
        elif abs(p) >= t and abs(r) < t:
            return "Schmidt et al. specific"
        else:
            return "Unchanged"
    
    merged["category"] = merged.apply(lambda row: assign_category(row, threshold), axis=1)
    
    # Compute Pearson correlation for concordant genes only if there are >= 2.
    concordant = merged[merged["category"] == "Concordant"]
    if len(concordant) > 1:
        r_val, _ = pearsonr(concordant["log2FC_treitz"], concordant["log2FC_schmidt"])
    else:
        r_val = np.nan  # Fewer than 2 points => Pearson r is NaN by definition.
    
    # Compute concordant up- and down-regulated subsets.
    concordant_up = merged[
        (merged["category"] == "Concordant") 
        & (merged["log2FC_treitz"] > 0) 
        & (merged["log2FC_schmidt"] > 0)
    ]
    concordant_down = merged[
        (merged["category"] == "Concordant") 
        & (merged["log2FC_treitz"] < 0) 
        & (merged["log2FC_schmidt"] < 0)
    ]
    
    # Create the concordant plot.
    plt.figure(figsize=(7, 6))
    cat_colors = {
        "Concordant": "purple",
        "Treitz et al. specific": "blue",
        "Schmidt et al. specific": "red",
        "Unchanged": "gray"
    }
    # For each category, scale marker sizes using the log of avg_exp (Schmidt data).
    for cat, color in cat_colors.items():
        subset = merged[merged["category"] == cat]
        sizes = np.log10(subset["avg_exp"] + 1) * 20  # Adjust scaling as needed
        plt.scatter(
            subset["log2FC_treitz"], subset["log2FC_schmidt"],
            color=color, alpha=0.35, edgecolors="none", s=sizes, label=cat
        )
    
    # # Annotate top five up-regulated concordant points.
    # top5_up = concordant_up.sort_values(by="log2FC_treitz", ascending=False).head(5)
    # for _, row in top5_up.iterrows():
    #     plt.annotate(
    #         row["gene"],
    #         xy=(row["log2FC_treitz"], row["log2FC_schmidt"]),
    #         xytext=(row["log2FC_treitz"] + 0.1, row["log2FC_schmidt"] + 0.1),
    #         fontsize=9, color="black"
    #     )
    
    # # Annotate top five down-regulated concordant points.
    # top5_down = concordant_down.sort_values(by="log2FC_treitz", ascending=True).head(5)
    # for _, row in top5_down.iterrows():
    #     plt.annotate(
    #         row["gene"],
    #         xy=(row["log2FC_treitz"], row["log2FC_schmidt"]),
    #         xytext=(row["log2FC_treitz"] - 0.1, row["log2FC_schmidt"] - 0.1),
    #         fontsize=9, color="black"
    #     )
    
    plt.axvline(0, color="gray", linestyle="--")
    plt.axhline(0, color="gray", linestyle="--")
    plt.xlabel("Protein Log2 Fold Change (Treitz et al.)")
    plt.ylabel("Protein Log2 Fold Change (Schmidt et al.)")
    plt.title("Concordance of proteomic datasets\nM9-Acetate versus M9-Glucose")
    plt.legend(frameon=False)
    plt.annotate(f"Pearson r = {r_val:.2f}", xy=(0.05, 0.95), xycoords="axes fraction",
                 fontsize=12, color="black")
    sns.despine()
    
    # Save the plot.
    plot_dir = batch_dir / full_config["datasets"]["treitz"]["output"]["plot_subdir"]
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_file = plot_dir / "treitz_schmidt_acetate_vs_glucose_concordant.png"
    plt.savefig(plot_file, dpi=300)
    plt.close()

    # Sort and save CSV files.
    # Remove any color column before saving CSVs (just in case).
    merged = merged.drop(columns=["color"], errors="ignore")

    # Sort the output CSVs (up and down) by log2FC_treitz and avg_exp.
    concordant_up = concordant_up.sort_values(
        by=["log2FC_treitz", "avg_exp"], ascending=[False, False]
    )
    concordant_down = concordant_down.sort_values(
        by=["log2FC_treitz", "avg_exp"], ascending=[True, False]
    )

    merged.to_csv(csv_dir / "treitz_schmidt_merged.csv", index=False)
    concordant_up.to_csv(csv_dir / "treitz_schmidt_concordant_up.csv", index=False)
    concordant_down.to_csv(csv_dir / "treitz_schmidt_concordant_down.csv", index=False)

    union_count = len(union_genes)
    print(
        f"[Treitz-Schmidt Concordant] Completed. Merged {merged_count} genes "
        f"({merged_count/union_count:.1%} of union)."
    )


if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)
    run_treitz_schmidt_concordant_pipeline(full_config)