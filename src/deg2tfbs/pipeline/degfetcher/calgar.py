"""
--------------------------------------------------------------------------------
<deg2tfbs project>
caglar.py

Module for processing and analyzing data from Calgar et al., which compared gene 
expression for lactate versus glucose media. The dataset is loaded from a CSV file 
in the dnadesign-data repository. This module filters the data to retain only the 
'lactateVSglucose' comparison, partitions the data by dataType ("mrna" or "protein"),
saves separate CSVs, and creates two plots:

1. An MA plot (BaseMean vs. Log2FoldChange) that distinguishes mRNA and protein 
   measurements using different point shapes and colors points red or green based 
   on a ±1.5 log₂FoldChange threshold.
2. A correlation plot comparing log₂FoldChange between protein and mRNA for genes 
   with both measurements. In this plot, only dashed reference lines at the origin 
   (x=0, y=0) are shown, and points are categorized as:
      - Purple: Concordant (both mRNA and protein pass the threshold in the same direction)
      - Blue: RNA specific (only mRNA passes the threshold)
      - Red: Protein specific (only protein passes the threshold)
      - Gray: Unchanged (neither passes the threshold or discordant)
   The plot title also displays the R² value computed on the concordant genes.

"The E. coli molecular phenotype under different growth conditions"
DOI: 10.1038/srep45303

Module Author(s): Eric South
Dunlop Lab
--------------------------------------------------------------------------------
"""
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import yaml
from scipy.stats import pearsonr

from deg2tfbs.pipeline.utils import load_dataset


def read_caglar_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the Calgar dataset using keys from the YAML config.

    Expected config structure:
      {
         "dataset_key": "calgar",
         "sheet_name": "41598_2017_BFsrep45303_MOESM59_",
         "usecols": [
            "gene_name",
            "testVSbase",
            "growthPhase",
            "growthPhase.y",
            "investigatedEffect",
            "dataType",
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "padj"
         ],
         "header": 0
      }
    Performs the following filtering steps:
      1. Retains only rows where 'testVSbase' equals "lactateVSglucose".
      2. Retains only rows where 'growthPhase.y' equals "Exp".
      3. Retains only rows where 'investigatedEffect' equals "batchNumberPLUScarbonSource".
    Also renames 'gene_name' to 'gene'.
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    df.rename(columns={"gene_name": "gene"}, inplace=True)
    # Filter for lactateVSglucose comparison.
    df = df[df["testVSbase"] == "lactateVSglucose"].copy()
    # Additional filters:
    df = df[df["growthPhase.y"] == "Exp"]
    df = df[df["investigatedEffect"] == "batchNumberPLUScarbonSource"]
    return df


def caglar_ma_plot(df: pd.DataFrame, threshold: float, plot_path: Path):
    """
    Generates an MA plot for the Calgar dataset.
    - X-axis: baseMean (log-scaled)
    - Y-axis: log2FoldChange
    Uses different point shapes to distinguish dataType ('mrna' vs 'protein') and 
    colors points red (if log2FoldChange ≥ threshold), green (if ≤ -threshold), 
    or lightgray otherwise.
    """
    plt.figure(figsize=(7, 6))
    marker_dict = {"mrna": "o", "protein": "^"}
    
    for dt, marker in marker_dict.items():
        subset = df[df["dataType"].str.lower() == dt]
        colors = np.where(
            subset["log2FoldChange"] >= threshold, "red",
            np.where(subset["log2FoldChange"] <= -threshold, "green", "lightgray")
        )
        plt.scatter(
            subset["baseMean"], subset["log2FoldChange"],
            marker=marker, c=colors, alpha=0.30, edgecolors="none", label=dt
        )
    
    plt.xscale("log")
    plt.axhline(threshold, color="gray", linestyle="--")
    plt.axhline(-threshold, color="gray", linestyle="--")
    plt.title("MA Plot: BaseMean vs. Log2 Fold Change (Lactate vs Glucose)")
    plt.xlabel("Base Mean")
    plt.ylabel("Log2 Fold Change")
    plt.legend(title="Data Type", frameon=False)
    sns.despine()
    plt.savefig(plot_path, dpi=300)
    plt.close()


def caglar_correlation_plot(df_mrna: pd.DataFrame, df_protein: pd.DataFrame,
                            threshold: float, plot_path: Path):
    """
    Generates a correlation plot comparing log2FoldChange between protein and mRNA.
    Merges the two data subsets on 'gene'. The x-axis represents RNA Log2 (Fold Change) 
    and the y-axis represents Protein log2(Fold Change). Genes are categorized as:
      - 'concordant': both RNA and protein surpass the threshold (in the same direction)
      - 'RNA specific': only RNA surpasses the threshold
      - 'Protein specific': only protein surpasses the threshold
      - 'unchanged': neither measurement surpasses the threshold (or discordant)
    A Pearson correlation is computed for the concordant genes and annotated on the plot.
    Only reference dashed lines at x=0 and y=0 are drawn.
    """
    merged = pd.merge(df_mrna, df_protein, on="gene", suffixes=("_mrna", "_protein"))
    
    def assign_category(row, t):
        r = row["log2FoldChange_mrna"]
        p = row["log2FoldChange_protein"]
        if (abs(r) >= t) and (abs(p) >= t) and (np.sign(r) == np.sign(p)):
            return "concordant"
        elif abs(r) >= t and abs(p) < t:
            return "RNA specific"
        elif abs(p) >= t and abs(r) < t:
            return "Protein specific"
        else:
            return "unchanged"
    
    merged["category"] = merged.apply(lambda row: assign_category(row, threshold), axis=1)
    
    # Compute Pearson correlation for concordant genes
    concordant = merged[merged["category"] == "concordant"]
    if len(concordant) > 1:
        r_val, p_val = pearsonr(concordant["log2FoldChange_mrna"], concordant["log2FoldChange_protein"])
    else:
        r_val = np.nan

    plt.figure(figsize=(7, 6))
    cat_colors = {
        "concordant": "purple",
        "RNA specific": "blue",
        "Protein specific": "red",
        "unchanged": "gray"
    }
    
    for cat, color in cat_colors.items():
        subset = merged[merged["category"] == cat]
        plt.scatter(
            subset["log2FoldChange_mrna"], subset["log2FoldChange_protein"],
            color=color, alpha=0.30, edgecolors="none", label=cat
        )
    
    # Draw reference dashed lines at the origin
    plt.axvline(0, color="gray", linestyle="--")
    plt.axhline(0, color="gray", linestyle="--")
    
    plt.xlabel("RNA Log2 (Fold Change)")
    plt.ylabel("Protein log2(Fold Change)")
    plt.title("Correlation Plot: RNA vs. Protein Log2 Fold Change")
    plt.legend(frameon=False)
    
    # Annotate Pearson correlation coefficient on an optimal position
    if not np.isnan(r_val):
        plt.annotate(f"Pearson r = {r_val:.2f}", xy=(0.05, 0.95), xycoords="axes fraction",
                     fontsize=12, color="black")
    
    sns.despine()
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    return merged


def run_caglar_pipeline(full_config: dict):
    """
    Loads the Calgar dataset, filters rows for the 'lactateVSglucose' comparison and 
    additional experimental conditions ('growthPhase.y' == "Exp" and 'investigatedEffect'
    == "batchNumberPLUScarbonSource"), partitions data by dataType, saves CSV outputs, and 
    creates two plots:
      1. An MA plot (BaseMean vs. Log2FoldChange) with distinct point shapes.
      2. A correlation plot comparing log2FoldChange between RNA and protein for genes
         with both measurements.
    
    Also outputs CSV files for concordant up-regulated and concordant down-regulated genes.
    
    Expected YAML config snippet for caglar:
      caglar:
        data:
          dataset_key: "calgar"
          sheet_name: "41598_2017_BFsrep45303_MOESM59_"
          usecols:
            - "gene_name"
            - "testVSbase"
            - "growthPhase"
            - "growthPhase.y"
            - "investigatedEffect"
            - "dataType"
            - "baseMean"
            - "log2FoldChange"
            - "pvalue"
            - "padj"
          header: 0
        thresholds:
          log2_fc_threshold: 1.5
        output:
          csv_subdir: "csvs"
          plot_subdir: "plots"
    """
    config_caglar = full_config["datasets"]["caglar"]
    df = read_caglar_data(config_caglar["data"])
    
    # Build output paths (assumed project structure similar to other pipelines)
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["pipeline"]["stages"]["degfetcher"]["root_dir"]
    batch_id = full_config["pipeline"]["stages"]["degfetcher"]["batch_id"]
    batch_dir = output_root / batch_id
    
    csv_dir = batch_dir / config_caglar["output"]["csv_subdir"]
    plot_dir = batch_dir / config_caglar["output"]["plot_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    
    # Partition data by dataType
    df_mrna = df[df["dataType"].str.lower() == "mrna"].copy()
    df_protein = df[df["dataType"].str.lower() == "protein"].copy()
    
    # Save partitioned CSVs
    mrna_csv_path = csv_dir / "caglar_lactateVSglucose_mrna.csv"
    protein_csv_path = csv_dir / "caglar_lactateVSglucose_protein.csv"
    df_mrna.to_csv(mrna_csv_path, index=False)
    df_protein.to_csv(protein_csv_path, index=False)
    
    # Retrieve threshold for log2 fold change
    threshold = config_caglar["thresholds"]["log2_fc_threshold"]
    assert isinstance(threshold, (int, float)), "Threshold must be numeric."
    
    # Create MA plot using the complete filtered dataset
    ma_plot_path = plot_dir / "caglar_lactateVSglucose_MA_plot.png"
    caglar_ma_plot(df, threshold, ma_plot_path)
    
    # Create correlation plot for genes measured by both methods and capture merged data
    correlation_plot_path = plot_dir / "caglar_lactateVSglucose_correlation_plot.png"
    merged = caglar_correlation_plot(df_mrna, df_protein, threshold, correlation_plot_path)
    
    # Output concordant up- and concordant down-regulated genes as separate CSVs
    concordant_up = merged[
        (merged["category"] == "concordant") &
        (merged["log2FoldChange_mrna"] > 0) &
        (merged["log2FoldChange_protein"] > 0)
    ]
    concordant_down = merged[
        (merged["category"] == "concordant") &
        (merged["log2FoldChange_mrna"] < 0) &
        (merged["log2FoldChange_protein"] < 0)
    ]
    
    concordant_up.to_csv(csv_dir / "caglar_lactateVSglucose_concordant_up.csv", index=False)
    concordant_down.to_csv(csv_dir / "caglar_lactateVSglucose_concordant_down.csv", index=False)
    
    print(f"[Calgar Pipeline] Completed. Processed {len(df)} entries for 'lactateVSglucose'.")


if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)
    run_caglar_pipeline(full_config)