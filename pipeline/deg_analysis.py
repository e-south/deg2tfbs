"""
--------------------------------------------------------------------------------
<deg2tfbs project>

Module: deg_analysis.py
Identify DEGs (differentially expressed genes)

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# deg_analysis.py
from deg2tfbs.analysis.datasets import ceroni, emani, mori, ...

def run_deg_analysis(dataset_name, df, config):
    if dataset_name == 'ceroni':
        return ceroni.run(df, config)
    elif dataset_name == 'emani':
        return emani.run(df, config)
    ...


#####

def identify_degs(df: pd.DataFrame,
                  log2_fc_col: str,
                  pval_col: str,
                  log2_fc_threshold: float,
                  pval_threshold: float,
                  output_csv: str = None,
                  plot_file: str = None) -> pd.DataFrame:
    """
    Identify differentially expressed genes based on user-defined thresholds,
    e.g., log2 fold change >= 1.0 and p-value <= 0.05.
    
    :param df: DataFrame with gene expression data (containing log2 fold-change, p-value).
    :param log2_fc_col: Column name for log2 fold-change data.
    :param pval_col: Column name for p-value data.
    :param log2_fc_threshold: Threshold for absolute log2 fold-change.
    :param pval_threshold: Threshold for p-value (or adjusted p-value).
    :param output_csv: Path to save CSV of DEGs. If None, skipping save.
    :param plot_file: Path to save a diagnostic volcano/MA plot. If None, skip plotting.
    :return: Filtered DataFrame of identified DEGs.
    """
    # Clean up the input DataFrame
    df = df.copy().dropna(subset=[log2_fc_col, pval_col])
    
    # Mark significant DEGs
    df["sig"] = (df[log2_fc_col].abs() >= log2_fc_threshold) & (df[pval_col] <= pval_threshold)
    degs = df[df["sig"]].copy()

    # Optionally plot
    if plot_file:
        plt.figure(figsize=(6,5))
        sns.scatterplot(
            x=df[log2_fc_col],
            y=-np.log10(df[pval_col]),
            hue=df["sig"],
            alpha=0.3,
            edgecolor=None
        )
        # Draw thresholds
        plt.axvline(x= log2_fc_threshold, color='red', ls='--')
        plt.axvline(x=-log2_fc_threshold, color='red', ls='--')
        # approximate pval threshold line
        plt.axhline(y=-np.log10(pval_threshold), color='blue', ls='--')
        
        plt.xlabel(f'{log2_fc_col} (log2-fold change)')
        plt.ylabel(f'-log10({pval_col})')
        plt.title('Volcano plot of DEGs')
        plt.legend(title='Significant DEG')
        plt.tight_layout()
        plt.savefig(plot_file)
        plt.close()

    # Optionally save CSV
    if output_csv:
        degs.to_csv(output_csv, index=False)
    
    return degs
