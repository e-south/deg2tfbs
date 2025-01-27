"""
--------------------------------------------------------------------------------
<deg2tfbs project>

Module: tf_association.py
Map DEGs to Transcription Factors (TFs)

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import pandas as pd

def map_degs_to_tfs(degs: pd.DataFrame,
                    gene_col: str,
                    tf_gene_interaction_df: pd.DataFrame,
                    gene_col_in_interaction: str = "Regulated Gene",
                    tf_col_in_interaction: str = "TF name",
                    output_csv: str = None) -> pd.DataFrame:
    """
    Map a list of DEGs to associated TFs by joining on a TF-Gene interaction table
    (e.g., from RegulonDB).
    
    :param degs: DataFrame with DEGs, including a column for gene IDs/names.
    :param gene_col: Column name in `degs` representing the gene (e.g. "gene_name").
    :param tf_gene_interaction_df: DataFrame of known TF->gene interactions.
    :param gene_col_in_interaction: Column name in TF-gene DataFrame for the gene.
    :param tf_col_in_interaction: Column name for the TF names in the interaction DataFrame.
    :param output_csv: If provided, path to save the resulting DataFrame.
    :return: DataFrame with DEGs joined to their associated TFs.
    """
    # Join the DEGs to TF relationships (left join to keep only DEGs)
    merged = degs.merge(
        tf_gene_interaction_df,
        left_on=gene_col,
        right_on=gene_col_in_interaction,
        how="left"
    )

    # Optionally drop rows with no TF match if you prefer
    # merged = merged.dropna(subset=[tf_col_in_interaction])

    # Save if needed
    if output_csv:
        merged.to_csv(output_csv, index=False)

    return merged
