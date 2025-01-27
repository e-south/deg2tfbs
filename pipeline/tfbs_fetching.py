"""
--------------------------------------------------------------------------------
<deg2tfbs project>

Module: tfbs_identification.py
Identify TFBS for the TFs associated with DEGs

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import pandas as pd

def find_tfbs_for_tfs(tfs_df: pd.DataFrame,
                      tf_col: str,
                      tfbs_db_df: pd.DataFrame,
                      tfbs_regulator_col: str = "regulator",
                      output_csv: str = None) -> pd.DataFrame:
    """
    Given a list of TFs (from deg->tf association), retrieve matching TFBS from a TFBS database.
    For instance, a RegulonDB or EcoCyc-based DF with columns: [regulator, sequence, start, end, ...]
    
    :param tfs_df: DataFrame containing TF info, e.g. the output from `map_degs_to_tfs`.
    :param tf_col: Column name in `tfs_df` that lists the TF name.
    :param tfbs_db_df: DataFrame of TFBS data from an external source, with a 'regulator' col.
    :param tfbs_regulator_col: The column in `tfbs_db_df` that denotes the TF's name.
    :param output_csv: If given, path to save the results.
    :return: DataFrame of all TFBS that match the TFs in `tfs_df`.
    """
    unique_tfs = tfs_df[tf_col].dropna().unique()
    tfbs_filtered = tfbs_db_df[tfbs_db_df[tfbs_regulator_col].isin(unique_tfs)].copy()
    
    if output_csv:
        tfbs_filtered.to_csv(output_csv, index=False)
    
    return tfbs_filtered
