"""
--------------------------------------------------------------------------------
<deg2tfbs project>
src/deg2tfbs/analysis/compare.py

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import pandas as pd
import logging
from pathlib import Path
from typing import List, Dict

logger = logging.getLogger(__name__)

def generate_pairing_dataframe(regs_reference: List[str], details: Dict) -> pd.DataFrame:
    """
    Create a DataFrame with one row per regulator. Columns include:
      - regulator
      - group1: binary (0/1) for the first dataset
      - group2: binary (0/1) for the second dataset
      - intersection: 1 if regulator is present in both
      - group1_only: 1 if only in group1
      - group2_only: 1 if only in group2
    """
    df = pd.DataFrame({
        "regulator": regs_reference,
        "group1": details["vector1"],
        "group2": details["vector2"],
        "intersection": details["intersection"],
        "group1_only": details["group1_only"],
        "group2_only": details["group2_only"]
    })
    return df

def print_pairing_summary(pair_name: str, details: Dict, total: int) -> str:
    """
    Print and return a human-readable summary for a pairing.
    """
    count1 = details["count_group1"]
    count2 = details["count_group2"]
    count_int = details["count_intersection"]
    pct1 = (count1 / total) * 100
    pct2 = (count2 / total) * 100
    pct_int = (count_int / total) * 100
    summary = (f"Pairing {pair_name}:\n"
               f"  Group1: {count1} regulators ({pct1:.1f}% of reference)\n"
               f"  Group2: {count2} regulators ({pct2:.1f}% of reference)\n"
               f"  Intersection: {count_int} regulators ({pct_int:.1f}% of reference)\n")
    print(summary)
    logger.info(summary)
    return summary

def save_pairing_dataframe(df: pd.DataFrame, output_file: Path) -> None:
    """
    Save the detailed pairing DataFrame to CSV.
    """
    df.to_csv(output_file, index=False)
    logger.info(f"Saved pairing details to {output_file}")
