"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tfbsfetcher/parsers/regulondb_tfbs_parser.py

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import re
from pathlib import Path
from typing import Dict, Set

import pandas as pd
from deg2tfbs.pipeline.utils import load_tfbs_file

class RegulonDBTFBSParser:
    """
    Example parser for:
      'TF-RISet.tsv' (RegulonDB 13)

    This parser reads a TSV file that starts with ~44 lines of comments. It then expects columns such as:
      - (4) regulatorName  -> e.g. 'AcrR'
      - (10) tfrsSeq       -> The binding site sequence (contains uppercase letters for the core site)
      - (20) confidenceLevel -> Confidence indication ('C' for confirmed, 'S' for strong, or 'W' for weak)

    Only rows where the confidence level is either "C" (confirmed) or "S" (strong) will be used.
    
    Returns:
      dict mapping transcription factor names (str) to a set of binding site sequences (set[str])
    """

    def parse(self) -> Dict[str, Set[str]]:
        tsv_path = load_tfbs_file("regulondb_13_tf_ri_set")

        # Skip ~44 lines of header comments
        skip_rows = 44
        df = pd.read_csv(tsv_path, sep='\t', header=0, skiprows=skip_rows, dtype=str)

        # Clean column names: remove numbered prefixes and normalize to lowercase.
        df.columns = (
            df.columns
            .str.replace(r'^\d+\)', '', regex=True)
            .str.strip()
            .str.lower()
        )

        # Ensure that the required columns exist.
        required_cols = ["regulatorname", "tfrsseq", "confidencelevel"]
        for rc in required_cols:
            if rc not in df.columns:
                raise ValueError(f"[RegulonDBTFBSParser] Missing column '{rc}' in {tsv_path.name}")

        # Clean the critical columns.
        df["regulatorname"] = df["regulatorname"].astype(str).str.lower().str.strip()
        df["tfrsseq"] = df["tfrsseq"].astype(str).str.strip()
        df["confidencelevel"] = df["confidencelevel"].astype(str).str.upper().str.strip()

        # Drop rows with empty transcription factor names or binding site sequences.
        df = df[(df["regulatorname"] != "") & (df["tfrsseq"] != "")]

        # Filter the DataFrame to only include rows with the allowed confidence levels.
        # Allowed values: "C" for confirmed and "S" for strong.
        allowed_confidences = {"C", "S"}
        df = df[df["confidencelevel"].isin(allowed_confidences)]

        def normalize_seq(seq: str) -> str:
            """Strictly preserves original uppercase segments only from a given sequence."""
            return "".join([c for c in seq if c.isupper() and c.isalpha()]).strip()

        # Apply the uppercase normalization to extract the binding site from tfrsseq.
        df["tfrsseq"] = df["tfrsseq"].apply(normalize_seq)

        # Build dictionary mapping TF names to a set of binding site sequences.
        tf2sites = {}
        for idx, row in df.iterrows():
            tf = row["regulatorname"]
            site_seq = row["tfrsseq"]
            if tf not in tf2sites:
                tf2sites[tf] = set()
            tf2sites[tf].add(site_seq)

        return tf2sites