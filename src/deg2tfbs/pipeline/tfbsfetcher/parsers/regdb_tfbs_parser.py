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

from deg2tfbs.pipeline.utils import load_tfbs_file
import pandas as pd

class RegulonDBTFBSParser:
    """
    Example parser for:
      'TF-RISet.tsv' (RegulonDB 13)
    Skip ~44 lines of header. Then columns:
      (4) regulatorName => e.g. 'AcrR'
      (10) tfrsSeq => might have uppercase plus flanking lowercases
         Keep uppercase portion only.

    Returns: dict[tf -> set(site_sequence)]
    """

    def parse(self) -> Dict[str, Set[str]]:
        tsv_path = load_tfbs_file("regulondb_13_tf_ri_set")

        # Skip ~44 lines of comments
        skip_rows = 44
        df = pd.read_csv(tsv_path, sep='\t', header=0, skiprows=skip_rows, dtype=str)

        # "4)regulatorName" => "regulatorName" after cleaning
        df.columns = (
            df.columns
            .str.replace(r'^\d+\)', '', regex=True)  # Remove numbered prefixes like "4)"
            .str.strip()                             
            .str.lower()                              
        )

        required_cols = ["regulatorname", "tfrsseq"]
        for rc in required_cols:
            if rc not in df.columns:
                raise ValueError(f"[RegulonDBTFBSParser] Missing column '{rc}' in {tsv_path.name}")

        df["regulatorname"] = df["regulatorname"].astype(str).str.lower().str.strip()
        df["tfrsseq"] = df["tfrsseq"].astype(str).str.strip()

        # Drop empties
        df = df[(df["regulatorname"] != "") & (df["tfrsseq"] != "")]

        def normalize_seq(seq: str) -> str:
            """Strictly preserves original uppercase segments only"""
            return "".join([c for c in seq if c.isupper() and c.isalpha()]).strip()

        df["tfrsseq"] = df["tfrsseq"].apply(normalize_seq)

        # Build dictionary
        tf2sites = {}
        for idx, row in df.iterrows():
            tf = row["regulatorname"]
            site_seq = row["tfrsseq"]
            if tf not in tf2sites:
                tf2sites[tf] = set()
            tf2sites[tf].add(site_seq)

        return tf2sites
