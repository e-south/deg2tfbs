"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tfbsfetcher/parsers/ecocyc_tfbs_parser.py

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import re
from pathlib import Path
from typing import Dict, Set

import pandas as pd

from deg2tfbs.pipeline.utils import load_tfbs_file


class EcoCycTFBSParser:
    """
    Example parser for:
      'SmartTable_All_Transcription_Factor_Binding_Sites.txt'
    Where the columns (header=3) might be:
       "Site", "Regulator", "Left", "Right", "Strand", ...
       "Sequence - DNA sequence"
    Parse them, apply a regex to the "Regulator" column 
    to remove trailing dash. Also uppercase the site sequence, ignoring any
    lines with empty or null regulator or site sequence.

    Returns a dict[tf -> set(site_sequence)]
    """

    _dash_cleanup_regex = re.compile(r'(.*?)\s*(?:-|dna-binding-site).*$')

    def parse(self) -> Dict[str, Set[str]]:
        tfbs_path = load_tfbs_file("ecocyc_28_tfbs_smart_table")
        
        df = pd.read_csv(tfbs_path, sep='\t', header=2, encoding='iso-8859-1')  # Skip first 3 lines
        required_cols = ["Regulator", "Sequence - DNA sequence"]
        for rc in required_cols:
            if rc not in df.columns:
                raise ValueError(f"[EcoCycTFBSParser] Missing '{rc}' in {tfbs_path.name}")

        # Drop empty
        df["Regulator"] = df["Regulator"].astype(str).str.strip()
        df["Sequence - DNA sequence"] = df["Sequence - DNA sequence"].astype(str).str.strip()
        df = df[(df["Regulator"] != "") & (df["Sequence - DNA sequence"] != "")]

        # Fix the "Regulator" to keep only the portion before dash
        # e.g. "AcrR-proflavin" => "AcrR"
        def clean_reg(r: str) -> str:
            m = self._dash_cleanup_regex.match(r)
            if m:
                return m.group(1).lower().strip()
            # else fallback
            return r.lower().strip()

        df["Regulator"] = df["Regulator"].apply(clean_reg)

        def normalize_seq(seq: str) -> str:
            """Extracts only uppercase DNA sequences from mixed-case input"""
            return "".join([c for c in seq if c.isupper() and c.isalpha()]).strip()

        df["Sequence - DNA sequence"] = df["Sequence - DNA sequence"].apply(normalize_seq)

        # Build dictionary
        tf2sites = {}
        for idx, row in df.iterrows():
            tf = row["Regulator"]
            site_seq = row["Sequence - DNA sequence"]
            if tf not in tf2sites:
                tf2sites[tf] = set()
            tf2sites[tf].add(site_seq)

        return tf2sites
