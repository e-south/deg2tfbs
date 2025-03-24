"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tffetcher/parsers/redb_parser.py

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import pandas as pd
from pathlib import Path
from typing import Dict, Set, Tuple
from collections import defaultdict

from deg2tfbs.pipeline.tffetcher import normalize_gene

class RegulonDBParser:
    """
    Class-based parser for RegulonDB's (Release: 13) regulatory network file:
    
        NetworkRegulatorGene.tsv

    Returns: dict[gene -> set of (regulator, polarity)]
    """
    def __init__(self, confidence_filter: bool = True):
        self.confidence_filter = confidence_filter

    def parse(self, tsv_path: Path) -> Dict[str, Set[Tuple[str, str]]]:
        if not tsv_path.is_file():
            raise FileNotFoundError(f"RegulonDB file not found: {tsv_path}")
        df = pd.read_csv(
            tsv_path,
            sep='\t',
            skiprows=24,
            header=0,
            dtype=str,
            encoding='utf-8-sig'
        )
        df.columns = (
            df.columns.str.strip()
            .str.replace(r'[)\s]', '', regex=True)
            .str.lower()
        )
        required_cols = [
            '2regulatorname',
            '5regulatedname',
            '6function',
            '7confidencelevel'
        ]
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            raise ValueError(f"Missing columns: {missing}")
        df = df[required_cols]
        df = df.dropna(subset=required_cols[:3])
        
        if self.confidence_filter:
            df['7confidencelevel'] = df['7confidencelevel'].astype(str).str.strip()
            df = df[df['7confidencelevel'].isin({'C', 'S'})]

        network = defaultdict(set)
        for _, row in df.iterrows():
            regulator = normalize_gene(row['2regulatorname'])
            gene = normalize_gene(row['5regulatedname'])
            func = str(row['6function']).strip()
            polarity = {'+': '+', '-': '-', '-+': '±'}.get(func, 'NA')
            network[gene].add((regulator, polarity))
        return dict(network)
