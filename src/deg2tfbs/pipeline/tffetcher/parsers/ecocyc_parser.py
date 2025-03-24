"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tffetcher/parsers/ecocyc_parser.py

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

from pathlib import Path
from typing import Dict, Set, Tuple
from collections import defaultdict
import re

from deg2tfbs.pipeline.tffetcher import normalize_gene

class EcoCycParser:
    """
    Class-based parser for EcoCyc's (Release: 28.5) regulatory network file:
    
        ECOLI-regulatory-network.txt

    Returns:
      A tuple: 
         - A dict mapping each regulated gene to a set of (regulator, polarity) tuples.
         - A set of all genes present in the network (including both regulators and regulatees).
    """
    _polarity_regex = re.compile(r'^((?:\+/-)|±|[+-])?')

    def parse(self, network_path: Path) -> Tuple[Dict[str, Set[Tuple[str, str]]], Set[str]]:
        if not network_path.is_file():
            raise FileNotFoundError(f"EcoCyc network file not found: {network_path}")
        network = defaultdict(set)
        all_genes = set()
        current_regulator = None
        line_number = 0

        with open(network_path, 'r') as f:
            for line in f:
                line_number += 1
                raw_line = line.rstrip('\n')
                stripped_line = raw_line.strip()

                if not stripped_line or stripped_line.startswith('#'):
                    continue

                if raw_line.startswith((' ', '\t')):
                    if not current_regulator:
                        raise ValueError(f"Target line without active regulator (line {line_number})")
                    self._process_target_line(stripped_line, current_regulator, network, all_genes, line_number)
                else:
                    if not stripped_line.endswith('*'):
                        raise ValueError(f"Invalid regulator format (line {line_number}): {raw_line}")
                    current_regulator = normalize_gene(stripped_line.rstrip('*').strip())
                    all_genes.add(current_regulator)
        return dict(network), all_genes

    def _process_target_line(self, line: str, current_reg: str,
                             network: defaultdict, all_genes: set, line_number: int) -> None:
        tokens = line.split()
        for token in tokens:
            match = self._polarity_regex.match(token)
            if not match:
                raise ValueError(f"Invalid token format '{token}' (line {line_number})")
            polarity_group = match.group(1)
            if polarity_group:
                gene_part = token[len(polarity_group):].lstrip()
            else:
                gene_part = token

            if polarity_group == '+':
                polarity = '+'
            elif polarity_group == '-':
                polarity = '-'
            elif polarity_group in ('+/-', '±'):
                polarity = '±'
            else:
                polarity = 'NA'
            gene = normalize_gene(gene_part.rstrip('*'))
            if not gene:
                raise ValueError(f"Empty gene name in '{token}' (line {line_number})")
            network[gene].add((current_reg, polarity))
            all_genes.add(gene)
