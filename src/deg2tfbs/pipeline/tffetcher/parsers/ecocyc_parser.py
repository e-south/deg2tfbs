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

class EcoCycParser:
    """
    Class-based parser for EcoCyc's (Release: 28.5) regulatory network file:
    
        ECOLI-regulatory-network.txt
    
    Returns: dict[gene -> set of (regulator, polarity)]
    """


    # Regex captures all valid polarities or none
    _polarity_regex = re.compile(r'^((?:\+/-)|±|[+-])?')

    def parse(self, network_path: Path) -> Dict[str, Set[Tuple[str, str]]]:
        if not network_path.is_file():
            raise FileNotFoundError(f"EcoCyc network file not found: {network_path}")
        reverse_network = defaultdict(set)
        current_regulator = None
        line_number = 0

        with open(network_path, 'r') as f:
            for line in f:
                line_number += 1
                raw_line = line.rstrip('\n')
                stripped_line = raw_line.strip()

                # Skip comments and empty lines
                if not stripped_line or stripped_line.startswith('#'):
                    continue

                if raw_line.startswith((' ', '\t')):
                    # This is a regulatee line for current_regulator
                    if not current_regulator:
                        raise ValueError(f"Target line without active regulator (line {line_number})")
                    self._process_target_line(stripped_line, current_regulator, reverse_network, line_number)
                else:
                    # This is a potential regulator line (should end with '*')
                    if not stripped_line.endswith('*'):
                        raise ValueError(f"Invalid regulator format (line {line_number}): {raw_line}")
                    # Strip trailing '*'
                    reg = stripped_line.rstrip('*').strip().lower()
                    current_regulator = reg

        return dict(reverse_network)

    
    def _process_target_line(self, line: str, current_reg: str,
                            network: defaultdict, line_number: int) -> None:
        """
        Process a target line from EcoCyc, capturing polarity or using 'NA' if not present.
        """
        tokens = line.split()
        for token in tokens:
            match = self._polarity_regex.match(token)
            if not match:
                raise ValueError(f"Invalid token format '{token}' (line {line_number})")
            
            # Extract polarity group and gene name
            polarity_group = match.group(1)
            if polarity_group:
                gene_part = token[len(polarity_group):].lstrip()
            else:
                gene_part = token  # No polarity prefix
            
            # Determine polarity
            if polarity_group == '+':
                polarity = '+'
            elif polarity_group == '-':
                polarity = '-'
            elif polarity_group in ('+/-', '±'):
                polarity = '±'
            else:
                polarity = 'NA'  # No recognized polarity
            
            gene = gene_part.rstrip('*').strip().lower()
            if not gene:
                raise ValueError(f"Empty gene name in '{token}' (line {line_number})")
            
            network[gene].add((current_reg, polarity))

