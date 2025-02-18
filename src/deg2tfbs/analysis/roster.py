"""
--------------------------------------------------------------------------------
<deg2tfbs project>
src/deg2tfbs/analysis/roster.py

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import List, Tuple, Dict

logger = logging.getLogger(__name__)

def extract_unique_regulators(csv_file: Path) -> List[str]:
    """
    Read a CSV file and return a sorted list of unique regulators from the 'tf' column.
    """
    df = pd.read_csv(csv_file)
    if "tf" not in df.columns:
        raise ValueError(f"'tf' column missing in {csv_file}")
    regs = sorted(df["tf"].dropna().unique().tolist())
    logger.info(f"Extracted {len(regs)} unique regulators from {csv_file}")
    return regs

def create_boolean_vector(regs_reference: List[str], csv_file: Path) -> np.array:
    """
    Return a binary vector (as a NumPy array) indicating the presence (1) or absence (0)
    of each regulator (in the order of regs_reference) in the CSV file.
    """
    df = pd.read_csv(csv_file)
    if "tf" not in df.columns:
        raise ValueError(f"'tf' column missing in {csv_file}")
    present = set(df["tf"].dropna().unique().tolist())
    vector = np.array([1 if reg in present else 0 for reg in regs_reference], dtype=np.int8)
    return vector

def build_reference_roster(mapping: Dict[str, Path], reference_key: str) -> Tuple[List[str], np.array]:
    """
    Given a mapping from tfbsbatch names to CSV file paths and a reference key,
    return a tuple (regs_reference, ref_vector) for the reference set.
    """
    if reference_key not in mapping:
        raise ValueError(f"Reference key {reference_key} not found in mapping")
    regs_reference = extract_unique_regulators(mapping[reference_key])
    ref_vector = create_boolean_vector(regs_reference, mapping[reference_key])
    return regs_reference, ref_vector

def build_pairing_roster(regs_reference: List[str], csv_file: Path) -> np.array:
    """
    Build a binary vector for a given tfbsbatch CSV using the provided regs_reference.
    Asserts that all regulators in the CSV are contained in the reference.
    """
    df = pd.read_csv(csv_file)
    if "tf" not in df.columns:
        raise ValueError(f"'tf' column missing in {csv_file}")
    regs_set = set(df["tf"].dropna().unique().tolist())
    missing = regs_set - set(regs_reference)
    if missing:
        raise ValueError(f"Regulators {missing} from {csv_file} are not in the reference set.")
    return create_boolean_vector(regs_reference, csv_file)

def build_matrix_from_rosters(roster_dict: Dict[str, np.array]) -> Tuple[np.array, List[str]]:
    """
    Given a dictionary mapping sample names to binary vectors, return a matrix (rows=samples)
    and a sorted list of sample names.
    """
    sample_names = sorted(roster_dict.keys())
    matrix = np.vstack([roster_dict[sample] for sample in sample_names])
    return matrix, sample_names

def compute_pairing_details(vector1: np.array, vector2: np.array) -> Dict:
    """
    Given two binary vectors (for a pairing), compute:
      - intersection: elementwise AND
      - group1_only: vector1 minus intersection
      - group2_only: vector2 minus intersection
    Also return counts.
    """
    intersection = vector1 & vector2
    group1_only = vector1 - intersection
    group2_only = vector2 - intersection
    details = {
        "vector1": vector1,
        "vector2": vector2,
        "intersection": intersection,
        "group1_only": group1_only,
        "group2_only": group2_only,
        "count_group1": int(vector1.sum()),
        "count_group2": int(vector2.sum()),
        "count_intersection": int(intersection.sum())
    }
    return details

def exclude_intersections(rosters: Dict[str, np.array]) -> Dict[str, np.array]:
    """
    For paired sources (determined by a common prefix before the last underscore),
    remove the intersection. For each group with more than one source, if a regulator is present
    in every source in the group, set that value to 0 in every source.
    """
    groups = {}
    for source in rosters.keys():
        parts = source.rsplit("_", 1)
        group_key = parts[0] if len(parts) == 2 else source
        groups.setdefault(group_key, []).append(source)
    
    updated_rosters = {}
    for group, sources in groups.items():
        if len(sources) > 1:
            intersection = rosters[sources[0]].copy()
            for src in sources[1:]:
                intersection = intersection & rosters[src]
            for src in sources:
                updated_rosters[src] = rosters[src] * (1 - intersection)
        else:
            updated_rosters[sources[0]] = rosters[sources[0]]
    return updated_rosters