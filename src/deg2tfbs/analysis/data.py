"""
--------------------------------------------------------------------------------
<deg2tfbs project>
src/deg2tfbs/analysis/data.py

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

from pathlib import Path
import logging
from typing import Dict

logger = logging.getLogger(__name__)

def get_tfbs_mapping_files(tfbs_dir: Path, file_pattern: str = "*.csv") -> Dict[str, Path]:
    """
    Scan the provided directory (tfbs_dir) for subdirectories whose names start with "tfbsbatch_"
    and return a mapping from each subdirectory name to the first CSV file found.
    
    Preconditions:
      - tfbs_dir must be a valid directory.
    
    Returns:
      - A dictionary mapping subdirectory names to CSV file paths.
    """
    assert tfbs_dir.is_dir(), f"tfbs_dir must be a directory: {tfbs_dir}"
    mapping: Dict[str, Path] = {}
    for subdir in tfbs_dir.iterdir():
        if subdir.is_dir() and subdir.name.startswith("tfbsbatch_"):
            # Look for a "csvs" subdirectory first.
            csv_dir = subdir / "csvs"
            if csv_dir.is_dir():
                csv_files = list(csv_dir.glob(file_pattern))
            else:
                csv_files = list(subdir.glob(file_pattern))
            if csv_files:
                mapping[subdir.name] = csv_files[0]
                logger.info(f"Found CSV in {subdir.name}: {csv_files[0]}")
            else:
                logger.warning(f"No CSV file found in {subdir}")
    return mapping
