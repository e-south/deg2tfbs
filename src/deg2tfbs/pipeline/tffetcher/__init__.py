"""
--------------------------------------------------------------------------------
<deg2tfbs project>

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

def normalize_gene(gene: str) -> str:
    """
    Normalize a gene name by lowercasing, stripping whitespace,
    and removing non-alphanumeric characters.
    """
    import re
    return re.sub(r'[^a-z0-9]', '', gene.lower().strip())
