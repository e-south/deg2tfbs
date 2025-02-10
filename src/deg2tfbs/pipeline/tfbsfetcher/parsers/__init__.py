"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tfbsfetcher/parsers/__init__.py

Any parser should implement the following transformation:

    parse(...) -> Dict[str, Set[str]]

This way tfbsfetcher.py can treat each parser identically.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

from .ecocyc_tfbs_parser import EcoCycTFBSParser
from .regdb_tfbs_parser import RegulonDBTFBSParser

def get_tfbs_parser(parser_name: str):
    """
    Return an instance of either EcoCycTFBSParser or RegulonDBTFBSParser
    depending on the key.
    """    
    parser_name = parser_name.lower()
    if parser_name.startswith("ecocyc"):
        return EcoCycTFBSParser()
    elif parser_name.startswith("regdb"):
        return RegulonDBTFBSParser()
    else:
        raise ValueError(f"Unknown parser: {parser_name}")
