"""
--------------------------------------------------------------------------------
<deg2tfbs project>
deg2tfbs/pipeline/tffetcher/parsers/__init__.py

Any parser should implement the following transformation:

    parse(...) -> Dict[str, Set[Tuple[str, str]]]

This way tffetcher.py can treat each parser identically.

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

from typing import Union
from .ecocyc_parser import EcoCycParser
from .regdb_parser import RegulonDBParser

def get_regulatory_parser(parser_name: str) -> Union[EcoCycParser, RegulonDBParser]:
    """
    Simple factory to return an instance of the requested parser.
    If you want to pass extra arguments (like confidence filter),
    you could parse them from the config here, or set them directly.
    """
    parser_name = parser_name.lower()
    if parser_name.startswith("ecocyc"):
        return EcoCycParser()
    elif parser_name.startswith("regulondb"):
        # For example, we can pass confidence_filter=True by default
        return RegulonDBParser(confidence_filter=True)
    else:
        raise ValueError(f"Unknown parser: {parser_name}")
