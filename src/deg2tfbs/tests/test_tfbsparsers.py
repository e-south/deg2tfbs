"""
--------------------------------------------------------------------------------
<deg2tfbs project>
tests/test_parsers.py

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import pandas as pd
import pytest
from deg2tfbs.pipeline.tfbsfetcher.parsers.regdb_tfbs_parser import RegulonDBTFBSParser

def test_regulondb_parser_nucleotides(monkeypatch):
    """
    Test that the RegulonDB parser returns site sequences containing only A, T, C, G.
    """
    # Create a dummy DataFrame to simulate file input.
    # The parser expects column names with numbered prefixes that it will clean.
    dummy_df = pd.DataFrame({
        "4)regulatorName": ["AcrR", "MarA"],
        "4)tfrsSeq": ["acgTACgg", "ATGcTAG"]
    })

    # Override pd.read_csv to return the dummy DataFrame regardless of the file path.
    def dummy_read_csv(*args, **kwargs):
        return dummy_df
    monkeypatch.setattr(pd, "read_csv", dummy_read_csv)
    
    # Monkeypatch the load_tfbs_file function so that it returns a dummy path.
    monkeypatch.setattr(
        "deg2tfbs.pipeline.tfbsfetcher.parsers.regdb_tfbs_parser.load_tfbs_file",
        lambda key: "dummy_path.tsv"
    )
    
    parser = RegulonDBTFBSParser()
    result = parser.parse()
    
    valid_nucleotides = set("ATCG")
    for tf, sites in result.items():
        for site in sites:
            assert site != "", "Site sequence should not be empty."
            # Ensure the site contains only uppercase A, T, C, G.
            assert set(site).issubset(valid_nucleotides), (
                f"Site '{site}' for TF '{tf}' contains invalid characters."
            )


#### Test for the EcoCyc Parser

from deg2tfbs.pipeline.tfbsfetcher.parsers.ecocyc_tfbs_parser import EcoCycTFBSParser

def test_ecocyc_parser_nucleotides(monkeypatch):
    """
    Test that the EcoCyc parser returns site sequences containing only A, T, C, G.
    """
    # Create a dummy DataFrame that simulates the file input.
    dummy_df = pd.DataFrame({
        "Regulator": ["AcrR-proflavin", "MarA-dna-binding-site"],
        "Sequence - DNA sequence": ["acgTACgg", "ATGcTAG"]
    })
    
    def dummy_read_csv(*args, **kwargs):
        return dummy_df
    monkeypatch.setattr(pd, "read_csv", dummy_read_csv)
    
    monkeypatch.setattr(
        "deg2tfbs.pipeline.tfbsfetcher.parsers.ecocyc_tfbs_parser.load_tfbs_file",
        lambda key: "dummy_path.txt"
    )
    
    parser = EcoCycTFBSParser()
    result = parser.parse()
    
    valid_nucleotides = set("ATCG")
    for tf, sites in result.items():
        for site in sites:
            assert site != "", "Site sequence should not be empty."
            # Ensure the site contains only uppercase A, T, C, G.
            assert set(site).issubset(valid_nucleotides), (
                f"Site '{site}' for TF '{tf}' contains invalid characters."
            )
