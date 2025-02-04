"""
--------------------------------------------------------------------------------
<deg2tfbs project>
tests/test_tfbsfetcher.py

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

# tests/test_tfbsfetcher.py

import pandas as pd
import pytest

from deg2tfbs.pipeline.tfbsfetcher.tfbsfetcher import (
    aggregate_deg_source,
    merge_tfbs_dicts
)

# -------------------------------------------------------------------
# Tests for aggregate_deg_source
# -------------------------------------------------------------------

def test_aggregate_deg_source():
    """
    Test that aggregate_deg_source splits hyphenated tokens,
    removes duplicates while preserving first occurrence order.
    """
    series = pd.Series([
        "sourceA_down",
        " sourceB_down ",
        "sourceA_down",
        "sourceC_down",
        "  "  # extra whitespace that should be ignored
    ])
    expected = "sourceA_down-sourceB_down-sourceC_down"
    result = aggregate_deg_source(series)
    assert result == expected, f"Expected '{expected}' but got '{result}'"


def test_aggregate_deg_source_empty():
    """
    Test that aggregate_deg_source returns an empty string when the series
    contains only empty strings or whitespace.
    """
    series = pd.Series(["", "   ", None, "   "])
    result = aggregate_deg_source(series)
    assert result == "", f"Expected empty string, got '{result}'"


def test_aggregate_deg_source_order():
    """
    Test that aggregate_deg_source preserves the order of tokens based on
    first occurrence even when duplicate tokens appear in different orders.
    
    For example:
      - The first value "s2_down-s1_down" yields tokens: ["s2_down", "s1_down"]
      - The second value "s1_down-s3_down" yields tokens: ["s1_down", "s3_down"]
      The final aggregated string should be "s2_down-s1_down-s3_down".
    """
    series = pd.Series(["s2_down-s1_down", "s1_down-s3_down"])
    expected = "s2_down-s1_down-s3_down"
    result = aggregate_deg_source(series)
    assert result == expected, f"Expected '{expected}', got '{result}'"


# -------------------------------------------------------------------
# Tests for deduplication logic
# -------------------------------------------------------------------

def deduplicate_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Simulate the deduplication block from tfbsfetcher.py.
    
    Uniqueness is defined by the tuple: (tf, tfbs). When duplicates are found,
    the following aggregation rules are applied:
    
      - gene: alphabetical order and underscore-joined.
      - deg_source: aggregated via aggregate_deg_source.
      - polarity: mode (most common value) of the group.
      - tfbs_source: alphabetical order and underscore-joined.
      - is_sigma_factor and is_global_regulator: first value.
    """
    group_cols = ["tf", "tfbs"]
    duplicate_mask = df.duplicated(subset=group_cols, keep=False)
    num_duplicates = duplicate_mask.sum()
    if num_duplicates:
        aggregation_rules = {
            "gene": lambda x: "_".join(sorted(set(x))),
            "deg_source": aggregate_deg_source,
            "polarity": lambda x: x.mode()[0] if not x.empty else "NA",
            "tfbs_source": lambda x: "_".join(sorted(set(x))),
            "is_sigma_factor": "first",
            "is_global_regulator": "first"
        }
        agg_rules = {col: rule for col, rule in aggregation_rules.items() if col not in group_cols and col in df.columns}
        df_clean = df.groupby(group_cols, as_index=False).agg(agg_rules)
    else:
        df_clean = df.drop_duplicates(subset=group_cols)
    return df_clean


def test_deduplication_logic():
    """
    Simulate a DataFrame with duplicate (tf, tfbs) pairs and test that
    after deduplication the output DataFrame has the expected aggregated values.
    
    For instance, suppose we have three rows for TF 'tf1' with binding site 'ATCG':
      - Row 1: gene='geneA', deg_source='s1_down', polarity='+', tfbs_source='source1'
      - Row 2: gene='geneB', deg_source=' s2_down' (with extra whitespace), polarity='+', tfbs_source='source2'
      - Row 3: gene='geneA', deg_source='s1_down', polarity='-', tfbs_source='source1'
    
    The aggregated values should be:
      - gene: "geneA_geneB" (sorted alphabetically)
      - deg_source: aggregate_deg_source(["s1_down", "s2_down", "s1_down"]) → "s1_down-s2_down"
      - polarity: mode of ["+", "+", "-"] → "+"
      - tfbs_source: sorted(set(["source1", "source2"])) → "source1_source2"
      - is_sigma_factor and is_global_regulator: "first" value, e.g. "no"
    """
    data = {
        "tf": ["tf1", "tf1", "tf1", "tf2"],
        "tfbs": ["ATCG", "ATCG", "ATCG", "CGTA"],
        "gene": ["geneA", "geneB", "geneA", "geneC"],
        "deg_source": ["s1_down", " s2_down", "s1_down", "s3_down"],
        "polarity": ["+", "+", "-", "+"],
        "tfbs_source": ["source1", "source2", "source1", "source3"],
        "is_sigma_factor": ["no", "no", "no", "yes"],
        "is_global_regulator": ["no", "no", "no", "yes"],
    }
    df = pd.DataFrame(data)
    
    # Apply deduplication to simulate the aggregation step.
    df_clean = deduplicate_df(df)
    
    # We expect two rows: one for ('tf1', 'ATCG') and one for ('tf2', 'CGTA').
    # For the tf1 row:
    # - gene: "geneA_geneB" (alphabetically sorted)
    # - deg_source: aggregate_deg_source(["s1_down", "s2_down", "s1_down"]) -> "s1_down-s2_down"
    # - polarity: mode of ["+", "+", "-"] -> "+"
    # - tfbs_source: "source1_source2"
    # - is_sigma_factor: "no" (first value)
    # - is_global_regulator: "no" (first value)
    
    # Extract the deduplicated row for tf1.
    tf1_row = df_clean[(df_clean["tf"] == "tf1") & (df_clean["tfbs"] == "ATCG")].iloc[0]
    
    assert tf1_row["gene"] == "geneA_geneB", f"Expected gene aggregation to be 'geneA_geneB', got '{tf1_row['gene']}'"
    assert tf1_row["deg_source"] == "s1_down-s2_down", f"Expected deg_source to be 's1_down-s2_down', got '{tf1_row['deg_source']}'"
    assert tf1_row["polarity"] == "+", f"Expected polarity to be '+', got '{tf1_row['polarity']}'"
    assert tf1_row["tfbs_source"] == "source1_source2", f"Expected tfbs_source to be 'source1_source2', got '{tf1_row['tfbs_source']}'"
    assert tf1_row["is_sigma_factor"] == "no", f"Expected is_sigma_factor to be 'no', got '{tf1_row['is_sigma_factor']}'"
    assert tf1_row["is_global_regulator"] == "no", f"Expected is_global_regulator to be 'no', got '{tf1_row['is_global_regulator']}'"
    
    # Also verify that the unique (tf, tfbs) pairs in the cleaned DataFrame are as expected.
    expected_keys = {("tf1", "ATCG"), ("tf2", "CGTA")}
    result_keys = set(df_clean[["tf", "tfbs"]].itertuples(index=False, name=None))
    assert result_keys == expected_keys, f"Expected keys {expected_keys}, got {result_keys}"


# -------------------------------------------------------------------
# Additional tests for merge_tfbs_dicts
# -------------------------------------------------------------------

def test_merge_tfbs_dicts():
    """
    Test that merge_tfbs_dicts performs a union across multiple dictionaries.
    """
    dict1 = {"tf1": {"ATCG", "CGTA"}, "tf2": {"GGGG"}}
    dict2 = {"tf1": {"TTTT"}, "tf3": {"CCCC"}}
    merged = merge_tfbs_dicts([dict1, dict2])
    
    assert merged["tf1"] == {"ATCG", "CGTA", "TTTT"}
    assert merged["tf2"] == {"GGGG"}
    assert merged["tf3"] == {"CCCC"}


def test_merge_tfbs_dicts_empty():
    """
    Test that merge_tfbs_dicts returns an empty dictionary when given an empty list.
    """
    result = merge_tfbs_dicts([])
    assert result == {}, f"Expected empty dict, got {result}"


def test_merge_tfbs_dicts_with_empty_dict():
    """
    Test that merge_tfbs_dicts correctly handles a list containing empty dictionaries.
    """
    dict1 = {}
    dict2 = {"tf1": {"ATCG"}}
    result = merge_tfbs_dicts([dict1, dict2])
    expected = {"tf1": {"ATCG"}}
    assert result == expected, f"Expected {expected}, got {result}"


def test_merge_tfbs_dicts_overlapping():
    """
    Test merge_tfbs_dicts with overlapping keys that have overlapping sets.
    """
    dict1 = {"tf1": {"ATCG"}}
    dict2 = {"tf1": {"CGTA"}, "tf2": {"GGGG"}}
    dict3 = {"tf1": {"ATCG", "TTTT"}, "tf2": {"CCCC"}}
    result = merge_tfbs_dicts([dict1, dict2, dict3])
    expected = {"tf1": {"ATCG", "CGTA", "TTTT"}, "tf2": {"GGGG", "CCCC"}}
    assert result == expected, f"Expected {expected}, got {result}"

