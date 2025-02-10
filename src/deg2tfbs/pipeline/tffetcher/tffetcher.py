"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tffetcher.py

Maps of DEGs -> TFs using one or more regulatory network resources (e.g., files
from EcoCyc, RegulonDB, etc.), according to the YAML config.

Read a CSV(s) from degfetcher containing columns:
  [gene, source, comparison]
  
Then loads either EcoCyc or RegulonDB 'regulatory network' datasets (or both),
parses them, locates instances of target genes and grabs corresponding regulators
from the same row, and produces a final CSV with columns (example):
    gene, regulator, polarity, source, is_global_regulator, is_sigma_factor, deg_source

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
import glob
import logging
import re
from pathlib import Path
from typing import Dict, Set, Tuple, List
from collections import defaultdict

import pandas as pd

# Import parser factory
from deg2tfbs.pipeline.tffetcher.parsers import get_regulatory_parser
from deg2tfbs.pipeline import utils

logging.basicConfig(level=logging.INFO)

# Flag master regulators; categorized here https://ecocyc.org/overviewsWeb/regOv.shtml?orgid=ECOLI#
MASTER_REGULATORS = {
    'flhd', 'flhc', 'fis', 'fhla', 'dksa', 'cysb', 'csgd', 'crp', 'cra', 'cpxr', 'barr', 'argr',
    'arca', 'srsr', 'spf', 'soxs', 'slya', 'ryhb', 'rob', 'rcsb', 'purr', 'pohp', 'phob', 'pdhr',
    'oxyr', 'ompr', 'nsrr', 'nhar', 'narp', 'narl', 'nagc', 'nac', 'mode', 'mara', 'lrp', 'lexa',
    'iscr', 'ihfb', 'ihfa', 'hns', 'glng', 'glar', 'gcvb', 'gadx', 'gade', 'fur', 'fnr',
}
# Flag sigma factors; categorized here https://ecocyc.org/ECOLI/NEW-IMAGE?type=ENZYME&object=Sigma-Factors
SIGMA_FACTORS = {
    'flia', 'feci', 'rpod', 'rpos', 'rpon', 'rpoh', 'rpoe'
}


def run_tffetcher_stage(config: dict) -> None:
    """
    Main entry point for the TF-fetcher stage.

    1) Reads CSV files from degfetcher (optionally filtered).
    2) Builds gene -> deg_source set (e.g. 'bie_down', 'kim_up', etc.).
    3) Loads and merges each requested regulatory network.
    4) For each (gene, regulator), merges polarisations and network sources,
       sets is_global_regulator / is_sigma_factor, and outputs a final CSV.

    Raises:
        ValueError, FileNotFoundError if any config keys or paths are invalid.
    """
    if "root_dir" not in config or "batch_id" not in config:
        raise ValueError("tffetcher config must have 'root_dir' and 'batch_id' keys")

    # Build the final output paths
    project_root = Path(__file__).parent.parent.parent
    
    # Output folder, e.g. deg2tfbs/pipeline/tffetcher
    root_dir = (project_root / Path(config['root_dir'])).resolve()
    batch_id = config["batch_id"]

    input_conf = config.get("input", {})
    deg_batch_id = input_conf.get("deg_batch_id")
    deg_csv_subdir = input_conf.get("deg_csv_subdir", "csvs")
    
    # Look for a group override
    group_def = input_conf.get("deg_csv_group")
    if group_def:
        # If a filter is provided, create a glob pattern that matches CSV files with that substring.
        if "filter" in group_def:
            filter_str = group_def["filter"]
            deg_includes = [f"*{filter_str}*.csv"]
        # If an explicit file list is provided, extract the file names.
        elif "files" in group_def:
            # Each file entry is expected to be a dict with at least a "file" key.
            deg_includes = [entry["file"] for entry in group_def.get("files", []) if "file" in entry]
        else:
            # If the group exists but has no recognized key, default to an empty list (or you could default to all CSVs).
            deg_includes = []
    else:
        # No group override provided: use the legacy key if it exists.
        deg_includes = input_conf.get("deg_csv_includes", [])
    
    if not deg_batch_id:
        raise ValueError("tffetcher config missing 'deg_batch_id' in 'input'")

    # Example: /path/to/deg2tfbs/pipeline/degfetcher/degbatch_20250130/csvs
    deg_csv_dir = (project_root / "pipeline" / "degfetcher" / deg_batch_id / deg_csv_subdir).resolve()

    if not deg_csv_dir.is_dir():
        raise FileNotFoundError(f"DEG CSV input dir not found: {deg_csv_dir}")

    sources_conf = config.get("sources", {}).get("regulatory_networks", {})
    if not sources_conf:
        raise ValueError("No 'regulatory_networks' configured under 'sources'")

    params_conf = config.get("params", {})
    network_strategy = params_conf.get("network_strategy", "union")

    # ------------------------------------------------------
    # STEP 1: Collect DEGs from CSVs, optionally filtering
    # ------------------------------------------------------
    logging.info(f"[tffetcher] Gathering genes from the following files: {', '.join(deg_includes)}")
    gene_to_deg_sources = collect_degs_with_sources(deg_csv_dir, deg_includes=deg_includes)
    num_files_considered = get_num_files_considered(deg_csv_dir, deg_includes)
    logging.info(f"[tffetcher] -> Found {len(gene_to_deg_sources)} unique genes from {num_files_considered} CSV file(s)\n")

    if not gene_to_deg_sources:
        logging.warning("[tffetcher] No DEGs found (possibly none matched deg_csv_includes?). Exiting.")
        return

    # Genes are simply the keys in the dictionary
    all_genes = set(gene_to_deg_sources.keys())

    # ------------------------------------------------------
    # STEP 2: Load & parse each enabled regulatory network
    # ------------------------------------------------------
    networks = []
    total_network_pairs = 0  # Track how many gene->regulator pairs found across networks

    for net_name, net_conf in sources_conf.items():
        enabled = net_conf.get("enabled", False)
        if not enabled:
            logging.info(f"[tffetcher] Skipping disabled network: {net_name}")
            continue

        parser_name = net_conf.get("parser")
        data_key = net_conf.get("path")  # e.g. "ecocyc_28_reg_network"

        if not parser_name:
            raise ValueError(f"Network '{net_name}' missing 'parser' key")
        if not data_key:
            raise ValueError(f"Network '{net_name}' missing 'path' key")

        # Use utils.DATA_FILES to find the actual path:
        try:
            network_file_path = utils.DATA_FILES[data_key]
        except KeyError:
            available_keys = ", ".join(utils.DATA_FILES.keys())
            raise ValueError(f"Invalid data_key '{data_key}'. Must be one of: {available_keys}")

        if not network_file_path.exists():
            raise FileNotFoundError(f"File not found for key '{data_key}': {network_file_path}")

        parser_obj = get_regulatory_parser(parser_name)

        logging.info(f"[tffetcher] Parsing '{parser_name}':")
        
        # This returns Dict[gene -> set of (regulator, polarity)]
        net_dict = parser_obj.parse(network_file_path)

        # Filter out any gene not in `all_genes`
        filtered_net_dict = {}
        for gene, regset in net_dict.items():
            if gene in all_genes:
                filtered_net_dict[gene] = regset

        # Log how many relevant genes remain
        relevant_pairs = sum(len(regs) for regs in filtered_net_dict.values())
        logging.info(f"[tffetcher] {net_name} => {len(filtered_net_dict)} input DEGs found, ~{relevant_pairs} gene-regulator edges after filtering\n")

        # Wrap it in a structure that also tracks the "source" (ecocyc or regdb)
        net_aug = defaultdict(lambda: defaultdict(lambda: {"polarity_set": set(), "sources": set()}))

        short_label = derive_short_network_label(net_name)

        for gene, regset in filtered_net_dict.items():
            for (regulator, polarity) in regset:
                net_aug[gene][regulator]["polarity_set"].add(polarity)
                net_aug[gene][regulator]["sources"].add(short_label)

        networks.append(net_aug)

    if not networks:
        logging.warning("[tffetcher] No regulatory networks enabled. Exiting.")
        return

    # ------------------------------------------------------
    # STEP 3: Merge results from different regulatory networks if more than one
    # ------------------------------------------------------
    merged_network = networks[0]
    for other_net in networks[1:]:
        if network_strategy == "union":
            merged_network = union_networks(merged_network, other_net)
        elif network_strategy == "intersect":
            merged_network = intersect_networks(merged_network, other_net)
        elif network_strategy == "ecocyc_only":
            # Skip merging with others
            pass
        elif network_strategy == "regulondb_only":
            # skip merging with others
            pass
        else:
            # default fallback is union
            merged_network = union_networks(merged_network, other_net)

    logging.info(f"[tffetcher] Merging strategy used: {network_strategy}")

    # A quick count of merged pairs:
    final_merged_pairs = 0
    for gene_dict in merged_network.values():
        final_merged_pairs += len(gene_dict)

    logging.info(f"[tffetcher] => After merging, there are ~{final_merged_pairs} total (gene-regulator) relationships")

    # ------------------------------------------------------
    # STEP 4: For each gene, gather its regulators from merged_network
    # Then produce final rows with required columns
    # ------------------------------------------------------
    rows = []
    for gene in sorted(all_genes):
        # Standardize gene
        gene_norm = gene.lower().strip()
        if gene_norm not in merged_network:
            continue

        # "deg_source" col listing all CSVs for that gene, e.g. 'bie_down-kim_down-sanchez_up'
        deg_source_str = "-".join(sorted(gene_to_deg_sources[gene]))

        # For each regulator of that gene
        for regulator, info in merged_network[gene_norm].items():
            polarity_str = unify_polarities(info["polarity_set"])
            net_sources = info["sources"]  # e.g. {'ecocyc_28', 'regdb_13'}

            # If more than one network, join them with _AND_ as requested
            if len(net_sources) > 1:
                source_str = "_AND_".join(sorted(net_sources))
            else:
                source_str = next(iter(net_sources))  # single element

            is_global = "yes" if is_master_regulator(regulator) else "no"
            is_sigma = "yes" if is_sigma_factor(regulator) else "no"

            rows.append({
                "gene": gene_norm,
                "regulator": regulator,
                "polarity": polarity_str,
                "source": source_str,
                "is_global_regulator": is_global,
                "is_sigma_factor": is_sigma,
                "deg_source": deg_source_str,
            })

    # ------------------------------------------------------
    # STEP 5: Write out results
    # ------------------------------------------------------
    out_dir = root_dir / batch_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "deg2tf_mapping.csv"

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv, index=False)

    logging.info(f"[tffetcher] => Writing final CSV with {len(df_out)} rows")


def collect_degs_with_sources(csv_dir: Path, deg_includes: List[str] = None) -> Dict[str, Set[str]]:
    """
    Collect a dictionary: gene -> set of deg_sources.
    For example, if a row in "bie_downregulated_degs.csv" has gene=ompR,
    add "bie_down" to the set for gene=ompR.

    If `deg_includes` is provided (e.g. ["*downregulated_degs.csv", "sanchez_vasquez_upregulated_degs.csv"]),
    only consider those CSV filenames.

    Returns:
        gene_to_sources: dict where each key is a gene,
                         and each value is e.g. {"bie_down", "kim_down", "sanchez_vasquez_up"}
    """
    gene_to_sources: Dict[str, Set[str]] = defaultdict(set)

    matching_files = []
    if deg_includes:
        # If user gave a list of patterns, gather all matching CSVs
        for pattern in deg_includes:
            matched = glob.glob(str(csv_dir / pattern))
            matching_files.extend(matched)
        if not matching_files:
            logging.warning(f"[tffetcher] deg_csv_includes was {deg_includes}, but no files matched!")
    else:
        # No filters => take all CSVs
        matching_files = glob.glob(str(csv_dir / "*.csv"))

    for csv_fp in sorted(set(matching_files)):
        filename = os.path.basename(csv_fp)
        # e.g. "bie_downregulated_degs.csv" => "bie_down"
        deg_label = derive_deg_label(filename)

        df = pd.read_csv(csv_fp)
        if "gene" not in df.columns:
            logging.warning(f"[tffetcher] Skipping {filename}, no 'gene' column.")
            continue

        for gene_val in df["gene"].dropna():
            g = gene_val.lower().strip()
            gene_to_sources[g].add(deg_label)

    return dict(gene_to_sources)


def get_num_files_considered(csv_dir: Path, deg_includes: List[str]) -> int:
    """
    Helper function to simply return how many CSVs ended up being included,
    to display in the logs. This is optional and purely for debug info.
    """
    if deg_includes:
        all_files = []
        for pattern in deg_includes:
            all_files.extend(glob.glob(str(csv_dir / pattern)))
        return len(set(all_files))
    else:
        return len(glob.glob(str(csv_dir / "*.csv")))


def union_networks(
    net_a: Dict[str, Dict[str, Dict[str, Set[str]]]],
    net_b: Dict[str, Dict[str, Dict[str, Set[str]]]]
) -> Dict[str, Dict[str, Dict[str, Set[str]]]]:
    """
    Merges two nested dicts of the form:
    gene -> regulator -> { "polarity_set": set(...), "sources": set(...) }

    For each gene/regulator pair, union the polarity_set, union the sources.
    """
    merged = {}
    for gene_a, regs_dict_a in net_a.items():
        merged[gene_a] = {}
        for reg_a, info_a in regs_dict_a.items():
            merged[gene_a][reg_a] = {
                "polarity_set": set(info_a["polarity_set"]),
                "sources": set(info_a["sources"])
            }

    for gene_b, regs_dict_b in net_b.items():
        if gene_b not in merged:
            merged[gene_b] = {}
        for reg_b, info_b in regs_dict_b.items():
            if reg_b not in merged[gene_b]:
                merged[gene_b][reg_b] = {
                    "polarity_set": set(),
                    "sources": set()
                }
            merged[gene_b][reg_b]["polarity_set"].update(info_b["polarity_set"])
            merged[gene_b][reg_b]["sources"].update(info_b["sources"])
    return merged


def intersect_networks(
    net_a: Dict[str, Dict[str, Dict[str, Set[str]]]],
    net_b: Dict[str, Dict[str, Dict[str, Set[str]]]]
) -> Dict[str, Dict[str, Dict[str, Set[str]]]]:
    """
    Return intersection of two networks in the same nested dict format.
    Keep only genes that appear in both; for each gene, keep only regulators that appear in both.
    Union their polarity and sources for those shared regulators.
    """
    merged = {}
    common_genes = set(net_a.keys()).intersection(set(net_b.keys()))
    for gene in common_genes:
        merged[gene] = {}
        regs_a = net_a[gene]
        regs_b = net_b[gene]
        common_regs = set(regs_a.keys()).intersection(set(regs_b.keys()))
        for reg in common_regs:
            pol_set = set(regs_a[reg]["polarity_set"]) | set(regs_b[reg]["polarity_set"])
            src_set = set(regs_a[reg]["sources"]) | set(regs_b[reg]["sources"])
            merged[gene][reg] = {
                "polarity_set": pol_set,
                "sources": src_set
            }
    return merged


def unify_polarities(polarity_set: Set[str]) -> str:
    """
    If there's exactly one polarity in the set, return it.
    If there are multiple distinct polarities (e.g. '+' and '-'),
    unify to '±'. If the set is empty or 'NA' is the only value, return 'NA'.
    """
    pol_clean = {p for p in polarity_set if p not in ["NA", None, ""]}
    if not pol_clean:
        return "NA"
    if len(pol_clean) == 1:
        return pol_clean.pop()
    # If there's conflicting polarities, use '±'
    return "±"


def derive_short_network_label(net_name: str) -> str:
    """
    Simple way to map net_name from config (e.g. 'ecocyc' or 'regulondb')
    to a label like 'ecocyc_28' or 'regdb_13'.
    """
    net_name = net_name.lower().strip()
    if "ecocyc" in net_name:
        return "ecocyc_28"
    elif "regulon" in net_name:
        return "regdb_13"
    return net_name


def derive_deg_label(filename: str) -> str:
    """
    Extract a short label from a CSV filename, e.g.:
      'bie_downregulated_degs.csv' -> 'bie_down'
      'kim_upregulated_degs.csv'   -> 'kim_up'
      'ceroni_downregulated_degs.csv' -> 'ceroni_down'
    If neither 'upregulated' nor 'downregulated' is in the filename,
    just strip off '.csv'.
    """
    base = filename.lower().replace('.csv', '')
    if 'upregulated' in base:
        return base.replace('_upregulated_degs', '_up')
    elif 'downregulated' in base:
        return base.replace('_downregulated_degs', '_down')
    else:
        return base


def is_master_regulator(regulator: str) -> bool:
    """
    Check if a regulator is in the MASTER_REGULATORS set, ignoring case.
    """
    reg_norm = regulator.lower().strip()
    return reg_norm in MASTER_REGULATORS


def is_sigma_factor(regulator: str) -> bool:
    """
    Check if a regulator is in the SIGMA_FACTORS set, ignoring case.
    """
    reg_norm = regulator.lower().strip()
    return reg_norm in SIGMA_FACTORS
