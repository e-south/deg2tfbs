"""
--------------------------------------------------------------------------------
<deg2tfbs project>
pipeline/tffetcher/tffetcher.py

This module maps differentially expressed genes (DEGs) to transcription factors (TFs)
using one or more regulatory network resources (e.g., EcoCyc, RegulonDB). It produces:
  - A consolidated deg2tf_mapping CSV showing the TF-DEG associations (filtered by the DEG CSV files)
  - An enrichment analysis of TFs based on the full (unfiltered) regulatory network.
    
Additional metadata is recorded for up- and down-regulated DEGs (as indicated by the 'deg_source' field),
but the enrichment test is performed on the union of all DEG files.
  
Output directories are consolidated under a single tfbatch_<date> directory.
  
Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
import glob
import logging
import re
from pathlib import Path
from typing import Dict, Set
from collections import defaultdict

import pandas as pd

from deg2tfbs.pipeline.tffetcher.parsers import get_regulatory_parser
from deg2tfbs.pipeline.tffetcher.regulator_utils import (
    is_master_regulator,
    is_sigma_factor,
    is_nucleoid_regulator
)
from deg2tfbs.pipeline.tffetcher import tfenrichment, normalize_gene
from deg2tfbs.pipeline import utils
from deg2tfbs.pipeline import filtering

logging.basicConfig(level=logging.INFO)

MASTER_REGULATORS = {
    'flhd', 'flhc', 'fis', 'fhla', 'dksa', 'cysb', 'csgd', 'crp', 'cra', 'cpxr',
    'barr', 'argr', 'arca', 'srsr', 'spf', 'soxs', 'slya', 'ryhb', 'rob', 'rcsb',
    'purr', 'pohp', 'phob', 'pdhr', 'oxyr', 'ompr', 'nsrr', 'nhar', 'narp', 'narl',
    'nagc', 'nac', 'mode', 'mara', 'lrp', 'lexa', 'iscr', 'ihf', 'ihfb', 'ihfa', 'hns',
    'glng', 'glar', 'gcvb', 'gadx', 'gade', 'fur', 'fnr',
}
NUCLEOID_REGULATORS = {"ihf", "ihfa", "ihfb", "fis", "h-ns", "dps", "dan"}
SIGMA_FACTORS = {'flia', 'feci', 'rpod', 'rpos', 'rpon', 'rpoh', 'rpoe'}

def run_tffetcher_stage(config: dict) -> None:
    """
    Main entry point for the TF-fetcher stage.
    
    Steps:
      1. Collect DEGs from CSV files, building a mapping (gene -> DEG source).
      2. Parse each regulatory network twice:
         a. A filtered network (only for genes present in the DEG CSV files) to build mapping rows.
         b. The full network (unfiltered or filtered using a merged gene set) to serve as the background for enrichment analysis.
      3. Merge the filtered networks (for mapping) and the full networks (for enrichment) separately.
      4. Construct mapping rows from the filtered (mapping) network and write a consolidated deg2tf_mapping CSV.
      5. Run enrichment analysis using the full (unfiltered or filtered) merged network as background.
      
    Note:
      The enrichment analysis calculates for each TF:
        - a: The number of DEGs (from the union of all DEG files) that are targets of the TF.
        - K: The total number of targets for the TF (from the full network).
        - Enrichment score: a/K.
        - Fisher's exact test p-value from a 2x2 table:
              [ a, K - a ]
              [ M - a, N - K - (M - a) ]
          where M = number of DEGs (intersection with full network) and N = total genes in the full network.
        - FDR (BH-corrected p-value) and a significance flag based on the user-specified threshold.
      
    The output is written to a consolidated tfbatch directory using the batch_id provided in the config.
    
    Raises:
      ValueError, FileNotFoundError for configuration or file errors.
    """
    if "root_dir" not in config or "batch_id" not in config:
        raise ValueError("tffetcher config must have 'root_dir' and 'batch_id' keys")
    
    project_root = Path(__file__).parent.parent.parent
    root_dir = (project_root / Path(config['root_dir'])).resolve()
    batch_id = config["batch_id"]
    
    input_conf = config.get("input", {})
    deg_batch_id = input_conf.get("deg_batch_id")
    deg_csv_subdir = input_conf.get("deg_csv_subdir", "csvs")
    
    if "deg_csv_groups" in input_conf:
        groups = input_conf["deg_csv_groups"]
        deg_files = []
        for group_name, group_def in groups.items():
            if "files" in group_def:
                for entry in group_def["files"]:
                    if "file" in entry:
                        deg_files.append(entry["file"])
            elif "filter" in group_def:
                deg_files.extend(glob.glob(str(Path(deg_csv_subdir) / f"*{group_def['filter']}*.csv")))
        deg_includes = list(set(deg_files))
    else:
        deg_includes = input_conf.get("deg_csv_includes", [])
    
    if not deg_batch_id:
        raise ValueError("tffetcher config missing 'deg_batch_id' in 'input'")
    
    deg_csv_dir = (project_root / "pipeline" / "degfetcher" / deg_batch_id / deg_csv_subdir).resolve()
    if not deg_csv_dir.is_dir():
        raise FileNotFoundError(f"DEG CSV input dir not found: {deg_csv_dir}")
    
    filter_gene_networks = input_conf.get("filter_gene_networks", False)
    if filter_gene_networks:
        merged_csv_name = input_conf.get("filter_gene_network_csv")
        if not merged_csv_name:
            raise ValueError("Filtering enabled but no 'filter_gene_network_csv' provided in config.")
        merged_csv_path = (project_root / "pipeline" / "degfetcher" / deg_batch_id / deg_csv_subdir / merged_csv_name).resolve()
        if not merged_csv_path.exists():
            raise FileNotFoundError(f"Merged CSV file for filtering not found: {merged_csv_path}")
        df_merged = pd.read_csv(merged_csv_path)
        if "gene" not in df_merged.columns:
            raise ValueError("Merged CSV file does not have a 'gene' column")
        merged_gene_set = set(normalize_gene(g) for g in df_merged["gene"].dropna())
        logging.info(f"[tffetcher] Filtering full regulatory network using merged gene set from {merged_csv_path}")
    else:
        merged_gene_set = None

    all_unfiltered_genes = set()
    
    sources_conf = config.get("sources", {}).get("regulatory_networks", {})
    if not sources_conf:
        raise ValueError("No 'regulatory_networks' configured under 'sources'")
    
    params_conf = config.get("params", {})
    network_strategy = params_conf.get("network_strategy", "union")
    
    logging.info(f"[tffetcher] Gathering genes from the following files: {', '.join(deg_includes)}")
    gene_to_deg_sources = collect_degs_with_sources(deg_csv_dir, deg_includes=deg_includes)
    num_files_considered = get_num_files_considered(deg_csv_dir, deg_includes)
    logging.info(f"[tffetcher] -> Found {len(gene_to_deg_sources)} unique DEGs from {num_files_considered} CSV file(s)\n")
    if not gene_to_deg_sources:
        logging.warning("[tffetcher] No DEGs found. Exiting.")
        return
    filtered_degs = set(gene_to_deg_sources.keys())
    
    mapping_networks = []
    enrichment_networks = []
    for net_name, net_conf in sources_conf.items():
        if not net_conf.get("enabled", False):
            logging.info(f"[tffetcher] Skipping disabled network: {net_name}")
            continue
        
        parser_name = net_conf.get("parser")
        data_key = net_conf.get("path")
        if not parser_name:
            raise ValueError(f"Network '{net_name}' missing 'parser' key")
        if not data_key:
            raise ValueError(f"Network '{net_name}' missing 'path' key")
        
        try:
            network_file_path = utils.DATA_FILES[data_key]
        except KeyError:
            available_keys = ", ".join(utils.DATA_FILES.keys())
            raise ValueError(f"Invalid data_key '{data_key}'. Must be one of: {available_keys}")
        
        if not network_file_path.exists():
            raise FileNotFoundError(f"File not found for key '{data_key}': {network_file_path}")
        
        parser_obj = get_regulatory_parser(parser_name)
        logging.info(f"[tffetcher] Parsing '{parser_name}':")
        if parser_name.lower().startswith("ecocyc"):
            net_dict, network_genes = parser_obj.parse(network_file_path)
            all_unfiltered_genes.update(network_genes)
        else:
            net_dict = parser_obj.parse(network_file_path)
            all_unfiltered_genes.update(net_dict.keys())
            for regset in net_dict.values():
                for regulator, _ in regset:
                    all_unfiltered_genes.add(normalize_gene(regulator))
        
        filtered_net_dict = {gene: regset for gene, regset in net_dict.items() if gene in filtered_degs}
        if merged_gene_set is not None:
            full_net_dict = {gene: regset for gene, regset in net_dict.items() if gene in merged_gene_set}
        else:
            full_net_dict = net_dict
        
        short_label = derive_short_network_label(net_name)
        
        net_aug_mapping = defaultdict(lambda: defaultdict(lambda: {"polarity_set": set(), "sources": set()}))
        for gene, regset in filtered_net_dict.items():
            for (regulator, polarity) in regset:
                net_aug_mapping[gene][regulator]["polarity_set"].add(polarity)
                net_aug_mapping[gene][regulator]["sources"].add(short_label)
        mapping_networks.append(net_aug_mapping)
        
        net_aug_enrichment = defaultdict(lambda: defaultdict(lambda: {"polarity_set": set(), "sources": set()}))
        for gene, regset in full_net_dict.items():
            for (regulator, polarity) in regset:
                net_aug_enrichment[gene][regulator]["polarity_set"].add(polarity)
                net_aug_enrichment[gene][regulator]["sources"].add(short_label)
        enrichment_networks.append(net_aug_enrichment)
    
    if merged_gene_set is not None:
        num_merged = len(merged_gene_set)
        num_present = len(merged_gene_set.intersection(all_unfiltered_genes))
        dropout_prop = 1.0 - (num_present / num_merged) if num_merged > 0 else 0
        logging.info(f"[tffetcher] Unfiltered networks contain {len(all_unfiltered_genes)} unique genes.")
        logging.info(f"[tffetcher] Merged gene set contains {num_merged} unique genes.")
        logging.info(f"[tffetcher] Proportion of merged genes found in networks: {(num_present/num_merged):.2%} (Dropout: {dropout_prop:.2%})")
        dropout_genes = sorted(merged_gene_set.difference(all_unfiltered_genes))
        out_dir = root_dir / batch_id
        out_dir.mkdir(parents=True, exist_ok=True)
        csvs_dir = out_dir / "csvs"
        csvs_dir.mkdir(parents=True, exist_ok=True)
        dropout_file = csvs_dir / "dropout_genes.txt"
        with open(dropout_file, "w") as f:
            for gene in dropout_genes:
                f.write(f"{gene}\n")
        logging.info(f"[tffetcher] Saved dropout genes ({len(dropout_genes)} genes) to {dropout_file}")
    
    if not mapping_networks:
        logging.warning("[tffetcher] No regulatory networks enabled. Exiting.")
        return
    
    merged_mapping_network = mapping_networks[0]
    for net in mapping_networks[1:]:
        if network_strategy == "union":
            merged_mapping_network = union_networks(merged_mapping_network, net)
        elif network_strategy == "intersect":
            merged_mapping_network = intersect_networks(merged_mapping_network, net)
        else:
            merged_mapping_network = union_networks(merged_mapping_network, net)
    
    merged_enrichment_network = enrichment_networks[0]
    for net in enrichment_networks[1:]:
        if network_strategy == "union":
            merged_enrichment_network = union_networks(merged_enrichment_network, net)
        elif network_strategy == "intersect":
            merged_enrichment_network = intersect_networks(merged_enrichment_network, net)
        else:
            merged_enrichment_network = union_networks(merged_enrichment_network, net)
    
    logging.info(f"[tffetcher] Merging strategy used: {network_strategy}")
    final_merged_pairs = sum(len(g_dict) for g_dict in merged_mapping_network.values())
    logging.info(f"[tffetcher] => After merging (mapping network), there are ~{final_merged_pairs} total (gene-regulator) relationships")
    
    rows = []
    for gene in sorted(filtered_degs):
        gene_norm = gene.lower().strip()
        if gene_norm not in merged_mapping_network:
            continue
        deg_source_str = "-".join(sorted(gene_to_deg_sources[gene]))
        for regulator, info in merged_mapping_network[gene_norm].items():
            polarity_str = unify_polarities(info["polarity_set"])
            net_sources = info["sources"]
            source_str = "_AND_".join(sorted(net_sources)) if len(net_sources) > 1 else next(iter(net_sources))
            is_global = "yes" if is_master_regulator(regulator.lower()) else "no"
            is_sigma = "yes" if is_sigma_factor(regulator.lower()) else "no"
            is_nucleo = "yes" if is_nucleoid_regulator(regulator.lower(), NUCLEOID_REGULATORS) else "no"
            rows.append({
                "gene": gene_norm,
                "regulator": regulator,
                "polarity": polarity_str,
                "source": source_str,
                "is_global_regulator": is_global,
                "is_sigma_factor": is_sigma,
                "is_nucleoid_regulator": is_nucleo,
                "deg_source": deg_source_str,
            })
    
    # Remove any transcription factors that do not contain binding sites.
    rows = filtering.filter_mapping_rows_by_tfbs(rows)
    
    out_dir = root_dir / batch_id
    out_dir.mkdir(parents=True, exist_ok=True)
    csvs_dir = out_dir / "csvs"
    csvs_dir.mkdir(parents=True, exist_ok=True)
    
    # Write the deg2tf_mapping CSV as before.
    out_csv = csvs_dir / "deg2tf_mapping.csv"
    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv, index=False)
    logging.info(f"[tffetcher] => Written deg2tf_mapping CSV with {len(df_out)} rows to {out_csv}")

    # Call the enrichment analysis and capture its returned DataFrame.
    df_enrich = tfenrichment.run_enrichment(merged_enrichment_network, rows, filtered_degs, out_dir, params_conf, nuc_reg=NUCLEOID_REGULATORS)

    # Save the TF enrichment summary CSV.
    enrich_csv_path = csvs_dir / "tf_enrichment_summary.csv"
    df_enrich.to_csv(enrich_csv_path, index=False)
    logging.info(f"[tfenrichment] TF enrichment summary saved: {enrich_csv_path}")


# --- Helper functions ---

def normalize_gene(gene: str) -> str:
    import re
    return re.sub(r'[^a-z0-9]', '', gene.lower().strip())

def collect_degs_with_sources(csv_dir: Path, deg_includes: list = None) -> Dict[str, Set[str]]:
    from collections import defaultdict
    gene_to_sources: Dict[str, Set[str]] = defaultdict(set)
    matching_files = []
    if deg_includes:
        for pattern in deg_includes:
            matched = glob.glob(str(csv_dir / pattern))
            matching_files.extend(matched)
        if not matching_files:
            logging.warning(f"[tffetcher] deg_csv_includes was {deg_includes}, but no files matched!")
    else:
        matching_files = glob.glob(str(csv_dir / "*.csv"))
    for csv_fp in sorted(set(matching_files)):
        filename = os.path.basename(csv_fp)
        deg_label = derive_deg_label(filename)
        df = pd.read_csv(csv_fp)
        if "gene" not in df.columns:
            logging.warning(f"[tffetcher] Skipping {filename}, no 'gene' column.")
            continue
        for gene_val in df["gene"].dropna():
            g = normalize_gene(gene_val)
            gene_to_sources[g].add(deg_label)
    return dict(gene_to_sources)

def get_num_files_considered(csv_dir: Path, deg_includes: list) -> int:
    import glob
    if deg_includes:
        all_files = []
        for pattern in deg_includes:
            all_files.extend(glob.glob(str(csv_dir / pattern)))
        return len(set(all_files))
    else:
        return len(glob.glob(str(csv_dir / "*.csv")))

def union_networks(net_a: Dict[str, Dict[str, Dict[str, Set[str]]]],
                   net_b: Dict[str, Dict[str, Dict[str, Set[str]]]]) -> Dict[str, Dict[str, Dict[str, Set[str]]]]:
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
                merged[gene_b][reg_b] = {"polarity_set": set(), "sources": set()}
            merged[gene_b][reg_b]["polarity_set"].update(info_b["polarity_set"])
            merged[gene_b][reg_b]["sources"].update(info_b["sources"])
    return merged

def intersect_networks(net_a: Dict[str, Dict[str, Dict[str, Set[str]]]],
                       net_b: Dict[str, Dict[str, Dict[str, Set[str]]]]) -> Dict[str, Dict[str, Dict[str, Set[str]]]]:
    merged = {}
    common_genes = set(net_a.keys()).intersection(set(net_b.keys()))
    for gene in common_genes:
        merged[gene] = {}
        regs_a = net_a[gene]
        regs_b = net_b[gene]
        common_regs = set(regs_a.keys()).intersection(set(regs_b.keys()))
        for reg in common_regs:
            pol_set = set(regs_a[reg]["polarity_set"]) | set(regs_b[reg]["sources"])
            src_set = set(regs_a[reg]["sources"]) | set(regs_b[reg]["sources"])
            merged[gene][reg] = {"polarity_set": pol_set, "sources": src_set}
    return merged

def unify_polarities(polarity_set: Set[str]) -> str:
    pol_clean = {p for p in polarity_set if p not in ["NA", None, ""]}
    if not pol_clean:
        return "NA"
    if len(pol_clean) == 1:
        return pol_clean.pop()
    return "±"

def derive_short_network_label(net_name: str) -> str:
    net_name = net_name.lower().strip()
    if "ecocyc" in net_name:
        return "ecocyc_28"
    elif "regulon" in net_name:
        return "regdb_13"
    return net_name

def derive_deg_label(filename: str) -> str:
    base = filename.lower().replace('.csv', '')
    if 'upregulated' in base:
        return base.replace('_upregulated_degs', '_up')
    elif 'downregulated' in base:
        return base.replace('_downregulated_degs', '_down')
    else:
        return base