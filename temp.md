I would like you to create new python modules (vazulka, kim, and zhang) that abide by pragmatic programming philosophy best practices, where the code is easier to change, design by contract, and intuitive with assertive programming.

The deliverable is to create new python modules that read data from a local repository, where the relative paths are provided in utils.py (which is below). Unlike the modules that are already developed (e.g., like sanchez_vasquez.py pasted below), some these new datasets need not generate plots or produce filtered datasets based on defined thresholds (but vazulka does need IQR). Rather, these genes are already labelled as either up or down regulated. For example, consider the below context:

"Comparative Omics (RNA‐seq, Absolute Proteomics)
Most of these datasets compare omics readouts between a single “target” and “reference” condition, enabling the identification of up- and down-regulated genes. Some studies provide full raw data, which lets us reproduce results, apply custom thresholds, and isolate differentially expressed genes ourselves (as is done in deg2tfbs). Other articles do not share raw data but instead list up- and down-regulated genes directly, in which case we simply import those gene sets into deg2tfbs to identify the associated transcription factors and their DNA-binding sites. These binding sites can then later be used in dnadesign."

We are currently creating modules in the dataloader directory within def2tfbs. Here is a readme:


### deg2tfbs

A pipeline for processing RNA-seq and proteomic datasets from *E. coli*. This tool identifies differentially expressed genes (DEGs) between conditions, maps them to transcription factors using RegulonDB, and retrieves their corresponding transcription factor binding sites (TFBSs).  

The pipeline steps (DEG analysis → TF mapping → TFBS identification).

### Usage
Obtain the raw data (e.g., from **dnadesign-data**) and place it so that the ```utils.py``` dictionary can locate each dataset via ```DATA_FILES[...]```.

Configure ```configs/example.yaml```:
- Adjust ```output.root_dir``` to point to ```outputs```.
- Change ```output.batch_identifier``` to e.g. ```batch_2023_08_02```.
- Modify each dataset’s ```data``` and ```thresholds``` as desired.

Run from the command line 
```bash
cd deg2tfbs
python main.py
```

###  Project Layout
```
deg2tfbs/
├── __init__.py
├── pipeline/
│   ├── __init__.py
│   ├── deg_analysis.py          # Step 1: Identify DEGs (WIP)
│   ├── tf_fetcher.py           # Step 2: Map DEGs to TFs (WIP)
│   ├── tfbs_fetcher.py         # Step 3: Identify TFBS (WIP)
│   └── dataloader/
│       ├── ceroni.py
│       ├── mori.py
│       ├── wu.py
│       ├── zhu.py
│       ├── emani.py
│       ├── schmidt.py
│       └── radzikowski.py
├── configs/
│   └── example.yaml
├── outputs/
│   └── batch_2023_08_02/ 
│       ├── csv
│       └── plots
├── main.py
└── README.md

```


### Data Source
Requires source data from [dnadesign-data](https://github.com/e-south/dnadesign-data). Update the config file to point to the correct data paths.

- **Data Source**: Utilize raw files from or another local copy.  
- **Intermediate Outputs**:  
  1. **DEGs** (CSV)  
  2. **DEGs + TF associations** (CSV)  
  3. **TFBS** for those TFs (CSV)  

## Installation

```bash
pip install -r requirements.txt

```

---

Here is some information about the datasets to create modules on (for example use this to populate top of module doc strings, but create this in the format of ceroni.py). Note that each csv outputted from each module should either/both <module_name>_upregulated_degs.csv, <module_name>_downregulated_degs.csv. Some of the datasets will only need to return one of these types of csvs, though, and that information is indicated below.

Required columns upon saving to csv are required_columns = ["gene", "source", "comparison"]
source is the all lower case name of the source, e.g., Vazulka is "vazulka", comparison is listed below

---

Vazulka et al.
Title: RNA-seq reveals multifaceted gene expression response to Fab production in Escherichia coli fed-batch processes with particular focus on ribosome stalling
DOI: 10.1186/s12934-023-02278-w
Association: Fab production
Comments: Characterized the gene expression response in E. coli BL21(DE3) and HMS174(DE3) to periplasmic Fab expression via fed-batch RNA‐seq.
Format: List of up- and down-regulated genes available.
TableS1.xlsx
sheet: 'Sheet1'
header=0
columns: 'Gene'	'baseMean'	'log2FoldChange'	'lfcSE'	'pvalue'	'padj'
perform IQR plot and thresholding on 'log2FoldChange' column
comparison: 'fab_production_2h_versus_control'

---

Kim et al.
Title: Heat-responsive and time-resolved transcriptome and metabolome analyses of Escherichia coli uncover thermo-tolerant mechanisms
DOI: 10.1038/s41598-020-74606-8
Association: Heat shock response
Comments: Applied RNA‐seq to capture early, middle, and late stages of heat stress (2 min–40 h), illuminating initiation, adaptation, and phenotypic plasticity phases in E. coli.
Format: List of up- and down-regulated genes available.
41598_2020_74606_MOESM1_ESM.xlsx
sheet: 'Table_S1
header=5
threshold=2.5
upregulated genes beyond threshold are red, downregulated are green
columns needed are 'Gene' (convert to 'gene' before saving as csv), '-30min', '1h'
comparison: '42C_versus_control'


---

Zhang et al.
Title: Heat-Shock Response Transcriptional Program Enables High-Yield and High-Quality Recombinant Protein Production in Escherichia coli
DOI: 10.1021/cb5004477
Association: Heat shock response
Comments: Demonstrated that a σ^32‐I54N HSR-like reprogrammed proteostasis network can boost soluble, folded, and functional recombinant proteins
Format: List of up-regulated genes (no down-) available.
TableS1.xlsx
sheet: 'Sheet1'
header: 0
columns: 'Gene'	'σ32 (I54N)'	'SD'	'σ32 (WT)'	'SD'	'Heat (42ºC)'	'SD'
change 'Gene' to 'gene' before saving as csv
select upregulated genes by isolating the 'σ32 (I54N)' row
there are no downregulated genes to save for this module
comparison: 'σ32-I54N_expression_versus_control'

---

here are the paths in utils.py:
"""
--------------------------------------------------------------------------------
<deg2tfbs project>
utils.py

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
from pathlib import Path
import pandas as pd

DNADESIGN_DATA = Path(__file__).parent.parent.parent.parent / 'dnadesign-data'

# A simple dictionary mapping keys to actual file paths:
DATA_FILES = {
    # Comparative Omics (RNA‐seq, Absolute Proteomics) Datasets
    # Most of these datasets compare omics readouts between a single “target” and “reference” condition, 
    # enabling the identification of up- and down-regulated genes. Some studies provide full raw data, 
    # which lets us reproduce results, apply custom thresholds, and isolate differentially expressed 
    # genes ourselves

    # Diverse Media Conditions
    'mori_ev8': DNADESIGN_DATA / 'primary_literature' /'Mori_et_al' / 'msb20209536-sup-0009-datasetev8.xlsx',
    'mori_ev9': DNADESIGN_DATA / 'primary_literature' /'Mori_et_al' / 'msb20209536-sup-0010-datasetev9.xlsx',
    'schmidt': DNADESIGN_DATA / 'primary_literature' / 'Schmidt_et_al' / 'table_s6.xlsx',
    
    # Stringent Response and Nutrient Limitation
    'durfee': DNADESIGN_DATA / 'primary_literature' / 'Durfee_et_al' / 'ST1_durfee_degs.xlsx',
    'fragoso_et_al': DNADESIGN_DATA / 'primary_literature' / 'Fragoso-Jimenez_et_al' / '12934_2022_1909_MOESM_ESM.xlsx',
    'franchini_a': DNADESIGN_DATA / 'primary_literature' / 'Franchini_et_al_a' / 'gene_list.xlsx',
    'franchini_b': DNADESIGN_DATA / 'primary_literature' / 'Franchini_et_al_b' / 'gene_list.xlsx',
    'gummesson': DNADESIGN_DATA / 'primary_literature' / 'Gummesson_et_al' / 'DataSheet_3.xlsx',
    'houser': DNADESIGN_DATA / 'primary_literature' / 'Houser_et_al' / 'pcbi.1004400.s004.csv',
    'lu': DNADESIGN_DATA / 'primary_literature' / 'Lu_et_al' / 'ST1_ST2_glyphosate_shock.csv',
    'sanchez_vasquez': DNADESIGN_DATA /'primary_literature' / 'Sanchez-Vazquez_et_al' / 'pnas.1819682116.sd01.xlsx',
    'wu': DNADESIGN_DATA / 'primary_literature' / 'Wu_et_al' / '41564_2022_1310_MOESM3_ESM.xlsx',
    'zhu': DNADESIGN_DATA / 'primary_literature' / 'Zhu_et_al' / '41467_2023_36254_MOESM3_ESM.xlsx',

    # Metabolic Burden
    'ceroni': DNADESIGN_DATA / 'primary_literature' / 'Ceroni_et_al' / '41592_2018_BFnmeth4635_MOESM4_ESM.xlsx',
    'rajacharya': DNADESIGN_DATA / 'primary_literature' / 'Rajacharya_et_al' / '41598_2024_63148_MOESM2_ESM.xlsx',
    
    # Membrane Stress and Fatty Acid Production
    'emani': DNADESIGN_DATA / 'primary_literature' / 'Emani_et_al' / 'Supplementary data.xlsx',
    'vazulka': DNADESIGN_DATA / 'primary_literature' / 'Vazulka_et_al' / 'TableS1.xlsx',
    
    # Antibiotic Stress
    'bie': DNADESIGN_DATA / 'primary_literature' / 'Bie_et_al' / 'spectrum.00317-23-s0002.xlsx',
    'deter': DNADESIGN_DATA / 'primary_literature' / 'Deter_et_al' / 'GSE156896_read_counts.csv',
    'radzikowski': DNADESIGN_DATA / 'primary_literature' / 'Radzikowski_et_al' / 'msb166998-sup-0002-tableev1.xlsx',
    
    # Heat Shock Response
    'kim': DNADESIGN_DATA / 'primary_literature' / 'Kim_et_al' / '41598_2020_74606_MOESM1_ESM.xlsx',
    'zhang': DNADESIGN_DATA / 'primary_literature' / 'Zhang_et_al' / 'TableS1.xlsx',
    
    # Phage Shock Response
    'jovanovic': DNADESIGN_DATA / 'primary_literature' / 'Jovanovic_et_al' / 'Supplementary_table_4_MG1655+PspA.xls',
    
    # Other Curated Literature
    'peebo': DNADESIGN_DATA / 'primary_literature' / 'Peebo_et_al' / 'gene_list.xlsx',

    # PRECISE-1K datasets from imodulonDB (gene sets following data-driven regulon categorization) (E. coli)
    'lamoureux_1k': DNADESIGN_DATA / 'primary_literature' / 'Lamoureux_et_al' / 'PRECISE_1K',
        
    # EcoCyc Pathway/Genome Database Full Regulatory Network
    'ecocyc_full_reg_network': DNADESIGN_DATA / 'EcoCyc' / 'ECOLI-regulatory-network.txt',
    'ecocyc_genes': DNADESIGN_DATA / 'RegulonDB_11' /'genes' / 'ecocyc_genes_and_synonyms.txt',
    
    # RegulonDB data sets
    'k_12_genome': DNADESIGN_DATA / 'RegulonDB_11' / 'K12_genome' / 'E_coli_K12_MG1655_U00096.3.txt',
    'regulondb_tf_set': DNADESIGN_DATA / 'RegulonDB_11' / 'tf_factors' / 'TFSet.csv',
    'regulondb_growth_condition_set': DNADESIGN_DATA / 'RegulonDB_11' /  'network_associations' / 'GCSet.txt',
    'regulondb_promoter_set': DNADESIGN_DATA / 'RegulonDB_11' / 'promoters' / 'PromoterSet.csv',
    'regulondb_promoter_predicted_set': DNADESIGN_DATA / 'RegulonDB_11' / 'promoters' / 'PromoterPredictionSet.csv',
    'regulondb_promoter_mendoza_race_set': DNADESIGN_DATA / 'RegulonDB_11' / 'promoters' / 'Promoter_from_RACE_Dataset.csv',
    'regulondb_promoter_mendoza_pyroseq_set': DNADESIGN_DATA / 'RegulonDB_11' / 'promoters' / 'Promoter_from_454_Dataset.csv',
    'regulondb_promoter_salgado_set': DNADESIGN_DATA / 'RegulonDB_11' / 'promoters' / 'ht_transcription_initiation_mapping_with_5_tri_or_monophosphate_enrichment_v3.0.csv',
    'regulondb_promoter_hernandez_set': DNADESIGN_DATA / 'primary_literature' / 'Hernandez et al' / 'positive2860.txt',
    'regulondb_non_promoter_hernandez_set': DNADESIGN_DATA / 'primary_literature' / 'Hernandez et al' / 'negative2860.txt',

    'regulondb_tfbs_prediction_medina_set': DNADESIGN_DATA / 'RegulonDB_11' / 'tf_binding_sites' / 'TFBSs_predictions_v3.txt',
    'regulondb_tfbs_prediction_hernandez_set': DNADESIGN_DATA / 'RegulonDB_11' / 'tf_binding_sites' / 'BindingSitePredictionSet.txt',
    'regulondb_tfbs_santos_set': DNADESIGN_DATA / 'RegulonDB_11' / 'tf_binding_sites' / 'BindingSiteSet.txt',
    'regulondb_tfbs_PSSM_set': DNADESIGN_DATA / 'RegulonDB_11' / 'tf_binding_sites' / 'RegulonDB_PSSM_v4.0' / 'results' / 'PSSM-Dataset-v4.0.txt',
    'regulondb_regulatory_interactions_set': DNADESIGN_DATA / 'RegulonDB_11' / 'network_associations' / 'RISet.txt',
    'regulondb_network_associations_tf_tf': DNADESIGN_DATA / 'RegulonDB_11' / 'network_associations' / 'network_tf_tf.txt',
    'regulondb_network_associations_tf_gene': DNADESIGN_DATA / 'RegulonDB_11' / 'network_associations' / 'network_tf_gene.txt',
    "target_set": DNADESIGN_DATA.parent.parent / 'decoydesigner' / 'data' / 'target_sets' / 'stress_and_growth_associated_tfs_20230603.csv',
    
    # High-throughput, functional promoter characterization datasets (E. coli)
    'thomason_tss_set': DNADESIGN_DATA / 'primary_literature' / 'Thomason_et_al' / 'zjb999093409sd1.xlsx',
    'johns_sequences': DNADESIGN_DATA / 'primary_literature' / 'Johns_et_al' / '41592_2018_BFnmeth4633_MOESM3_ESM.xlsx',
    'johns_expression_metrics': DNADESIGN_DATA / 'primary_literature' / 'Johns_et_al' / '41592_2018_BFnmeth4633_MOESM4_ESM.xlsx',
    'johns_expression_across_conditions': DNADESIGN_DATA / 'primary_literature' / 'Johns_et_al' / '41592_2018_BFnmeth4633_MOESM5_ESM.xlsx',
    'sun_yim': DNADESIGN_DATA / 'primary_literature' / 'Sun_Yim_et_al' / 'msb198875-sup-0002-sdatafig1.xlsx',
}


def load_dataset(dataset_key, sheet_name=None, usecols=None, header=0, skiprows=None):
    """
    Minimal example that loads a dataset from the known path in DATA_FILES.
    If it's an Excel, we read via pd.read_excel; if CSV, pd.read_csv, etc.
    """
    file_path = DATA_FILES[dataset_key]
    if not file_path.exists():
        raise FileNotFoundError(f"File not found for key {dataset_key}: {file_path}")

    # Fix: Treat both .xlsx and .xls as Excel
    if file_path.suffix in [".xlsx", ".xls"]:
        df = pd.read_excel(
            file_path, 
            sheet_name=sheet_name,
            usecols=usecols,
            header=header,
            skiprows=skiprows
        )
    elif file_path.suffix == ".csv":
        df = pd.read_csv(file_path, usecols=usecols, header=header)
    else:
        raise ValueError(f"Unsupported file format for {file_path.suffix}")

    return df

I also want you to update the config data file:
# configs/example.yaml

output:
  root_dir: "pipeline/deg_fetcher"
  batch_identifier: "batch_20250130"

# Ceroni et al. (DOI: 10.1038/nmeth.4635)
ceroni:
  data:
    dataset_key: "ceroni"
    sheet_name: "Mean FPKM"
    usecols: "A:P"
    header: 0
    skiprows: 3

  thresholds:
    log2_fc_threshold: 2.5
    time_point: 15
    drop_zeros: true
    plasmid_pairs:
      - [ "pLys-M1", "pLys" ]
      - [ "pPD864-LacZ", "pPD864" ]
      - [ "pSB1C3-H3", "pSB1C3" ]

  output:
    csv_subdir: "csv"
    plot_subdir: "plots"

# Sanchez‐Vazquez et al. DOI: 10.1073/pnas.1819682116
sanchez_vasquez:
  data:
    dataset_key: "sanchez_vasquez"
    sheet_name: "Data"
    usecols:
      - "Gene"
      - "1+2+ 5 min"
      - "1+2+ 5 min Category"
      - "1+2+ 10 min"
      - "1+2+ 10 min Category"
    header: 0
  thresholds:
    fc_col: "1+2+ 5 min" 
  output:
    csv_subdir: "csv"
    plot_subdir: "plots"

...


---
here is sanchez_vasquez.py for reference, in terms of part of my requires to have a load function and a pipeline function associated with each module:
"""
--------------------------------------------------------------------------------
<deg2tfbs project>
sanchez_vasquez.py

Module for reading Sanchez‐Vazquez et al., which performed RNA‐seq of E. coli 
with and without ppGpp‐binding sites on RNAP, characterizing gene expression 
changes at 5–10 min.

The module identifies extreme up- and down-regulated genes using IQR-based outlier 
detection applied to the '1+2+ 5 min' column from the RNA-deg dataset.

"Genome-wide effects on Escherichia coli transcription from ppGpp binding to 
its two sites on RNA polymerase"
DOI: 10.1073/pnas.1819682116

Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
from pathlib import Path

import yaml
import pandas as pd
import numpy as np  # Add this import
import seaborn as sns
import matplotlib.pyplot as plt

from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_sanchez_vasquez_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Sanchez Vasquez dataset using keys from the YAML config.
    
    Expected columns: 'Gene', '1+2+ 5 min Category', '1+2+ 10 min Category', etc.
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data.get("sheet_name"),
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    return df


def sanchez_distribution_plot(
    df: pd.DataFrame,
    log2_fc_col: str,
    up_mask,  # Change from up_idx to up_mask
    down_mask,  # Change from down_idx to down_mask
    output_path: Path
):
    sns.set_style("ticks")
    plt.figure(figsize=(5,7))

    # Create color map using direct boolean masks
    color_map = np.select(
        [up_mask, down_mask],
        ["red", "green"],
        default="gray"
    )

    xvals = np.random.uniform(-0.2, 0.2, size=len(df))
    plt.scatter(xvals, df[log2_fc_col], c=color_map, alpha=0.5, linewidth=0.5)

    plt.title("Sanchez-Vazquez et al.\nLog2 Fold Change (5 min) with IQR Filtering")
    plt.ylabel("log2(Fold Change)")
    plt.xticks([], [])
    
    sns.despine(top=True, right=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()



def run_sanchez_vasquez_pipeline(full_config: dict):
    """
    Reads the data, 
    Up => '1+2+ 5 min Category' == 'B'
    Down => '1+2+ 5 min Category' == 'A'

    comparison => "relA_overexpression_versus_control"
    source => "sanchez_vasquez"
    """
    config_sv = full_config.get("sanchez_vasquez", None)
    if config_sv is None:
        print("[Sanchez-Vazquez Pipeline] No 'sanchez_vasquez' config found. Skipping.")
        return

    # Load
    df = read_sanchez_vasquez_data(config_sv["data"])

    assert "Gene" in df.columns, "Expected 'Gene' in Sanchez-Vazquez data"
    cat_col = "1+2+ 5 min Category"
    assert cat_col in df.columns, f"Expected '{cat_col}' in Sanchez-Vazquez data"

    # Build output path
    project_root = Path(__file__).parent.parent.parent
    output_root = project_root / full_config["output"]["root_dir"]
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    csv_dir = batch_dir / config_sv["output"]["csv_subdir"]
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = batch_dir / config_sv["output"]["plot_subdir"]
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Plot distribution
    log2_fc_col = "1+2+ 5 min"  # Verify column name matches your data
    df[log2_fc_col] = pd.to_numeric(df[log2_fc_col], errors='coerce')
    df = df.dropna(subset=[log2_fc_col, cat_col])
    
    # Calculate IQR bounds for outlier detection
    q1, q3 = df[log2_fc_col].quantile([0.25, 0.75])
    iqr = q3 - q1
    upper_bound = q3 + 1.5*iqr
    lower_bound = q1 - 1.5*iqr
    
    up_filter = (df[log2_fc_col] > upper_bound)
    down_filter = (df[log2_fc_col] < lower_bound)
    df_up = df[up_filter].copy()
    df_down = df[down_filter].copy()

    # Plot
    plot_path = plot_dir / "sanchez_vasquez_iqr.png"
    sanchez_distribution_plot(
        df=df,
        log2_fc_col=log2_fc_col,
        up_mask=up_filter,  # Pass the boolean mask directly
        down_mask=down_filter,  # Pass the boolean mask directly
        output_path=plot_path
    )
    
    # Validation asserts
    assert up_filter.sum() > 0, "No upregulated genes passing both IQR and category filters"
    assert down_filter.sum() > 0, "No downregulated genes passing both filters"
    
    up_clean = pd.DataFrame({
        "gene": df_up["Gene"],
        "source": "sanchez_vasquez",
        "comparison": "relA_overexpression_versus_control"
    }).drop_duplicates()

    down_clean = pd.DataFrame({
        "gene": df_down["Gene"],
        "source": "sanchez_vasquez",
        "comparison": "relA_overexpression_versus_control"
    }).drop_duplicates()

    # Save
    up_csv = csv_dir / "sanchez_vasquez_upregulated_degs.csv"
    down_csv = csv_dir / "sanchez_vasquez_downregulated_degs.csv"
    up_clean.to_csv(up_csv, index=False)
    down_clean.to_csv(down_csv, index=False)

    print(f"[Sanchez-Vazquez et al. Pipeline] Completed. Gathered list of reported DEGs across 1 condition pair with IQR-based outlier detection: {len(up_clean)} up, {len(down_clean)} down.")

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_sanchez_vasquez_pipeline(full_config)


---
reference this ceroni.py plot for plotting the thresholds for Kim et al.

"""
--------------------------------------------------------------------------------
<deg2tfbs project>
ceroni.py

Module for loading and analyzing data described in Ceroni et al., which 
combined RNA-seq with in vivo assays to identify transcriptional changes that
occur in E. coli when burdensome synthetic constructs are expressed.

The module isolates up- and down-regulated genes based on a user-defined log2 fold 
change threshold and saves an MA plot (average expression vs. log2 Fold Change).

"Burden-driven feedback control of gene expression"
DOI: 10.1038/nmeth.4635

Module Author(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

import os
from pathlib import Path

import yaml
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from deg2tfbs.pipeline.dataloader.utils import load_dataset

def read_ceroni_data(config_data: dict) -> pd.DataFrame:
    """
    Reads the sourced Ceroni dataset using keys from the YAML config.
    
    E.g. config:
      {
        "dataset_key": "ceroni",
        "sheet_name": "Mean FPKM",
        "usecols": "A:P",
        "header": 0,
        "skiprows": 3
      }
    """
    df = load_dataset(
        dataset_key=config_data["dataset_key"],
        sheet_name=config_data["sheet_name"],
        usecols=config_data.get("usecols"),
        header=config_data.get("header", 0),
        skiprows=config_data.get("skiprows", None)
    )
    
    # Specify which plasmids (i.e., strains) and time points we're interested in
    plasmid_names = ['pSB1C3-H3','pLys-M1','pSB1C3','pLys','pSB1C3-Lux','pPD864-LacZ','pPD864']
    time_points = ['15','60'] # minutes
    headers = ['gene_id','gene_name'] + [f"{p}_{t}" for p in plasmid_names for t in time_points]
    df.columns = headers
    
    # Tidy up the DataFrame
    df_melted = pd.melt(
        df,
        id_vars=['gene_id', 'gene_name'],
        var_name='Plasmid_Time',
        value_name='DE_of_Mean_FPKM'
    )
    df_melted[['Plasmid','Time']] = df_melted['Plasmid_Time'].str.rsplit('_', n=1, expand=True)
    df_melted['Time'] = df_melted['Time'].astype(int)

    return df_melted


def ceroni_ma_plot(
    df,
    plasmid1,
    plasmid2,
    threshold=2.0,
    epsilon=1e-5,
    plot_path=None,
    time_point=None
):
    """
    Creates and saves an MA-plot comparing strains carrying plasmid1 
    and plasmid2 from the Ceroni dataset.
    """
    # Set up the plot
    sns.set_theme(style="ticks")
    df['log2_ratio'] = np.log2((df[plasmid1] + epsilon) / (df[plasmid2] + epsilon))
    df['average_expression'] = (df[plasmid1] + df[plasmid2]) / 2

    # Up- and down-regulated genes are colored red and green, respectively
    colors = np.where(
        df['log2_ratio'] >= threshold, 'red',
        np.where(df['log2_ratio'] <= -threshold, 'green', 'lightgray')
    )
    # Plot 
    plt.figure(figsize=(6,5))
    plt.scatter(df['average_expression'], df['log2_ratio'], c=colors, alpha=0.25, edgecolors='none')
    plt.xscale('log')
    plt.axhline(threshold, color='gray', linestyle='--')
    plt.axhline(-threshold, color='gray', linestyle='--')

    # Title and subtitle
    if time_point is not None:
        plt.title(f"Ceroni et al.\nlog2({plasmid1} / {plasmid2}); {time_point} mins")
    else:
        plt.title(f"Ceroni et al.\nlog2({plasmid1} / {plasmid2})")

    plt.xlabel("Average Expression")
    plt.ylabel("Log2 Fold Change")
    sns.despine()

    if plot_path:
        plt.savefig(plot_path, dpi=150)
    plt.close()


def ceroni_thresholding(df, plasmid1, plasmid2, threshold=2.0, time_point=15, drop_zeros=False):
    """
    For a given pair of strains carrying different plasmids, filter df by time_point, 
    pivot the data, compute the log₂ ratio between average FPKM values, and classify 
    genes as upregulated or downregulated.

    Returns two DataFrames: upregulated and downregulated.
    """
    df_sub = df[(df['Time'] == time_point) & (df['Plasmid'].isin([plasmid1, plasmid2]))].copy()
    pivoted = df_sub.pivot_table(
        index=['gene_id','gene_name'], 
        columns='Plasmid',
        values='DE_of_Mean_FPKM'
    )

    # Drop rows with zeros in either plasmid
    if drop_zeros:
        pivoted = pivoted[(pivoted[plasmid1]!=0) & (pivoted[plasmid2]!=0)]

    # Compute log2 ratio and average expression
    pivoted['log2_ratio'] = np.log2((pivoted[plasmid1]+1e-5)/(pivoted[plasmid2]+1e-5))
    pivoted['average_expression'] = (pivoted[plasmid1]+pivoted[plasmid2])/2
    pivoted.reset_index(inplace=True)

    # Figure out upregulated, downregulated
    up = pivoted[pivoted['log2_ratio'] >= threshold]
    down = pivoted[pivoted['log2_ratio'] <= -threshold]

    return up, down, pivoted


def run_ceroni_pipeline(full_config: dict):
    """
    full_config is the entire YAML content, containing:
      - output: { root_dir: ..., batch_identifier: ... }
      - ceroni: { data, thresholds, output subdirs, etc. }
    """

    # Extract the "ceroni" config subset
    config_ceroni = full_config["ceroni"]

    # Load data
    df_melted = read_ceroni_data(config_ceroni["data"])

    # Build the final output paths
    project_root = Path(__file__).parent.parent.parent
    
    # Output folder, e.g. deg2tfbs/pipeline/deg_fetcher
    output_root = project_root / full_config["output"]["root_dir"]

    # Append the batch identifier
    batch_id = full_config["output"]["batch_identifier"]
    batch_dir = output_root / batch_id

    # Define subfolders for "csv" and "plots" values from config
    csv_dir = batch_dir / config_ceroni["output"]["csv_subdir"]
    plot_dir = batch_dir / config_ceroni["output"]["plot_subdir"]

    # Create the directories if they don't exist
    csv_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Retrieve thresholds
    plasmid_pairs = config_ceroni["thresholds"]["plasmid_pairs"]
    threshold = config_ceroni["thresholds"]["log2_fc_threshold"]
    time_point = config_ceroni["thresholds"]["time_point"]
    drop_zeros = config_ceroni["thresholds"]["drop_zeros"]
    time_point = config_ceroni["thresholds"]["time_point"]  # e.g. 15, 60

    all_up = []
    all_down = []
    
    for (p1, p2) in plasmid_pairs:
        up, down, pivoted = ceroni_thresholding(
            df_melted,
            p1, p2,
            threshold=threshold,
            time_point=time_point,
            drop_zeros=drop_zeros
        )

        required_columns = ["gene", "source", "thresholds", "comparison"]
        
        # Tidy up the DataFrames of up- and down-regulated genes
        up_clean = up.copy()
        up_clean["gene"] = up_clean["gene_name"]
        up_clean["source"] = "ceroni"
        up_clean["thresholds"] = threshold
        up_clean["comparison"] = f"{p1}_versus_{p2}"
        up_clean = up_clean[required_columns]
        all_up.append(up_clean)

        down_clean = down.copy()
        down_clean["gene"] = down_clean["gene_name"] 
        down_clean["source"] = "ceroni"
        down_clean["thresholds"] = threshold
        down_clean["comparison"] = f"{p1}_versus_{p2}"
        down_clean = down_clean[required_columns]
        all_down.append(down_clean)

        # Save MA plot: fold changes (M values) against average expression (A values)
        plot_path = plot_dir / f"ceroni_{p1}_versus_{p2}.png"
        ceroni_ma_plot(pivoted, p1, p2, threshold=threshold, plot_path=plot_path)

    # Combine all up and down DEGs
    up_all = pd.concat(all_up, ignore_index=True).drop_duplicates()
    down_all = pd.concat(all_down, ignore_index=True).drop_duplicates()

    up_all.to_csv(csv_dir / "ceroni_upregulated_degs.csv", index=False)
    down_all.to_csv(csv_dir / "ceroni_downregulated_degs.csv", index=False)

    print(f"[Ceroni et al. Pipeline] Completed. Identified DEGs across {len(plasmid_pairs)} condition pairs at log2 ≥ {threshold}: {len(up_all)} up, {len(down_all)} down.")  

if __name__ == "__main__":
    config_path = Path(__file__).parent.parent.parent / "configs" / "example.yaml"
    with open(config_path, "r") as f:
        full_config = yaml.safe_load(f)

    run_ceroni_pipeline(full_config)
