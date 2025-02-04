## **deg2tfbs**

**deg2tfbs** is a pipeline to derive transcription factor binding sites (TFBSs) from comparative RNA-seq and proteomic data, focusing on *E. coli* but extendable to other organisms. The process identifies differentially expressed genes (DEGs) between experimental conditions, maps them to upstream transcription factors (TFs) using [**RegulonDB**](https://regulondb.ccg.unam.mx/) or [**EcoCyc**](https://ecocyc.org/) databases, and then retrieves cognate transcription factor binding sites (TFBSs) if available.  

![Alt text](images/pipeline.png)

## **Installation**

You can manage your Python dependencies using [conda](https://docs.conda.io/) (or [mamba](https://mamba.readthedocs.io/), a fast drop-in replacement for conda). 

**Using Conda/Mamba:**

1. **Create and Activate an Environment:**

   ```bash
   conda create -n deg2tfbs_env python
   conda activate deg2tfbs_env
   ```

2. **Install Dependencies:**

   If you are using conda, you can install packages from the `requirements.txt` or use pip:
   
   ```bash
   cd deg2tfbs
   conda install --file requirements.txt
   # or
   pip install -r requirements.txt
   ```

## **Pipeline Steps**

1. **degfetcher** *(Step 1: Isolate DEGs)*  
   - Loads comparative omics datasets from the [**dnadesign-dna**](https://github.com/e-south/dnadesign-data) repository.
   - Produces tidy CSV outputs, such as `ceroni_upregulated_degs.csv`, containing columns:  
     - **gene**: DEG identifier.
     - **source**: The source dataset (e.g., "ceroni").
     - **comparison**: The experimental context defining the target vs. reference condition.

      For example:
      | gene  | source  | thresholds | comparison                  |
      |-------|--------|------------|------------------------------|
      | groS  | ceroni | 2.5        | pLys-M1_versus_pLys          |
      | gadY  | ceroni | 2.5        | pLys-M1_versus_pLys          |
      | azuC  | ceroni | 2.5        | pLys-M1_versus_pLys          |
      | yadL  | ceroni | 2.5        | pPD864-LacZ_versus_pPD864    |
      | rclC  | ceroni | 2.5        | pPD864-LacZ_versus_pPD864    |
      | rclR  | ceroni | 2.5        | pPD864-LacZ_versus_pPD864    |


2. **tffetcher** *(Step 2: Map DEGs to TFs)*  
   - Reads the CSV outputs, saved in batches, from **degfetcher**.
   - Loads **regulatory network**-type datasets (e.g., from **RegulonDB** or **EcoCyc**).  
   - Fetches TFs that reportedly regulate these DEGs (the “regulatees”).  
   - Produces a tidy CSV output, `deg2tf_mapping`, containing columns: 
     - **gene** - DEG identifier.
     - **regulator**: TFs reported to regulated target gene.
     - **polarity**: Reported "type" of regulation imposed, if available.
     - **source**: Record of which regulatory interaction datasets were used.
     - **is_global_regulator**: Boolean indicating whether the regulator is classified as a global regulator by EcoCyc.
     - **is_sigma_factor**: Boolean indicating whether the regulator is classified as a sigma factor by EcoCyc.
     - **deg_source**: Indicates the source datasets in which this gene appears, as processed by degfetcher.

      For example:
      | gene  | regulator | polarity | source                     | is_global_regulator | is_sigma_factor | DEG Source       |
      |-------|-----------|----------|----------------------------|---------------------|-----------------|------------------|
      | aaea  | crp       | +        | ecocyc_28_AND_regdb_13     | yes                 | no              | houser_up        |
      | aaea  | aaer      | +        | ecocyc_28                  | no                  | no              | houser_up        |
      | aaeb  | crp       | +        | ecocyc_28_AND_regdb_13     | yes                 | no              | houser_up        |
      | abga  | nac       | -        | ecocyc_28_AND_regdb_13     | yes                 | no              | ceroni_up        |
      | acca  | accd      | -        | ecocyc_28                  | no                  | no              | houser_down      |
      | acca  | rpod      | +        | ecocyc_28                  | no                  | yes             | houser_down      |
      | adia  | adiy      | +        | ecocyc_28_AND_regdb_13     | no                  | no              | houser_up-lu_up  |


3. **tfbsfetcher** *(Step 3: Map TFs to TFBSs)*  
   - Reads the CSV outputs, saved in batches, from **tffetcher**.
   - Loads *TFBS data* from a resource (e.g. RegulonDB `.txt` or `.csv`).  
   - For each TF identified in step 2, fetch the corresponding binding site(s).  
   - Saves a final CSV output, `tf2tfbs_mapping`, containing the following information:
  
      | Column              | Description |
      |---------------------|--------------------------------------------------------------------------------------------------|
      | **tf**              | The transcription factor name (lowercase)                                                        |
      | **tfbs**            | The TF binding site sequence (normalized to uppercase letters)                                   |
      | **gene**            | The DEG (or gene) associated with this TF mapping                                                |
      | **deg_source**      | A hyphen-delimited string indicating the source DEG datasets. For example, `houser_down-schmidt_down-durfee_down-bie_up` shows that multiple DEG datasets contributed to this mapping.                                                                                             |
      | **polarity**        | The regulatory polarity (e.g., "-" for repression, "+" for activation)                           |
      | **tfbs_source**     | A string indicating which TFBS resource(s) supplied the binding site (e.g., "regdb" or "ecocyc") |
      | **is_sigma_factor** | Boolean flag (e.g., "no") indicating whether the TF is a sigma factor                            |
      | **is_global_regulator** | Boolean flag (e.g., "no") indicating whether the TF is a global regulator                    |

   ***Warning:*** In ```tfbsfetcher.py```, the TFBS duplication check is stringent and does not account for cases where binding sites differ by ±1 nucleotide at either end.


#### Directory Layout
```text
deg2tfbs/
├── README.md
├── __init__.py
├── main.py                         # CLI entry point
├── configs/                        # User-defined configurations
│   └── example.yaml                # Customize to process different DEGs and retrieve different TFBSs
└── pipeline/
   ├── utils.py                    # Dictionary of paths to datasets (from dnadesign-data)
   ├── degfetcher/                 # Step 1: Isolate DEGs
   │   ├── __init__.py 
   │   ├── <dataset>_module.py     # Each omics dataset has its own respective module
   │   └── degbatch_<date>/        # Batch of DEGs retrieved from degfetcher in a given run
   │       ├── csvs                
   │       └── plots
   ├── tffetcher/                  # Step 2: Map DEGs to TFs
   │   ├── __init__.py    
   │   ├── tffetcher.py            # Coordinates TF retrieval for a given DEG batch
   │   ├── parsers/                
   │   │   ├── __init__.py
   │   │   ├── ecocyc_parser.py
   │   │   ├── regdb_parser.py     
   │   │   └── ...                 
   │   └── tfbatch_<date>/         
   │       └── deg2tf_mapping.csv  
   └── tfbsfetcher/                # Step 3: Map TFs to TFBSs
   │   ├── __init__.py    
   │   ├── tfbsfetcher.py          # Coordinates TFBS retrieval for a given set of TFs
   │   ├── parsers/                
   │   │   ├── __init__.py
   │   │   ├── ecocyc_tfbs_parser.py
   │   │   ├── regdb_tfbs_parser.py    
   │   │   └── ...                 
   │   └── tfbsbatch_<date>/  
```

## **Usage Example**

1. Clone the [**dnadesign-data**](https://github.com/e-south/dnadesign-data) repository to access a curated set of comparative omics datasets. Placing it as a sibling directory to **deg2tfbs** enables **degfetcher** to generate custom DEG tables from these sources. 

2. Update `configs/mycustomparams.yaml` with desired pointers to omics datasets, thresholds, and I/O paths.  
   ```yaml
   # deg2tfbs/configs/example.yaml

   pipeline:
   name: "default"

   stages:
      degfetcher:
         root_dir: "pipeline/degfetcher"
         batch_id: "degbatch_20250130"  # Where to save isolated DEGs
         modules:                       # Which source datasets to reference
            - ceroni
            - mori

      tffetcher:
         root_dir: "pipeline/tffetcher"
         batch_id: "tfbatch_20250130"
         input:
         deg_batch_id: "degbatch_20250130" # Points to a degfetcher batch
         deg_csv_subdir: "csvs"            
            deg_csv_includes:              # Optional: "Only fetch TFs from these lists"
               - "bie_downregulated_degs.csv"
               - "sanchez_vasquez_upregulated_degs.csv"
         sources:
         regulatory_networks:
            ecocyc:
               path: "ecocyc_28_reg_network" # Reference to utils.DATA_FILES
               enabled: true
               parser: "ecocyc_network_v28-5"
            myothersource: ...
         params:
         network_strategy: "union"     # How to merge results from different parsers
         include_master_regulators: true
         include_sigma_factors: true

      tfbsfetcher:
         root_dir: "pipeline/tfbsfetcher"
         batch_id: "tfbsbatch_20250130"
         input:
            tf_batch_id: "tfbatch_20250130" # Points to a tffetcher batch
         sources:
            binding_sites:
            ecocyc:
               path: "ecocyc_28_tfbs_smart_table"  # Reference to utils.DATA_FILES
               ecocyc_motifs: true
            regdb:
               path: "regulondb_13_tf_ri_set"
               regulondb_pssm: false
         params:  # extend as needed
            pass
      ```

3. After configuring your `mycustomparams.yaml`, run the pipeline as follows:
   ```bash
   cd deg2tfbs
   python main.py # Make sure that the bottom of this module references your config file.
   ```
   - **degfetcher** will process tables of genes, isolate DEGs, and then save them as tidy CSV files.
   
   - **tffetcher** will reference upstream CSV files to generate a ```deg2tf_mapping.csv``` file.
    
   - **tfbsfetcher** will reference a ```deg2tf_mapping.csv``` file to generate a ```tf2tfbs_mapping.csv``` file.


### **Extendability**

The **degfetcher**, **tffetcher**, and **tfbsfetcher** steps are designed to be extendable. You can add your own dataset ingestion modules for **degfetcher** or develop custom parser modules for **tffetcher** and **tfbsfetcher**.

To use your own custom input genes, create a CSV with required columns: 'gene', 'source', 'comparison'. Then place it in a "batch" subdirectory within **degfetcher**, and point to it in your custom config file. Alternatively, if you added a new gene table in a cloned **dnadesign-data** repository, then update the ```utils.py``` dictionary in ```pipeline``` to ensure datasets can be found via ```DATA_FILES[...]```.

When creating a custom parser, make sure it conforms to one of the following interfaces:

- For **tffetcher** parsers:  
  Implement: gene, (TF, polarity)
  ```python
  parse(...) -> Dict[str, Set[Tuple[str, str]]]
  ```
  This signature allows **tffetcher.py** to process new **regulatory network** data the same.

- For **tfbsfetcher** parsers:  
  Implement: TF, Set(TFBS_1, TFBS_2, TFBS_3, ...)
  ```python
  parse(...) -> Dict[str, Set[str]]
  ```
  This signature allows **tfbsfetcher.py** to process new **TFBS** data the same.

---

### Data Source
**deg2tfbs** is designed to reference data from [**dnadesign-data**](https://github.com/e-south/dnadesign-data). Update the config file to point to the correct data paths.
