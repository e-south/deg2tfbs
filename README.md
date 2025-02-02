## **deg2tfbs**

A pipeline for processing RNA-seq and proteomic datasets, focusing on *E. coli* but extendable to other organisms. This tool identifies differentially expressed genes (DEGs) between conditions,maps them to upstream transcription factors using *RegulonDB* or *EcoCyc*, and retrieves their corresponding transcription factor binding sites (TFBSs).  

- [figure]

## **Pipeline Steps**

1. **degfetcher** *(Step 1: Isolate DEGs)*  
   - Loads comparative omics datasets (e.g., `ceroni.py`, `mori.py`) from the ```dnadesign-dna``` repository.
   - Produces tidy CSV outputs, e.g. `ceroni_upregulated_degs.csv`, containing columns:  
     - **gene** – the gene name or locus  
     - **source** – the dataset (e.g., `"ceroni"`)  
     - **comparison** – the target vs. reference condition (e.g., `"pLys-M1_vs_pLys"`)  
     - **thresholds** – optional numeric or descriptive threshold used  

2. **tf_fetcher** *(Step 2: Map DEGs to TFs)*  
   - Reads the CSV outputs, saved in batches, from **degfetcher**.
   - Loads a *regulatory interaction* dataset(s) (e.g., from **EcoCyc** or **RegulonDB**).  
   - **Identifies** which transcription factors (TFs) regulate these DEGs (the “regulatees”).  
   - Outputs a CSV of all `(gene, regulator, source, comparison, regulatorsource)` rows.  
     - **example**: `dnaK, rpoH, ceroni, pLys-M1_vs_pLys, ecocyc_v28.5`  

3. **tfbs_fetcher** *(Step 3: Identify TF binding sites)*  
   - Loads TFBS data from a resource (e.g. RegulonDB `.txt` or `.csv`).  
   - For each TF identified in step 2, fetch the corresponding binding site(s).  
   - Saves a final CSV that can be used downstream (e.g., in **dnadesign**).  

## Project Layout

```text
deg2tfbs/
├── README.md
├── __init__.py
├── main.py                         # CLI entry point
├── configs/                        # User-defined variables live here
│   └── example.yaml
│   
└── pipeline/                       
    ├── degfetcher/                 # Step 1: Run modules to ingest and fetch DEGs from datasets
    │   ├── __init__.py 
    │   ├── <dataset>_module.py     # Each dataset has its own respective module
    │   ├── <dataset>_module.py
    │   ├── ...
    │   └── degbatch_<date>/        # Batch of DEGs collected from degfetcher datasets
    │       └── csvs                # In tidy format (req. cols: 'gene', 'source', and 'comparison')
    │       └── plots               # Some modules generate plots to showcase DEGs from non-DEGs
    │
    ├── tffetcher/                  # Step 2: 
    │   ├── __init__.py    
    │   ├── tffetcher.py
    │   ├── parsers/
    │   │   ├── __init__.py         # get_regulatory_parser factory
    │   │   ├── ecocyc_parser.py    # EcoCycParser class
    │   │   ├── regdb_parser.py     # RegulonDBParser class
    │   │   └── ...                 # Extendable to reference your own DEG->TF mapping file 
    │   │
    │   └── tfbatch_<date>/         # Batch of DEG-TF pairings collected from tffetcher datasets
    │       └── csvs                # In tidy format (req. cols: XXX)
    │       └── plots  
    │
    ├── tfbsfetcher/                # Step 3: 
    │
    └── utils.py                    # Dictionary of paths to datasets in the dnadesign-data repository
```


## **Usage**

1. **Load** raw data from **dnadesign-data** or create your own module in ```degfetcher``` to prepare custom lists of DEGs.
2. **Update** `configs/example.yaml` with desired datasets, thresholds, and output paths.  
3. **Run** from the command line:  
   ```bash
   cd deg2tfbs
   python main.py
   ```

   this will generate csv and perhaps a plot.
  - **[plot example]**

4. The **degfetcher** modules produce CSVs in `outputs/.../csv/`. 

#### Example Output from `degfetcher`
  ```csv
  gene,source,thresholds,comparison
  dnaK,ceroni,2.5,pLys-M1_versus_pLys
  yaaU,ceroni,2.5,pLys-M1_versus_pLys
  ...
  ```

5. Then, a **tffetcher** module can read those CSVs to compute the `(gene, TF)` pairs.

Extensibility: If you add more networks (e.g., a custom “MyNetworkParser”), you only need to put an additional entry in your YAML + a parser class, and tffetcher will pick it up.
Confidence filters: If you want to filter out low-confidence RegulonDB interactions, that can be handled inside the parser class.

```csv
gene,tf,source,comparison,tf_source
dnaK,rpoH,ceroni,pLys-M1_versus_pLys,ecocyc_v27.5
```

## Usage

1. **Install** dependencies **(WIP)**:  
   ```bash
   pip install -r requirements.txt
   ```  
2. If you're interesting in loading custom datasets, **configure** your `configs/example.yaml` to set input data paths, thresholds, and network references. you can also add data to dnadesign-data and thenObtain the raw data (e.g., from **dnadesign-data**) and place it so that the ```utils.py``` dictionary can locate each dataset via ```DATA_FILES[...]```.
      - Configure ```configs/example.yaml```:
      - Adjust ```output.root_dir``` to point to ```outputs```.
      - Change ```output.batch_identifier``` to e.g. ```batch_2023_08_02```.
      - Modify each dataset’s ```data``` and ```thresholds``` as desired.
3. **Run** the pipeline:  
   ```bash
   cd deg2tfbs
   python main.py
   ```  
   - The **degfetcher** modules produce CSVs like `ceroni_upregulated_degs.csv`. 
      - generates tidy DataFrames of up- and down-regulated genes,  
      - <dataset-name>_<typ-of-regulation>_degs.csv

      e.g., the head of ```ceroni_upregulated_degs.csv```
      ```csv
      gene	source	thresholds	comparison
      dnaK	ceroni	2.5	pLys-M1_versus_pLys
      yaaU	ceroni	2.5	pLys-M1_versus_pLys
      yadK	ceroni	2.5	pLys-M1_versus_pLys
      mmuP	ceroni	2.5	pLys-M1_versus_pLys
      insE1	ceroni	2.5	pLys-M1_versus_pLys
      mhpE	ceroni	2.5	pLys-M1_versus_pLys
      ```

      each batch contains .csv files indicated as either, where some source have both, some have just one:
      <dataset-name>_upregulated_degs.csv
      <dataset-name>_downregulated_degs.csv


   - A **tf_fetcher** module merges these CSVs with a regulatory network, producing `tf_associations.csv`.  

    The pipeline steps (DEG analysis → TF mapping → TFBS identification).


### Data Source
Requires source data from [dnadesign-data](https://github.com/e-south/dnadesign-data). Update the config file to point to the correct data paths.



### Some list of all the available conditions
Heterologous Burden
DEGs from RNA-seq
Heterologous expression versus empty vector
'pLys-M1’ versus 'pLys’ (Ceroni et al., 2018)
'pPD864-LacZ’ versus 'pPD864’ (Ceroni et al., 2018)
'pSB1C3-H3’ versus 'pSB1C3’ (Ceroni et al., 2018)
DiMAs from RNA-seq
Heterologous expression versus empty vector
M9-glucose with and without MBP expression (Lamoureux et a., 2021)
M9-glucose with and without BRCA expression (Lamoureux et a., 2021)

Membrane Integrity
DEGs from RNA-seq
Protein expression with and without periplasmic secretion tag (Emani et al., 2023)
DiMAs from RNA-seq
Exposure to cell wall-targeting antibiotic
M9-glucose with and without Ceftriaxone (4.1mg/mL) (Lamoureux et a., 2021)
Deletion of regulators known to be involved in membrane response
LB MG1655 versus MG1655ΔcpxR Batch Culture (Lamoureux et a., 2021)
LB MG1655 versus MG1655ΔbaeR Batch Culture (Lamoureux et a., 2021)
LB MG1655 versus MG1655ΔrcsB Batch Culture (Lamoureux et a., 2021)

Starvation (Stringent Response)
DEGs from RNA-seq
Glucose uptake rate low and high (Fragoso-Jimenez et al., 2022) 
Grown in minimal media versus two-weeks starving (Houser et al., 2015)
Gene expression with and without heightened ppGpp (Sanchez-Vasquez et al., 2019)
DEGs from Comparative Proteomics
High synthesis priority in slow- and fast-growing cells (Peebo et al., 2015)
Persisters and non-persisters (Radzikowski et al., 2016)
Chemostat_µ=0.20 versus Chemostat_µ=0.5 (Schmidt et al., 2016)
Carbon limited MOPS versus Carbon Rich MOPs (Wu et al., 2023)
Overexpression or no expression of RelA (Zhu et al., 2023)
LB log-phase versus LB stationary (Mori et al., 2021)
DiMAs from RNA-seq
Deletion of regulators known to be involved in stringent response
Fed batch BW25113 versus BW25113ΔrpoS (6h) (Lamoureux et a., 2021)


