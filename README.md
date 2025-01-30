## **deg2tfbs**

A pipeline for processing RNA-seq and proteomic datasets, focusing on *E. coli* but extendable to other organisms. This tool identifies differentially expressed genes (DEGs) between conditions,maps them to upstream transcription factors using *RegulonDB* or *EcoCyc*, and retrieves their corresponding transcription factor binding sites (TFBSs).  

- [figure]

## **Pipeline Steps**

1. **degfetcher** *(Step 1: Identify DEGs)*  
   - Loads comparative omics datasets (e.g., `ceroni.py`, `mori.py`) from the ```dnadesign-dna``` repository.
   - Produces tidy CSV outputs, e.g. `ceroni_upregulated_degs.csv`, containing columns:  
     - **gene** – the gene name or locus  
     - **source** – the dataset (e.g., `"ceroni"`)  
     - **comparison** – the target vs. reference condition (e.g., `"pLys-M1_vs_pLys"`)  
     - **thresholds** – optional numeric or descriptive threshold used  

2. **tffetcher** *(Step 2: Map DEGs to TFs)*  
   - Reads the CSV outputs from **degfetcher**—the user may specify one or more.  
   - Loads a *regulatory interaction* dataset (e.g., from **EcoCyc** or **RegulonDB**).  
   - **Identifies** which transcription factors (TFs) regulate these DEGs (the “regulatees”).  
   - Outputs a CSV of all `(gene, TF, source, comparison, TF_source)` rows.  
     - **example**: `dnaK, rpoH, ceroni, pLys-M1_vs_pLys, ecocyc_v27.5`  

3. **tfbs_fetcher** *(Step 3: Identify TF binding sites)*  
   - Loads TFBS data from a resource (e.g. RegulonDB `.txt` or `.csv`).  
   - For each TF identified in step 2, fetch the corresponding binding site(s).  
   - Saves a final CSV that can be used downstream (e.g., in **dnadesign**).  

## Project Layout

```text
deg2tfbs/
├── __init__.py
├── main.py                         # Top-level API
├── configs/                        # User-defined variables live here
│   └── example.yaml
├── pipeline/                       
│   ├── degfetcher/                 # Contains modules to ingest and fetch DEGs from omics datasets
│   │   ├── <dataset>_module.py     # Each dataset has its own respective module
│   │   ├── <dataset>_module.py
│   │   ├── ...
│   │   └── degbatch_<date>/        # Batch of DEGs collected from a given degfetcher run
│   │       └── csvs                # In tidy format (req. cols: 'gene', 'source', and 'comparison')
│   │       └── plots               # Some modules generate plots to showcase DEGs from non-DEGs
│   ├── tf_fetcher/
│   │   ├── utils.py      
│   │   └── tf_fetcher.py           # Applies network of choice to identify regulators
│   └── tfbs_fetcher/

├── outputs/
│   └── ...
└── README.md
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
