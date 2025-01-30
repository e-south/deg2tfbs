
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
│       ├── radzikowski.py
│       ├── bie.py
│       ├── deter.py
│       ├── jovanovic.py
│       ├── rajacharya.py
│       ├── durfee.py
│       ├── gummesson.py
│       ├── houser.py
│       ├── lu.py
│       ├── sanchez_vasquez.py
│       ├── vazulka.py
│       ├── kim.py
│       └── zhang.py
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
