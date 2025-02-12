deg2tfbs/
├── LICENSE
├── README.md
├── environment.yml
├── images
├── pyproject.toml
└── src/
    └── deg2tfbs/
        ├── __init__.py
        ├── analysis/
        │   ├── __init__.py
        │   └── tfbsbatch_plots/
        │       ├── __init__.py
        │       ├── config_utils.py
        │       ├── data.py
        │       ├── roster.py
        │       ├── compare.py
        │       ├── umap_clustering.py
        │       └── main.py
        ├── configs/
        │   └── example.yaml
        ├── main.py
        └── pipeline/
            └── tfbsfetcher/
                ├── __init__.py
                ├── tfbsfetcher.py
                └── … (other modules and tfbatch_* directories)
