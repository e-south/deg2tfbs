
### Overview
Pipeline for differential expression analysis (DEGs), transcription factor mapping, and TFBS identification.


The pipeline steps (DEG analysis → TF mapping → TFBS identification).


# Data
Requires the raw data from ```dnadesign-data```. Update the config file to point to the correct data paths.

```
deg2tfbs/
├── analysis/
│   ├── deg_analysis.py
│   ├── tf_association.py
│   ├── tfbs_identification.py
│   └── ...
├── pipeline.py
├── requirements.txt
├── README.md
└── ...
```