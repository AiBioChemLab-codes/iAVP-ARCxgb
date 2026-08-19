# iAVP-ARCXGB

iAVP-ARCXGB is a computational framework for antiviral peptide identification that combines eight handcrafted sequence-descriptor groups, ArcFace-derived embeddings, and an XGBoost classifier.

This repository contains the processed benchmark dataset, web-server deployment code, and the archived model artifacts associated with the manuscript experiments.

## Repository structure

```text
data/                  Processed training and test FASTA files
ml_model/              Model files used by the web-server deployment
utils/                 Handcrafted feature extraction code
reproducibility/       Archived manuscript model, features, optimization records, and evaluation script
docs/                  Data and usage notes
examples/              Example FASTA input
iAVP-ARCxgb.py         Streamlit application
```

## Installation

```bash
pip install -r requirements.txt
```

## Run the web application locally

```bash
streamlit run iAVP-ARCxgb.py
```

Input sequences must contain standard amino-acid characters only (`ACDEFGHIKLMNPQRSTVWY`) and have a length of 3-50 residues. Up to 200 sequences are processed per submission.

Example input:

```text
>peptide_1
IFWDCWAPEEPACQDFLGAMIH
>peptide_2
LSELDDRADALQAGASQFETSAAKLKRKYWWKN
```

## Dataset

The processed dataset contains 1,334 peptide sequences:

- Training set: 1,067 sequences (243 AVPs and 824 non-AVPs)
- Test set: 267 sequences (49 AVPs and 218 non-AVPs)

The exact sequences used in the study are provided in `data/train.fasta` and `data/test.fasta`. FASTA headers use `|1` for AVP and `|0` for non-AVP labels. Additional information is provided in `data/README.md` and `docs/DATA_PROVENANCE.md`.

## Manuscript reproduction

The model artifacts corresponding to the manuscript experiments are stored under `reproducibility/`.

To reproduce the archived manuscript evaluation from the released feature matrix:

```bash
python reproducibility/reproduce_reported_results.py
```

To regenerate the 2,125 handcrafted features directly from `data/test.fasta` before evaluation:

```bash
python reproducibility/reproduce_reported_results.py --from-fasta
```

The expected manuscript metrics are ACC = 0.873, MCC = 0.568, Precision = 0.660, F1 = 0.646, Sn = 0.633, Sp = 0.927, auROC = 0.837, auPRC = 0.558, and BA = 0.780 after rounding to three decimal places.

See `REPRODUCIBILITY.md` for the file map and reproduction workflow.

## Manuscript model and deployment model

The archived manuscript model is stored under `reproducibility/paper_model/`. The model used by the web application is stored separately under `ml_model/iAVP_ARCfaceXGB/`. These serialized model artifacts correspond to different saved model instances and should not be used interchangeably when reproducing the manuscript results.

## Web server

The online prediction server is available at:

https://iservers.aibiochem.net/iAVP-ARCxgb

A short usage guide is provided in `docs/WEB_SERVER_USAGE.md`.

## Feature extraction attribution

The handcrafted feature extraction implementation is based in part on descriptor routines and data used by iLearnPlus/iFeature. Please see `NOTICE.md` for attribution and the relevant references.

## Citation

If you use this repository, please cite the associated iAVP-ARCXGB manuscript and the software/data sources listed in the manuscript.
