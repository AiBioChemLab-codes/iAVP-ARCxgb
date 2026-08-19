# Processed dataset

This directory contains the exact processed peptide sequences used in the iAVP-ARCXGB study.

## Files

- `train.fasta`: 1,067 sequences, including 243 AVPs and 824 non-AVPs.
- `test.fasta`: 267 sequences, including 49 AVPs and 218 non-AVPs.

## FASTA header format

Each header ends with the binary class label:

```text
>sequence_identifier|1
PEPTIDESEQUENCE
```

- `1`: antiviral peptide (AVP)
- `0`: non-antiviral peptide

## Dataset construction

The study integrated peptide records from AVPdb, HIPdb, APD3, CAMPR3, DBAASP, DRAMP, LAMP2, YADAMP, and published AVP datasets cited in the manuscript. The assembled sequences were standardized, duplicate records were removed, sequences longer than 50 amino acids were excluded, and CD-HIT was used at a 30% sequence-similarity threshold before the final 1,334-sequence benchmark dataset was formed.

The exact processed sequences used for model development and evaluation are supplied here so that downstream feature generation and model evaluation can use the same sequence set.
