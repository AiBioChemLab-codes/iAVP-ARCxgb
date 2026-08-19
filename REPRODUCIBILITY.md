# Reproducibility

## Overview

The `reproducibility/` directory contains the processed feature matrices, ArcFace embeddings, saved manuscript-model artifacts, recorded Hyperopt trials, cross-validation summaries, and the evaluation script associated with the reported iAVP-ARCXGB results.

## File map

### Features

- `features/train_handcrafted_2125.csv.gz`: 2,125-dimensional handcrafted features for the training set.
- `features/test_handcrafted_2125.csv.gz`: 2,125-dimensional handcrafted features for the test set.
- `features/train_arcface_embedding_256.csv.gz`: 256-dimensional ArcFace embeddings for the training set.
- `features/test_arcface_embedding_256.csv.gz`: 256-dimensional ArcFace embeddings for the test set.

### Manuscript model

- `paper_model/arcface_encoder_state_dict.pt`: ArcFace encoder state dictionary used to generate the archived manuscript embeddings.
- `paper_model/scaler.joblib`: scaler used for the manuscript feature matrix.
- `paper_model/xgb_model.joblib`: XGBoost classifier associated with the reported manuscript result.
- `paper_model/best_params.csv`: archived selected parameter record, including the operating threshold rounded in the original record.

### Results

- `results/reported_test_results.csv`: archived metric values corresponding to the manuscript result.
- `results/cv_folds.csv`: five-fold cross-validation records for the selected configuration.
- `results/cv_mean.csv`: cross-validation summary.

### Optimization

- `optimization/hyperopt_trials.csv.gz`: merged table of 3,000 archived Hyperopt trial records.
- `optimization/README.md`: column and parameter notes for the optimization archive.

## Reproducing the reported evaluation

From the repository root, run:

```bash
python reproducibility/reproduce_reported_results.py
```

The script performs the following steps:

1. Loads the 2,125-dimensional handcrafted feature matrix.
2. Applies the archived feature scaler.
3. Generates the 256-dimensional ArcFace embeddings with the manuscript encoder.
4. Concatenates the standardized handcrafted features and ArcFace embeddings into a 2,381-dimensional representation.
5. Generates prediction probabilities with the archived XGBoost classifier.
6. Reconstructs the recorded operating-threshold procedure.
7. Calculates the evaluation metrics and checks them against `results/reported_test_results.csv`.

To regenerate the handcrafted features directly from the FASTA sequences, run:

```bash
python reproducibility/reproduce_reported_results.py --from-fasta
```

The regenerated feature matrix is passed through the same saved scaler, ArcFace encoder, and XGBoost classifier.

## Expected output

The reproduced values should agree with the archived results before rounding:

| Metric | Value |
|---|---:|
| ACC | 0.8726591760 |
| MCC | 0.5684321877 |
| Sensitivity | 0.6326530612 |
| Specificity | 0.9266055046 |
| Precision | 0.6595744681 |
| F1 | 0.6458333333 |
| auROC | 0.8373900019 |
| auPRC | 0.5575281842 |
| BA | 0.7796292829 |

The archived parameter CSV stores the operating threshold rounded to `0.4234547`. The reproduction script recalculates the threshold at full floating-point precision from the saved prediction scores so that the archived metric values are reproduced exactly.

## Web deployment

The web-server model files under `ml_model/iAVP_ARCfaceXGB/` are maintained separately from the manuscript model files in this directory. Use `reproducibility/paper_model/` for reproduction of the manuscript evaluation.

## File integrity

SHA-256 checksums for the released reproducibility artifacts are listed in `reproducibility/SHA256SUMS.txt`.
