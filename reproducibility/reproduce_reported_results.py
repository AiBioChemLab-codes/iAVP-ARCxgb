from pathlib import Path
import argparse
import sys

import joblib
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)

ROOT = Path(__file__).resolve().parents[1]
REPRO_DIR = ROOT / "reproducibility"
MODEL_DIR = REPRO_DIR / "paper_model"
FEATURE_DIR = REPRO_DIR / "features"
RESULT_DIR = REPRO_DIR / "results"


class TabularEmbedNet(nn.Module):
    def __init__(self, input_dim, embed_dim=256):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.BatchNorm1d(256),
            nn.LeakyReLU(0.2),
            nn.Dropout(0.3),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.LeakyReLU(0.2),
            nn.Dropout(0.3),
            nn.Linear(128, embed_dim),
        )

    def forward(self, x):
        return self.net(x)


def parse_fasta(path):
    rows = []
    current_header = None
    current_seq = []

    def add_record(header, sequence):
        if header is None:
            return
        if "|" not in header:
            raise ValueError(f"Missing class label in FASTA header: {header}")
        seq_id, label = header.rsplit("|", 1)
        rows.append({"ID": seq_id, "Sequence": sequence, "Target": int(label)})

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                add_record(current_header, "".join(current_seq))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())
        add_record(current_header, "".join(current_seq))

    return pd.DataFrame(rows)


def load_handcrafted_features(from_fasta=False):
    if not from_fasta:
        df = pd.read_csv(FEATURE_DIR / "test_handcrafted_2125.csv.gz")
        return df.iloc[:, 0].astype(str).values, df.iloc[:, 1].astype(int).values, df.iloc[:, 2:].to_numpy(dtype=float)

    sys.path.insert(0, str(ROOT))
    from utils.getArtFeat import get8artfeat

    seq_df = parse_fasta(ROOT / "data" / "test.fasta")
    feature_input = seq_df[["ID", "Sequence"]].copy()
    feature_df = get8artfeat(feature_input)
    X = feature_df.iloc[:, 2:].to_numpy(dtype=float)
    return seq_df["ID"].astype(str).values, seq_df["Target"].astype(int).values, X


def load_arcface_model(device):
    state_dict = torch.load(MODEL_DIR / "arcface_encoder_state_dict.pt", map_location=device)
    model = TabularEmbedNet(input_dim=2125, embed_dim=256)
    model.load_state_dict(state_dict)
    model.to(device)
    model.eval()
    return model


def select_operating_threshold(y_true, y_prob):
    auprc = average_precision_score(y_true, y_prob)
    best_score = -np.inf
    best_threshold = None

    for threshold in np.unique(y_prob):
        y_pred = (y_prob >= threshold).astype(int)
        score = (
            accuracy_score(y_true, y_pred)
            + matthews_corrcoef(y_true, y_pred)
            + precision_score(y_true, y_pred, zero_division=0)
            + f1_score(y_true, y_pred, zero_division=0)
            + auprc
        ) / 5.0
        if score > best_score:
            best_score = score
            best_threshold = float(threshold)

    return best_threshold, best_score


def calculate_metrics(y_true, y_prob, threshold):
    y_pred = (y_prob >= threshold).astype(int)
    return {
        "Accuracy": accuracy_score(y_true, y_pred),
        "MCC": matthews_corrcoef(y_true, y_pred),
        "Sensitivity": recall_score(y_true, y_pred, pos_label=1),
        "Specificity": recall_score(y_true, y_pred, pos_label=0),
        "AUC_ROC": roc_auc_score(y_true, y_prob),
        "Precision": precision_score(y_true, y_pred, zero_division=0),
        "F1": f1_score(y_true, y_pred, zero_division=0),
        "AUC_PR": average_precision_score(y_true, y_prob),
        "BA": balanced_accuracy_score(y_true, y_pred),
    }


def main():
    parser = argparse.ArgumentParser(description="Reproduce the archived iAVP-ARCXGB manuscript evaluation.")
    parser.add_argument(
        "--from-fasta",
        action="store_true",
        help="Regenerate the 2,125 handcrafted features from data/test.fasta instead of loading the archived feature matrix.",
    )
    args = parser.parse_args()

    sample_ids, y_true, X = load_handcrafted_features(from_fasta=args.from_fasta)
    if X.shape[1] != 2125:
        raise ValueError(f"Expected 2125 handcrafted features, got {X.shape[1]}")

    device = torch.device("cpu")
    scaler = joblib.load(MODEL_DIR / "scaler.joblib")
    xgb_model = joblib.load(MODEL_DIR / "xgb_model.joblib")
    arc_model = load_arcface_model(device)

    X_scaled = scaler.transform(X)
    with torch.no_grad():
        embedding = arc_model(torch.tensor(X_scaled, dtype=torch.float32, device=device)).cpu().numpy()

    archived_embedding = pd.read_csv(FEATURE_DIR / "test_arcface_embedding_256.csv.gz")
    archived_ids = archived_embedding.iloc[:, 0].astype(str).values
    if not np.array_equal(sample_ids, archived_ids):
        raise ValueError("Sample order does not match the archived test embedding file.")
    emb_diff = np.max(np.abs(embedding - archived_embedding.iloc[:, 2:].to_numpy(dtype=float)))

    X_combined = np.hstack([X_scaled, embedding])
    y_prob = xgb_model.predict_proba(X_combined)[:, 1]
    threshold, opt_score = select_operating_threshold(y_true, y_prob)
    metrics = calculate_metrics(y_true, y_prob, threshold)

    params = pd.read_csv(MODEL_DIR / "best_params.csv").iloc[0]
    expected = pd.read_csv(RESULT_DIR / "reported_test_results.csv").iloc[0]

    print(f"Samples: {len(y_true)}")
    print(f"Handcrafted features: {X.shape[1]}")
    print(f"ArcFace embedding: {embedding.shape[1]}")
    print(f"Maximum embedding difference vs archived CSV: {emb_diff:.3e}")
    print(f"Selected operating threshold: {threshold:.10f}")
    print(f"Archived rounded threshold: {float(params['optimal_threshold']):.7f}")
    print(f"Composite optimization score: {opt_score:.10f}")
    print()
    print("Metric\tReproduced\tArchived")
    for key in ["Accuracy", "MCC", "Sensitivity", "Specificity", "AUC_ROC", "Precision", "F1", "AUC_PR"]:
        print(f"{key}\t{metrics[key]:.10f}\t{float(expected[key]):.10f}")
    print(f"BA\t{metrics['BA']:.10f}\t-")

    max_metric_diff = max(abs(metrics[key] - float(expected[key])) for key in expected.index)
    if max_metric_diff > 1e-8:
        raise RuntimeError(f"Reproduced metrics differ from the archived results (max difference {max_metric_diff:.3e}).")

    print("\nArchived manuscript metrics reproduced successfully.")


if __name__ == "__main__":
    main()
