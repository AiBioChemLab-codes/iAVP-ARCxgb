# Archived Hyperopt trials

`hyperopt_trials.csv.gz` contains the 3,000 recorded optimization trials associated with the manuscript experiment. The original individual trial CSV files were combined into one gzip-compressed table without changing their recorded values.

The table includes the candidate ArcFace/XGBoost settings, evaluation metrics, composite optimization score, and operating threshold recorded for each trial.

The archived best record is trial 1020. Its reported settings include:

- ArcFace embedding dimension: 256
- ArcFace epochs: 1
- XGBoost learning rate: 0.10992458218499501
- max_depth: 6
- min_child_weight: 5
- subsample: 0.7665346220918834
- colsample_bytree: 0.933033839018152
- reg_alpha: 3.277393264080812e-07
- reg_lambda: 7.788300007047023e-05
- scale_pos_weight: 3.390946488103101
- n_estimators: 200
- random_state: 42
- archived rounded threshold: 0.4234547

The recorded composite score is the arithmetic mean of ACC, MCC, Precision, F1, and auPRC.
