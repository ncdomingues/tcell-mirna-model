# -*- coding: utf-8 -*-
"""
Machine-learning fate prediction (2026 addition).

Two questions, using two very differently-shaped datasets:

1. PYTHON MODEL (rich data): across the 400-run asynchronous ensembles already
   used in run_analysis.py, can a classifier predict a trajectory's eventual
   Th-phenotype from its own EARLY (partial) state, before it has converged?
   This uses genuine per-run trajectory data -- every run's full sweep-by-sweep
   history is available because we wrote the simulator ourselves.

2. "R MODEL" (GINsim/MaBoSS, 2019): the original thesis never published
   per-steady-state node vectors, only condition-level aggregate percentages
   (Table 4a/4b) -- so there is no per-sample FEATURE data to train on at all,
   only a handful of condition-level records. We reconstruct the closest
   defensible dataset from those published percentages (see
   reconstruct_r_model_dataset() -- every number is clearly derived, nothing
   invented) and run the same classification task on it, using "condition"
   (miRNA x rule-variant x source model) as the only available feature.

The comparison is therefore not just "which model predicts differentiation
better" but "what kind of downstream analysis does each implementation's
published output even support" -- a genuine, and fairly stark, difference
between a from-scratch reimplementation that keeps every run's raw trajectory
and a 2019 analysis that only ever reported aggregated percentages.
"""
import os
import json
import random
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

from rules_data import RULES, ALIASES
from logic_engine import LogicalModel
from run_analysis import build_condition_models, base_initial_state, MASTER_TFS, WATCHLIST

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "output")
FIG_DIR = os.path.join(HERE, "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

N_SWEEPS = 10
N_RUNS = 400
SEED = 42
CAPTURE_SWEEPS = [1, 2, 3, 5, 7, 10]  # 10 = fully converged (upper-bound sanity check)

MASTER_TF_NODES = list(MASTER_TFS.keys())
UPSTREAM_NODES = [n for n in WATCHLIST if n not in MASTER_TF_NODES]


def phenotype_label(final_state):
    active = [MASTER_TFS[tf] for tf in MASTER_TF_NODES if final_state.get(tf, 0) >= 1]
    if len(active) == 0:
        return "None"
    if len(active) == 1:
        return active[0]
    return "Mixed"


# ---------------------------------------------------------------- 1. Python-model dataset
def collect_python_dataset():
    """Re-simulate all 4 conditions, this time keeping each run's own sweep-by-
    sweep history (async_run already returns it; run_analysis.ensemble_run
    discards it after averaging). Same rules, same initial conditions, same
    seed scheme as run_analysis.py / main_results_with_ci.py, so results are
    the same steady-state distributions already reported in Results 3.3.4."""
    conditions = build_condition_models()
    rows = []
    for label, (model, overrides) in conditions.items():
        init = base_initial_state(model)
        init.update(overrides)
        rng = random.Random(SEED)
        n_steps = N_SWEEPS * len(model.rule_nodes)
        for run in range(N_RUNS):
            final_state, history = model.async_run(init, n_steps, rng)
            pheno = phenotype_label(final_state)
            row = {"condition": label, "run": run, "phenotype": pheno}
            for sweep in CAPTURE_SWEEPS:
                snap = history[min(sweep, len(history) - 1)]
                for node in WATCHLIST:
                    row[f"sw{sweep}__{node}"] = int(snap.get(node, 0) >= 1)
            rows.append(row)
    return pd.DataFrame(rows)


MIN_CLASS_COUNT = 10


def collapse_rare_classes(y, min_count=MIN_CLASS_COUNT):
    """Classes with too few members (e.g. Th9/SPI1-PU1, which this ensemble
    only ever hit once in 1600 runs) can't be stratified or meaningfully
    learned. Collapse everything below min_count into one 'Rare' bucket; if
    that bucket itself ends up with under 2 members (stratified splitting
    needs >=2 of every class), drop those rows outright rather than crash --
    the exclusion is reported, not silently swallowed."""
    counts = y.value_counts()
    rare = set(counts[counts < min_count].index)
    if not rare:
        return y, pd.Series([True] * len(y), index=y.index)
    collapsed = y.apply(lambda v: f"Rare (<{min_count} samples)" if v in rare else v)
    rare_bucket_n = (collapsed == f"Rare (<{min_count} samples)").sum()
    if rare_bucket_n >= 2:
        return collapsed, pd.Series([True] * len(y), index=y.index)
    keep = ~y.isin(rare)
    print(f"  dropping {(~keep).sum()} run(s) from ultra-rare classes "
          f"{sorted(set(y[~keep]))} (fewer than 2 total, cannot be split/stratified)")
    return y[keep], keep


def run_python_ml(df):
    """For each capture sweep, train/evaluate two feature sets:
    'full' (all WATCHLIST nodes) and 'upstream' (excludes the master TFs
    themselves) -- the latter asks whether fate is legible before any master
    regulator has committed."""
    results = []
    y_all, keep = collapse_rare_classes(df["phenotype"])
    for sweep in CAPTURE_SWEEPS:
        for feat_set, nodes in [("full", WATCHLIST), ("upstream_only", UPSTREAM_NODES)]:
            cols = [f"sw{sweep}__{n}" for n in nodes]
            X = df.loc[keep, cols]
            y = y_all
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=0.25, random_state=0, stratify=y
            )
            clf = RandomForestClassifier(n_estimators=300, random_state=0, min_samples_leaf=2)
            clf.fit(X_train, y_train)
            pred = clf.predict(X_test)
            acc = accuracy_score(y_test, pred)

            baseline = DummyClassifier(strategy="most_frequent", random_state=0)
            baseline.fit(X_train, y_train)
            base_acc = accuracy_score(y_test, baseline.predict(X_test))

            results.append({
                "sweep": sweep, "sweep_frac": round(sweep / N_SWEEPS, 2),
                "feature_set": feat_set, "n_features": len(cols),
                "accuracy": round(acc, 4), "majority_class_baseline": round(base_acc, 4),
                "n_train": len(X_train), "n_test": len(X_test),
            })
            if sweep == max(CAPTURE_SWEEPS) - (max(CAPTURE_SWEEPS) % N_SWEEPS) and feat_set == "upstream_only":
                pass  # placeholder, importances handled separately below
    return pd.DataFrame(results)


def feature_importance_at(df, sweep, nodes, feat_set_name):
    cols = [f"sw{sweep}__{n}" for n in nodes]
    y, keep = collapse_rare_classes(df["phenotype"])
    X = df.loc[keep, cols]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=0, stratify=y)
    clf = RandomForestClassifier(n_estimators=300, random_state=0, min_samples_leaf=2)
    clf.fit(X_train, y_train)
    pred = clf.predict(X_test)
    importances = pd.Series(clf.feature_importances_, index=nodes).sort_values(ascending=False)
    cm = confusion_matrix(y_test, pred, labels=sorted(y.unique()))
    report = classification_report(y_test, pred, labels=sorted(y.unique()), zero_division=0, output_dict=True)
    return importances, cm, sorted(y.unique()), report


# ---------------------------------------------------------------- 2. Reconstructed R-model dataset
def reconstruct_r_model_dataset():
    """Back out approximate integer steady-state counts from the published
    Table 4a/4b percentages (Naldi 2010 and Abou-Jaoude 2015 models, original
    2019 thesis). Only cells the original table actually reported are used --
    every blank cell in the source (miR-34c-5p and 'both' conditions have no
    RORC/FOXP3/Mixed/SSs(n) entries for the Naldi model; only TH1/TH2 and
    SSs(n) are reported for those two conditions in the Abou-Jaoude model) is
    left out rather than filled in. counts = round(pct/100 * SSs(n)) per
    master TF, independently (the original table reports per-TF activation
    frequency across a condition's enumerated steady states, not a single
    mutually-exclusive label -- same convention as MASTER_TFS/phenotype_table
    in run_analysis.py) -- so this reconstruction has the same shape as the
    Python-side ensemble output, just tiny and condition-only (no per-sample
    node-state features exist in the source to reconstruct)."""
    records = []

    naldi_full = [  # condition, variant, SSs(n), {TF: pct}
        ("no miRNA", "A", 12, {"TBX21": 42, "GATA3": 17, "RORC": 42, "FOXP3": 33}),
        ("no miRNA", "B", 9, {"TBX21": 42, "GATA3": 17, "RORC": 42, "FOXP3": 33}),
        ("miR-155-5p only", "A", 2, {"TBX21": 44, "GATA3": 22, "RORC": 78, "FOXP3": 44}),
        ("miR-155-5p only", "B", 2, {"TBX21": 44, "GATA3": 11, "RORC": 56, "FOXP3": 44}),
    ]
    for cond, variant, n_ss, pcts in naldi_full:
        for tf, pct in pcts.items():
            n_active = round(pct / 100 * n_ss)
            for i in range(n_ss):
                records.append({
                    "source_model": "Naldi_2010", "variant": variant, "condition": cond,
                    "master_tf": tf, "active": 1 if i < n_active else 0,
                })

    abou_ss = {"no miRNA": 12, "miR-155-5p only": 4, "miR-34c-5p only": 4, "both": 4}
    abou_full_pcts = {"no miRNA": {"TBX21": 50, "GATA3": 33, "RORC-STAT3+": 8,
                                     "SPI1": 33, "FOXP3": 33, "STAT3_IL22": 8}}
    abou_partial_pcts = {  # only TH1/TH2 reported for these three conditions
        "miR-155-5p only": {"TBX21": 50, "GATA3": 50},
        "miR-34c-5p only": {"TBX21": 50, "GATA3": 50},
        "both": {"TBX21": 50, "GATA3": 50},
    }
    for cond, n_ss in abou_ss.items():
        pcts = abou_full_pcts.get(cond) or abou_partial_pcts[cond]
        for tf, pct in pcts.items():
            n_active = round(pct / 100 * n_ss)
            for i in range(n_ss):
                records.append({
                    "source_model": "Abou-Jaoude_2015", "variant": "-", "condition": cond,
                    "master_tf": tf, "active": 1 if i < n_active else 0,
                })

    return pd.DataFrame(records)


def run_r_model_ml(df_r):
    """Only feature available is the condition label itself (source model x
    variant x miRNA condition x which master TF is being asked about) -- there
    is no per-sample node-state data in the source to use as a real feature.
    This is deliberately a near-degenerate classification problem; the point
    is to make that limitation explicit and quantified, not to hide it."""
    df = df_r.copy()
    X = pd.get_dummies(df[["source_model", "variant", "condition", "master_tf"]])
    y = df["active"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=0, stratify=y
    )
    clf = RandomForestClassifier(n_estimators=300, random_state=0, min_samples_leaf=2)
    clf.fit(X_train, y_train)
    acc = accuracy_score(y_test, clf.predict(X_test))

    baseline = DummyClassifier(strategy="most_frequent", random_state=0)
    baseline.fit(X_train, y_train)
    base_acc = accuracy_score(y_test, baseline.predict(X_test))

    return {
        "n_total_records": len(df), "n_train": len(X_train), "n_test": len(X_test),
        "n_features": X.shape[1], "accuracy": round(acc, 4),
        "majority_class_baseline": round(base_acc, 4),
    }


# ---------------------------------------------------------------- main
def main():
    print("Simulating Python-model dataset (this re-runs the same 4x400 ensembles)...")
    df_py = collect_python_dataset()
    df_py.to_csv(os.path.join(OUT_DIR, "ml_python_dataset.csv"), index=False)
    print(f"  {len(df_py)} runs collected.")

    print("Training classifiers across capture sweeps...")
    py_results = run_python_ml(df_py)
    py_results.to_csv(os.path.join(OUT_DIR, "ml_python_accuracy_by_sweep.csv"), index=False)
    print(py_results.to_string(index=False))

    importances, cm, labels, report = feature_importance_at(df_py, sweep=3, nodes=UPSTREAM_NODES, feat_set_name="upstream_only")
    importances.to_csv(os.path.join(OUT_DIR, "ml_feature_importance_sweep3_upstream.csv"))
    with open(os.path.join(OUT_DIR, "ml_confusion_matrix_sweep3_upstream.json"), "w", encoding="utf-8") as f:
        json.dump({"labels": labels, "matrix": cm.tolist(), "report": report}, f, indent=2)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for feat_set, sub in py_results.groupby("feature_set"):
        ax.plot(sub["sweep_frac"], sub["accuracy"], marker="o", label=feat_set)
    ax.axhline(py_results["majority_class_baseline"].iloc[0], color="grey", linestyle="--",
               linewidth=1, label="majority-class baseline")
    ax.set_xlabel("fraction of trajectory elapsed (sweep / 10)")
    ax.set_ylabel("held-out accuracy (7-class phenotype)")
    ax.set_title("How early can fate be predicted? RandomForest accuracy vs. trajectory progress\n"
                  "(Python model, 1600 runs across 4 miRNA conditions)", fontsize=11, loc="left")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "ml_accuracy_vs_sweep.png"), dpi=170)
    plt.close(fig)

    print("\nReconstructing R-model (GINsim/MaBoSS 2019) dataset from published percentages...")
    df_r = reconstruct_r_model_dataset()
    df_r.to_csv(os.path.join(OUT_DIR, "ml_r_model_reconstructed_dataset.csv"), index=False)
    print(f"  {len(df_r)} reconstructed pseudo-records "
          f"(condition-level only -- no per-sample node states exist in the source).")

    r_results = run_r_model_ml(df_r)
    with open(os.path.join(OUT_DIR, "ml_r_model_results.json"), "w", encoding="utf-8") as f:
        json.dump(r_results, f, indent=2)
    print(r_results)

    summary = {
        "python_model": {
            "n_samples": len(df_py),
            "n_features_available_per_sample": len(WATCHLIST),
            "best_accuracy": float(py_results["accuracy"].max()),
            "best_config": py_results.loc[py_results["accuracy"].idxmax()][["sweep", "feature_set"]].to_dict(),
            "majority_class_baseline": float(py_results["majority_class_baseline"].iloc[0]),
        },
        "r_model_reconstructed": r_results,
    }
    with open(os.path.join(OUT_DIR, "ml_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print("\nSaved all ML outputs to", OUT_DIR, "and", FIG_DIR)


if __name__ == "__main__":
    main()
