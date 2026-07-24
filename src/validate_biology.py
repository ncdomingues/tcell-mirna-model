# -*- coding: utf-8 -*-
"""
Biological validation, independent of the miRNA question entirely.

The base differentiation model this rebuild extends (Naldi et al. 2010) was
originally validated by showing it reproduces the textbook mapping from
Th-polarising cytokine exposure to T-helper master-regulator commitment:

    IL12            -> TBX21 (Th1)
    IL4             -> GATA3 (Th2)
    TGFB + IL6       -> RORC  (Th17)
    TGFB (alone)     -> FOXP3 (iTreg)

This script runs the same test against the digitized rule table (with the
miRNAs entirely absent, so this only tests the transcription + parser +
simulator, not the miRNA hypothesis), by holding each cytokine node fixed as
an exogenous input rather than letting the cell produce it autocrinely, and
checking that the expected master regulator responds.

It also runs the thesis's own positive control: miR-155-5p is the
best-characterised immune miRNA, so its simulated effect direction can be
checked against established literature (Th1-promoting, Th2-repressing)
before trusting anything the model says about the much less characterised
miR-34c-5p.
"""
import os
import sys
import json
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rules_data import RULES, ALIASES
from logic_engine import LogicalModel

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "output")
os.makedirs(OUT_DIR, exist_ok=True)

N_SWEEPS = 10
N_RUNS = 300
SEED = 11

MASTER_TFS = {"TBX21": "Th1", "GATA3": "Th2", "RORC": "Th17", "FOXP3": "iTreg"}
LIGANDS = ["IL12", "IL4", "TGFB", "IL6"]
BASE_INPUTS = ("APC-Antigen", "CD3", "CD45RA", "CGC", "CREBBP")


# ---------------------------------------------------------- Part A: Th-polarisation
def build_polarisation_model():
    """No-miRNA network with the 4 polarising cytokines promoted to fixed
    exogenous inputs (rather than autocrine outputs), so their level can be
    set directly to mimic an in-vitro Th-polarising cocktail."""
    excluded = {"miR-34c-5p", "miR-155-5p"} | set(LIGANDS)
    rules = [r for r in RULES if r[0] not in excluded]
    return LogicalModel(rules, ALIASES)


def run_polarisation_condition(model, forced_ligands):
    init = {n: 0 for n in model.all_nodes}
    for inp in BASE_INPUTS:
        init[inp] = 1
    for lig in LIGANDS:
        init[lig] = 1 if lig in forced_ligands else 0
    finals, traj = model.ensemble_run(init, n_sweeps=N_SWEEPS, n_runs=N_RUNS, seed=SEED)
    n = len(finals)
    return {tf: round(100.0 * sum(1 for f in finals if f.get(tf, 0) >= 1) / n, 1) for tf in MASTER_TFS}


def part_a_th_polarisation():
    model = build_polarisation_model()
    conditions = {
        "baseline (TCR+IL2 only)": frozenset(),
        "+IL12 (expect Th1)": frozenset({"IL12"}),
        "+IL4 (expect Th2)": frozenset({"IL4"}),
        "+TGFB+IL6 (expect Th17)": frozenset({"TGFB", "IL6"}),
        "+TGFB only (expect iTreg)": frozenset({"TGFB"}),
    }
    rows = {}
    for label, ligs in conditions.items():
        rows[label] = run_polarisation_condition(model, ligs)
    df = pd.DataFrame(rows).T.rename(columns=MASTER_TFS)
    return df


def evaluate_part_a(df):
    """For each polarising condition, check the expected master TF increased
    versus baseline and is the top-scoring TF in that row."""
    baseline = df.loc["baseline (TCR+IL2 only)"]
    checks = [
        ("+IL12 (expect Th1)", "Th1"),
        ("+IL4 (expect Th2)", "Th2"),
        ("+TGFB+IL6 (expect Th17)", "Th17"),
        ("+TGFB only (expect iTreg)", "iTreg"),
    ]
    results = []
    for row_label, expected_tf in checks:
        row = df.loc[row_label]
        increased = row[expected_tf] > baseline[expected_tf]
        is_top = row[expected_tf] == row.max()
        results.append({
            "condition": row_label,
            "expected": expected_tf,
            "baseline_pct": baseline[expected_tf],
            "condition_pct": row[expected_tf],
            "increased_vs_baseline": increased,
            "is_top_response_in_row": is_top,
            "pass": increased and is_top,
        })
    return pd.DataFrame(results)


# ---------------------------------------------------------- Part B: miR-155-5p positive control
def part_b_mir155_positive_control(main_results_csv):
    df = pd.read_csv(main_results_csv, index_col=0)
    no_mir = df["no miRNA"]
    mir155 = df["miR-155-5p only"]
    checks = [
        ("TBX21", "Th1 up (miR-155-5p is Th1-promoting; O'Connell et al. 2010)", "up"),
        ("GATA3", "Th2 down (miR-155-5p represses the GATA3/Th2 axis)", "down"),
    ]
    results = []
    for node, note, direction in checks:
        before, after = no_mir[node], mir155[node]
        ok = (after > before) if direction == "up" else (after < before)
        results.append({
            "node": node, "expectation": note,
            "no_miRNA_pct": before, "miR155_only_pct": after,
            "matches_literature": ok,
        })
    return pd.DataFrame(results)


def main():
    print("=== Part A: Th-polarisation test (no miRNA, base biology only) ===")
    df_a = part_a_th_polarisation()
    print(df_a)
    df_a.to_csv(os.path.join(OUT_DIR, "validation_th_polarisation.csv"))

    print("\n=== Part A: pass/fail ===")
    checks_a = evaluate_part_a(df_a)
    print(checks_a.to_string(index=False))
    checks_a.to_csv(os.path.join(OUT_DIR, "validation_th_polarisation_checks.csv"), index=False)

    print("\n=== Part B: miR-155-5p positive control (vs. main steady_state_comparison.csv) ===")
    main_csv = os.path.join(OUT_DIR, "steady_state_comparison.csv")
    df_b = part_b_mir155_positive_control(main_csv)
    print(df_b.to_string(index=False))
    df_b.to_csv(os.path.join(OUT_DIR, "validation_mir155_control.csv"), index=False)

    summary = {
        "th_polarisation_all_pass": bool(checks_a["pass"].all()),
        "th_polarisation_n_pass": int(checks_a["pass"].sum()),
        "th_polarisation_n_total": int(len(checks_a)),
        "mir155_all_match": bool(df_b["matches_literature"].all()),
    }
    with open(os.path.join(OUT_DIR, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== Summary ===")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
