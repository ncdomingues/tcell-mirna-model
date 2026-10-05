# -*- coding: utf-8 -*-
"""
Convergence check: does the ensemble size (400 runs, used throughout the main
analysis) actually give a stable estimate, or would more runs change the
reported percentages? Re-runs the "both miRNAs" condition at increasing
ensemble sizes and tracks how each master-TF estimate (with its Wilson 95%
CI) moves, for the 7 Th master regulators.
"""
import os
import sys
import json
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rules_data import RULES, ALIASES
from logic_engine import LogicalModel
from stats_utils import wilson_ci

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "output")
FIG_DIR = os.path.join(HERE, "..", "figures")

MASTER_TFS = {"TBX21": "Th1", "GATA3": "Th2", "RORC": "Th17", "FOXP3": "iTreg",
              "BCL6": "Tfh", "AHR": "Th22", "SPI1-PU1": "Th9"}
RUN_SIZES = [50, 100, 200, 400, 800, 1600]
N_SWEEPS = 10
SEED = 42


def base_initial_state(model):
    state = {n: 0 for n in model.all_nodes}
    for inp in ("APC-Antigen", "CD3", "CD45RA", "CGC", "CREBBP"):
        if inp in state:
            state[inp] = 1
    return state


def main():
    model = LogicalModel(RULES, ALIASES)
    init = base_initial_state(model)

    rows = []
    for n_runs in RUN_SIZES:
        finals, _traj = model.ensemble_run(init, n_sweeps=N_SWEEPS, n_runs=n_runs, seed=SEED)
        n = len(finals)
        for tf, label in MASTER_TFS.items():
            k = sum(1 for f in finals if f.get(tf, 0) >= 1)
            pct, lo, hi = wilson_ci(k, n)
            rows.append({"n_runs": n_runs, "node": tf, "subtype": label,
                         "pct": pct, "ci_lo": lo, "ci_hi": hi, "ci_halfwidth": round((hi - lo) / 2, 1)})
        print(f"n_runs={n_runs} done")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_DIR, "robustness_convergence.csv"), index=False)
    print(df[df["n_runs"].isin([50, 400, 1600])].to_string(index=False))

    at_400 = df[df["n_runs"] == 400].set_index("node")["ci_halfwidth"]
    at_1600 = df[df["n_runs"] == 1600].set_index("node")["ci_halfwidth"]
    summary = {
        "run_sizes_tested": RUN_SIZES,
        "max_ci_halfwidth_at_400_runs_pct": float(at_400.max()),
        "max_ci_halfwidth_at_1600_runs_pct": float(at_1600.max()),
        "max_estimate_shift_400_to_1600_pct": round(float(
            (df[df.n_runs == 1600].set_index("node")["pct"] - df[df.n_runs == 400].set_index("node")["pct"]).abs().max()
        ), 1),
    }
    with open(os.path.join(OUT_DIR, "robustness_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))

    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=True)
    for ax, (tf, label) in zip(axes.flat, MASTER_TFS.items()):
        sub = df[df.node == tf]
        ax.errorbar(sub.n_runs, sub.pct, yerr=[sub.pct - sub.ci_lo, sub.ci_hi - sub.pct],
                     marker="o", capsize=3, color="#1F6F63")
        ax.set_xscale("log")
        ax.set_title(f"{tf} ({label})", fontsize=10)
        ax.grid(alpha=0.3)
    axes.flat[-1].axis("off")
    for ax in axes[1, :4]:
        ax.set_xlabel("ensemble size (n runs, log scale)")
    axes[0, 0].set_ylabel("% active (95% CI)")
    axes[1, 0].set_ylabel("% active (95% CI)")
    fig.suptitle("Convergence check: master-TF estimate vs. ensemble size ('both miRNAs' condition)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "robustness_convergence.png"), dpi=170, facecolor="#F3F4EF")
    plt.close(fig)
    print("\nFigure written to figures/robustness_convergence.png")


if __name__ == "__main__":
    main()
