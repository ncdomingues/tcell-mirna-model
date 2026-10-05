# -*- coding: utf-8 -*-
"""Re-runs the main 4-condition analysis and emits every watched node's
activation frequency WITH a Wilson 95% confidence interval, plus the
miR-155-5p positive control as a proper table -- addressing the fact that
run_analysis.py's CSVs report bare point estimates with no uncertainty."""
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_analysis import run_all, WATCHLIST, MASTER_TFS
from stats_utils import wilson_ci

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "output")

ORDER = ["no miRNA", "miR-155-5p only", "miR-34c-5p only", "both"]


def main():
    conditions, results = run_all()

    all_nodes = list(dict.fromkeys(list(MASTER_TFS.keys()) + WATCHLIST))
    rows = []
    for node in all_nodes:
        row = {"node": node}
        for cond in ORDER:
            finals = results[cond]["finals"]
            n = len(finals)
            k = sum(1 for f in finals if f.get(node, 0) >= 1)
            pct, lo, hi = wilson_ci(k, n)
            row[f"{cond}_pct"] = pct
            row[f"{cond}_ci"] = f"{lo}-{hi}"
        rows.append(row)
    df = pd.DataFrame(rows).set_index("node")
    df.to_csv(os.path.join(OUT_DIR, "full_results_with_ci.csv"))
    print(f"Full results with 95% CI written for {len(all_nodes)} nodes x {len(ORDER)} conditions")
    print(df.head(10).to_string())

    # miR-155-5p positive control, with CI, as a standalone table
    control_rows = []
    for node, label, direction in [
        ("TBX21", "Th1-promoting (O'Connell et al. 2010)", "up"),
        ("GATA3", "Th2-repressing", "down"),
    ]:
        finals_none = results["no miRNA"]["finals"]
        finals_155 = results["miR-155-5p only"]["finals"]
        n = len(finals_none)
        k0 = sum(1 for f in finals_none if f.get(node, 0) >= 1)
        k1 = sum(1 for f in finals_155 if f.get(node, 0) >= 1)
        p0, l0, h0 = wilson_ci(k0, n)
        p1, l1, h1 = wilson_ci(k1, n)
        matches = (p1 > p0) if direction == "up" else (p1 < p0)
        control_rows.append({
            "node": node, "expectation": label,
            "no_miRNA_pct": p0, "no_miRNA_ci": f"{l0}-{h0}",
            "miR155_pct": p1, "miR155_ci": f"{l1}-{h1}",
            "matches_literature": matches,
        })
    df_control = pd.DataFrame(control_rows)
    df_control.to_csv(os.path.join(OUT_DIR, "mir155_control_with_ci.csv"), index=False)
    print("\n", df_control.to_string(index=False))


if __name__ == "__main__":
    main()
