# -*- coding: utf-8 -*-
"""Exports RULES/ALIASES to TSV so the R port loads byte-identical rule text
(avoids any risk of manual-transcription drift between the two engines)."""
import csv
import os
from rules_data import RULES, ALIASES

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "r_analysis", "data")
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(OUT_DIR, "rules.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["node", "value", "formula", "description"])
    for node, value, formula, desc in RULES:
        w.writerow([node, value, formula, desc])

with open(os.path.join(OUT_DIR, "aliases.tsv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["alias", "canonical"])
    for a, c in ALIASES.items():
        w.writerow([a, c])

print(f"Wrote {len(RULES)} rules and {len(ALIASES)} aliases to {OUT_DIR}")
