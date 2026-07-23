# Revisiting the T-cell model in Python

A from-scratch Python rebuild of a piece of my PhD's logical (Boolean/
multi-valued) model of CD4+ T-cell activation and differentiation, extended
with miR-34c-5p and miR-155-5p regulation -- originally built in GINsim/R for
my 2019 thesis at the University of Lisbon ("sncRNA regulatory networks in
T cell activation and viral response").

The original rule table (76 nodes, 89 formulas) is transcribed verbatim from
the thesis into `src/rules_data.py`. See `MODEL_NOTES.md` for exactly which
modelling choices were made in the rebuild, and where results diverge from
the 2019 findings.

## Structure

```
src/
  rules_data.py     digitized rule table (verbatim from the thesis)
  logic_engine.py   tokenizer/parser/evaluator + asynchronous simulator
  run_analysis.py   builds the 4 miRNA conditions, runs the ensemble,
                     writes tables to output/ and figures to figures/
output/              CSV tables + run metadata (JSON)
figures/              network_graph.png, trajectories.png, phenotype_bars.png
MODEL_NOTES.md        modelling decisions & caveats
```

## Running it

```bash
pip install -r requirements.txt
python src/run_analysis.py
```

Takes under a minute; regenerates everything in `output/` and `figures/`.

## What it does

1. Parses the thesis's own logical notation (`&` AND, `|` OR, `!` NOT,
   `NODE:k` threshold reference) with a small recursive-descent parser.
2. Builds 4 versions of the network: no miRNA regulation, miR-34c-5p only,
   miR-155-5p only, and both -- by including/excluding their transcriptional
   rules.
3. Runs 400 independent asynchronous simulations per condition from a
   TCR + IL2-stimulated naive CD4+ T-cell initial state.
4. Compares master transcription factor activation frequency (Th1/Th2/Th17/
   iTreg/Tfh/Th22/Th9 markers) across the four conditions, and plots node
   activity trajectories over time.

## Dependencies

`numpy`, `pandas`, `matplotlib`, `networkx` -- see `requirements.txt`. No
external logical-modelling software (GINsim, MaBoSS, BoolNet) is required;
the simulation engine is a self-contained ~150-line reimplementation.
