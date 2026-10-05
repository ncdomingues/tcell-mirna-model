# Revisiting the T-cell model in Python

A from-scratch Python rebuild of a piece of my PhD's logical (Boolean/
multi-valued) model of CD4+ T-cell activation and differentiation, extended
with miR-34c-5p and miR-155-5p regulation -- originally built in GINsim/R for
my 2019 thesis at the University of Lisbon ("sncRNA regulatory networks in
T cell activation and viral response").

The original rule table (76 nodes, 89 formulas) is transcribed verbatim from
the thesis into `src/rules_data.py`. See `MODEL_NOTES.md` for exactly which
modelling choices were made in the rebuild, where results diverge from the
2019 findings, and `VALIDATION.md` for what's actually been checked (23 unit
tests + a biology-only sanity check independent of the miRNA question).

## Structure

```
src/
  rules_data.py       digitized rule table (verbatim from the thesis)
  logic_engine.py     tokenizer/parser/evaluator + asynchronous simulator
  run_analysis.py     builds the 4 miRNA conditions, runs the ensemble,
                       writes tables to output/ and figures to figures/
  validate_biology.py Th-polarisation + miR-155-5p literature sanity checks
tests/
  test_logic_engine.py  23 unit tests: parser, evaluator, digitization integrity
walkthrough.ipynb      annotated, executed notebook walking through every
                       piece of code above, cell by cell, with live output
output/                CSV tables + run/validation metadata (JSON)
figures/                network_graph.png, trajectories.png, phenotype_bars.png
MODEL_NOTES.md          modelling decisions & caveats
VALIDATION.md           what's been validated, and how
```

## Running it

```bash
pip install -r requirements.txt
python -m unittest tests.test_logic_engine -v   # 23 tests, code correctness
python src/run_analysis.py                       # main 4-condition comparison
python src/validate_biology.py                   # biological sanity checks
```

Each takes well under a minute; the latter two regenerate everything in
`output/` and `figures/`.

To read the annotated walkthrough (explains the parser, the multi-valued
semantics, the simulator, and the 4-condition experiment with live code and
output), open `walkthrough.ipynb` in Jupyter, or view it directly on GitHub
(it renders with all outputs already saved).

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
5. `ml_predict.py`: trains a classifier to predict eventual phenotype from a
   trajectory's own partial (not-yet-converged) state, and separately asks
   whether the same task is even answerable from what the original 2019
   R/GINsim analysis ever published (it isn't -- see the script's docstring).

## Dependencies

`numpy`, `pandas`, `matplotlib`, `networkx`, `scikit-learn` -- see
`requirements.txt`. No external logical-modelling software (GINsim, MaBoSS,
BoolNet) is required; the simulation engine is a self-contained ~150-line
reimplementation.

## R re-analysis (`r_analysis/`)

An independent second implementation of the same engine, in R + ggplot2,
used to (1) cross-validate the Python engine against an unrelated codebase,
and (2) properly re-run two real historical experiments found in the
original 2016-2017 `ThesisNDPhD/Modelling/` working directory -- the
miR-34c-5p transcription-factor hypothesis testing and the Th1/Th2/Th17/Treg
stimulation-condition sweeps -- this time with proper ensembles and
confidence intervals. See `r_analysis/R_ANALYSIS.md` for the full write-up
and findings. Requires R (>=4.0) with `tidyverse`/`ggplot2`/`jsonlite`.
