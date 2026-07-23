# Modelling notes

Honest documentation of the choices, ambiguities, and deliberate simplifications
made in rebuilding the logical model. The goal of this project was fidelity to
the *published rule table*, not a bit-for-bit reproduction of the 2019 GINsim
steady-state numbers -- those used a different simulation engine (GINsim +
MaBoSS + R) and, in places, an exact/approximate steady-state search rather
than ensemble asynchronous sampling. Where this rebuild's numbers diverge from
the thesis's Table 3/4, that is expected and discussed in the case study.

## What was transcribed exactly

The full 76-node, 89-formula rule table (`src/rules_data.py`) is a verbatim
transcription of Table 5 in `Results_PhDthesis.docx`, including its base
models (Saez-Rodriguez et al. 2007 for proximal TCR signalling; Naldi et al.
2010 and Abou-Jaoudé et al. 2015 for Th differentiation) and the miR-34c-5p /
miR-155-5p extensions.

## Multi-valued threshold semantics

For nodes with more than one target value (e.g. `NFKB`, `TP53`, `GATA3`,
`TBX21`, `miR-34c-5p`), a node takes the **highest** value whose formula
evaluates true (checked in decreasing order), otherwise 0. This is the
standard reading of GINsim-style logical parameter tables and matches every
example in the table (e.g. `ZAP70=2` requires `LCK & !miR-34c-5p`, while
`ZAP70=1` requires only `LCK` -- i.e. miR-34c-5p partially, not fully,
represses ZAP70).

## Operator precedence

The thesis's legend defines `&` = AND, `|` = OR, `!` = NOT, but not explicit
precedence. This rebuild uses `NOT` > `AND` > `OR`, i.e. exactly how Python
itself binds `not`/`and`/`or` -- the natural reading of `A & B | C & !D` as
`(A&B) | (C & !D)`.

## Two formulas preserved "as ambiguous"

- **SIRT1**: `FOXO3 | (FOXO3 & CREB) & !TP53 & !miR-155-5p`. Under NOT>AND>OR
  precedence this parses as `FOXO3 | ((FOXO3&CREB) & !TP53 & !miR-155-5p)`,
  meaning SIRT1 turns on whenever FOXO3 alone is active, regardless of TP53 --
  even though the thesis text says "TP53 inhibits it." This may be a missing
  parenthesis in the original table. It was **not** silently corrected;
  it is transcribed and evaluated exactly as published.
- **MAF**: `STAT5 | STAT3 & !(TGFB & miR-155-5p)` similarly only gates the
  STAT3 route with the TGFB/miR-155-5p brake, not the STAT5 route, even
  though the surrounding text implies both.

## Operationalising "TCR + IL2 stimulation"

The thesis's in vitro protocol is anti-CD3 + anti-CD28 stimulation of naive
CD4+ T cells with exogenous IL2. In the rule table, `IL2` is itself a
regulated node (auto-produced via NFAT/NFKB once TCR signalling begins) and
there is no separate "exogenous IL2" input. Rather than inventing a new node,
this rebuild sets the receptor/co-receptor/coactivator inputs
(`APC-Antigen`, `CD3`, `CD45RA`, `CGC`, `CREBBP`) to 1 for the full run and
lets IL2 (and everything downstream) evolve through the model's own rules --
consistent with how the thesis describes the model's built-in
NFAT/NFKB -> IL2 -> STAT5 feed-forward loop.

## Simulation method

The original work used GINsim's asynchronous non-deterministic update scheme
plus MaBoSS for reachability/steady-state analysis. This rebuild uses its own
from-scratch asynchronous engine (`src/logic_engine.py`): each of 400
independent trajectories updates one randomly-chosen rule-governed node per
step (order reshuffled every sweep), for 10 full sweeps, starting from the
all-OFF state plus the fixed inputs above. "Activation frequency" per node =
the fraction of the 400 final states in which that node's value was >= 1.
This is a Monte Carlo *approximation* of the steady-state distribution, not
an exhaustive attractor search -- a deliberate, documented trade-off to keep
the rebuild dependency-free and inspectable in pure Python.

## Where results diverge from the 2019 thesis

The original discussion emphasises a **Th17/iTreg** shift when both miRNAs
are present. This rebuild instead finds **miR-155-5p** as the dominant driver
of phenotype change on its own (suppressing Th2/Tfh/Th22/Th9 master
regulators, raising Th1/iTreg), with **miR-34c-5p**'s individual effect much
closer to the no-miRNA baseline -- but the *combination* of both miRNAs
suppresses the Th2/Tfh/Th22/Th9 branch further than miR-155-5p alone, echoing
the thesis's own hypothesis of a cooperative interaction between the two
miRNAs, even though the specific destination phenotype differs. Likely
sources of the discrepancy: the asynchronous ensemble-average method used
here versus GINsim/MaBoSS's steady-state enumeration in 2019, and possibly
transcription artifacts in a couple of formulas (see above). This divergence
is reported honestly in the case study rather than tuned away.
