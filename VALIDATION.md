# Validation

Two separate questions, tested separately:

1. **Is the code correct?** Does the parser/evaluator/simulator actually implement
   the notation and semantics it claims to? -- `tests/test_logic_engine.py`
2. **Is the digitized network biologically sane?** Independent of the miRNA
   hypothesis entirely: does the transcribed rule table reproduce textbook
   CD4+ T-helper biology? -- `src/validate_biology.py`

Run both yourself:

```bash
python -m unittest tests.test_logic_engine -v
python src/validate_biology.py
```

## 1. Code correctness -- 23/23 tests pass

`tests/test_logic_engine.py` covers:

- **Tokenizer**: hyphenated identifiers (`AKT1-PKB`, `miR-34c-5p`), threshold
  suffixes (`IL2R:2`) stay attached as one token, whitespace is ignored.
- **Parser precedence**: `NOT > AND > OR` is verified against the AST shape
  directly (not just the boolean result), including on a real formula from
  the table (`IL2`'s value-2 rule) and on the one formula flagged in
  `MODEL_NOTES.md` as ambiguously grouped (`SIRT1`) -- locking in the
  documented reading so a future refactor can't silently change it.
- **Evaluator**: AND/OR/NOT, bare references (truthy at >=1) vs. threshold
  references (`NODE:k`), missing nodes defaulting to 0, alias resolution
  (`AKT` -> `AKT1-PKB`).
- **Multi-valued threshold semantics**: a node takes the highest satisfied
  value, on a minimal synthetic 2-rule example.
- **Structural integrity of the digitization itself**: exactly 89 rule rows
  and 76 rule-governed nodes (matching Table 5), exactly the 5 expected
  input nodes (`APC-Antigen`, `CD3`, `CD45RA`, `CGC`, `CREBBP`), both miRNAs
  confirmed rule-governed (not accidentally treated as inputs), every rule
  parses without error.

## 2. Biological validation

### Part A -- classic Th-polarisation test (no miRNA at all)

The base model this rebuild extends (Naldi et al. 2010) was originally
validated by showing that exposing it to known Th-polarising cytokines
produces the correct T-helper master regulator. This rebuild's rule table
was tested the same way, with the two miRNAs entirely removed from the
network -- this section tests nothing about miR-34c-5p or miR-155-5p, only
whether the transcribed core network behaves like real T-cell biology.

Each cytokine was held as a fixed exogenous input (300-run ensemble, 10
sweeps, on top of baseline TCR+IL2 stimulation):

| Condition | Expected | Result |
|---|---|---|
| + IL12 | Th1 (TBX21) | **Pass** -- TBX21 0% &rarr; 64%, clear top responder |
| + IL4 | Th2 (GATA3) | **Pass** -- GATA3 42.7% &rarr; 47%, top responder |
| + TGFB + IL6 | Th17 (RORC) | **Pass** -- RORC 0% &rarr; 73.7%, top responder |
| + TGFB alone | iTreg (FOXP3) | **Fail** -- FOXP3 *drops* 30.3% &rarr; 19.3%; RORC dominates instead (73.7%) |

**3 of 4 pass cleanly.** The iTreg failure has a specific, traceable cause,
not a mysterious one: in the digitized `STAT3` rule,

```
STAT3 = 1 if (TGFBR | IL6R | IL23R | IL1R) & !SIRT1 & !miR-155-5p & !miR-34c-5p
```

`TGFBR` alone is one of four alternative activators of `STAT3` (they're
OR'd), so supplying TGFB by itself is already sufficient to drive the
STAT3 -> RORC (Th17) axis in this network, with no additional requirement
for IL6. Structurally, this specific integrative model doesn't cleanly
separate "TGFB alone" from "TGFB + IL6" the way a dedicated Th17-vs-iTreg
model would. This isn't a transcription error -- it's the literal rule as
published in Table 5 -- and it's consistent with how the thesis itself
describes this particular model: it was purpose-built to simulate one
specific condition (TCR + IL2 stimulation with autocrine cytokine release),
not to be pushed with arbitrary exogenous cytokine cocktails outside that
design envelope. Reported here rather than hidden.

### Part B -- miR-155-5p positive control

Before trusting anything the model predicts about the poorly-characterised
miR-34c-5p, the thesis used the immune system's best-studied miRNA,
miR-155-5p, as a sanity check against known literature. Same check, run
against this rebuild's main 400-run, 4-condition results:

| Node | Expectation | no miRNA | miR-155-5p only | Match? |
|---|---|---|---|---|
| TBX21 | Th1-promoting (O'Connell et al. 2010) | 27.0% | 36.2% | **Yes** |
| GATA3 | Th2-repressing | 44.0% | 15.5% | **Yes** |

Both directions match established literature.

## Bottom line

The code does what it says (23/23 unit tests). The digitized network
reproduces 3 of 4 textbook Th-polarisation outcomes with a fully explained
structural cause for the fourth, and both directions of the miR-155-5p
positive control match literature. That's a reasonable basis for trusting
what Section 8 of `walkthrough.ipynb` reports about miR-34c-5p -- with the
iTreg caveat above carried forward honestly rather than smoothed over.
