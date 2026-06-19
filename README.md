# local-bending-probe

**An empirical demonstration that local sequence alone does not contain enough
information to predict protein backbone bending — which helps explain why
successful structure-prediction models (AlphaFold and kin) depend on long-range
interactions.**

This began as an attempt to predict mutation-induced backbone bending from
*local* sequence — a lightweight, interpretable alternative to heavy
structure-prediction models. **It failed.** That failure is the result, and the
more interesting outcome:

```
Bending is real, reproducible, and distributed                          (Gate 1)
        │
        ▼
Local sequence cannot predict it
   AUC ≈ 0.52 (cores), ≈ 0.59 (loops, powered) — every CI touches 0.5,
   evaluated on completely unseen protein families                      (Gates 2–4)
        │
        ▼
Reintroducing crude NON-local (3D contact) information shifts it back in
   the predicted direction — core subset 0.48 → 0.58
   (suggestive evidence for the mechanism — NOT a certified predictor)   (Gate 5)
```

The experiment effectively asks **"where does the information that determines
bending live?"** and answers: **not primarily in local sequence.** Remove the
non-local information and predictive signal collapses to chance; reintroduce
even a crude version of it and the signal moves back — and it moves back
specifically in the protein core, where packing and tertiary contacts dominate,
exactly as protein physics predicts.

> **This is not a failed AlphaFold competitor.** It is a controlled experiment
> that *locates the information* — and in doing so illustrates why successful
> structure-prediction models rely on long-range and evolutionary signal rather
> than local sequence motifs.

## Why a 0.52 is a result, not a non-result

A predictor scoring AUC 0.52 sounds like nothing — until you ask *how* it was
measured. **0.52 under rigorous, leakage-controlled evaluation on unseen protein
families is more informative than 0.80 under questionable splits**, because the
first is close to the truth. The controls are the contribution:

- **136,961** local windows from **568 non-redundant chains**
- **Family-level holdout** + **≤30% sequence-identity culling** — no fragment
  leakage; tested on entirely unseen families (where most biological-ML results
  quietly lose their performance)
- **Bootstrap 90% confidence intervals** on every estimate — and the discipline
  of reading them: a CI touching 0.5 is reported as *"not distinguishable from
  chance,"* never "weak signal"
- **Six independent falsification gates**, each with an explicit go/no-go
- Crystallographic-artifact QC, empirically-derived noise floors, and a
  metric-correction (the 5.7× differencing error-cancellation that fixed an
  over-conservative test criterion)

Full numbers — **every estimate as point + 90% CI + n** — are in
**[RESULTS.md](RESULTS.md)**.

## The falsification ladder

| Gate | Question | Result | Verdict |
|---|---|---|---|
| 0 | Is the bending metric sound? | exact 0–170°; hairpin reads 170° not 10° | ✅ sound |
| 1 | Does the bending signal exist? | 29% of mutations move >2σ, across 38 residues | **GO\*** |
| 2 | Can local sequence predict it? | AUC 0.524, Spearman 0.08 (n=72) | ✗ no |
| 2b | What's locally addressable? | residual 7%, mostly artifact | ≈ 0 |
| 3 | Signal in loops (T4L)? | AUC 0.70 but **CI [0.49, 0.89]**, n=10 | underpowered |
| 4 | Replicates when powered? | **AUC 0.591, CI [0.49, 0.69]**, n=33 | not distinguishable from chance |
| 5 | Does 3D context rescue it? | cores 0.48→0.58 | mechanism ✓ (suggestive) / predictor ✗ |

\* **GO with a caveat:** 29% is plausibly a *ceiling, not a floor* — T4 lysozyme
is among the most mutation-tolerant proteins known and is core-biased, and sparse-
window floors inflate the fraction. Details in [RESULTS.md](RESULTS.md).

## What this demonstrates (for someone reviewing the repo)

```
Hypothesis  →  designed falsification tests  →  hypothesis failed  →  explained why
```

That chain demonstrates **experimental design and statistical honesty**, not just
model-training. Many people can train a model to 0.59; far fewer design an
experiment that meaningfully answers *where the information lives* — and then
report the negative answer without inflating it. That discipline is the point of
this repo.

## Relation to AlphaFold (a framing, not a benchmark)

There is deliberately **no head-to-head accuracy number against AlphaFold** — it
would be a category error. AlphaFold predicts *global 3D structure from
sequence + MSA*; this predicts *local bending change from local sequence*.
Different inputs, outputs, task. The honest connection is conceptual: by removing
non-local information and watching predictive signal collapse, this project is the
**negative-space proof** of why AlphaFold's reliance on long-range/coevolutionary
information is necessary rather than incidental.

## Reproduce

```bash
pip install numpy scipy scikit-learn matplotlib   # no gemmi/BioPython needed
python bending_metric.py            # Gate 0: metric self-test
python feasibility_t4l.py           # Gate 1: data feasibility (downloads T4L PDBs)
python gate2_model_feasibility.py   # Gate 2: local model
python mover_composition.py         # Gate 2b: what's locally addressable
python loop_gate.py                 # Gate 3: loops (T4L)
python powered_loop_gate.py         # Gate 4: powered loop replication
python gate3_3d.py                  # Gate 5: 3D contact-environment ingredient
```
Structures are fetched from RCSB and cached locally; runs are CPU-only,
minutes-scale.

## Limitations (honest, up front)

Small mover counts throughout (10–72 → wide CIs); single-family bias (T4L cores,
a human-lysozyme-dominated loop pool, SNase excluded); a crude 3D feature
(contact-AA composition only, no geometry/energy); single-sequence inputs (no MSA
profiles); a Cα-only bending definition; observational PDB mutants rather than
designed. None are hidden — they bound exactly what the conclusions can claim.

## Files

| File | Role |
|---|---|
| `bending_metric.py` | Canonical Cα bending metric (the one reusable component) |
| `feasibility_t4l.py` | Gate 1 — data feasibility |
| `gate2_model_feasibility.py` | Gate 2 — local prediction model |
| `mover_composition.py` | Gate 2b — mover composition |
| `loop_gate.py` | Gate 3 — loop-specific test |
| `powered_loop_gate.py` | Gate 4 — powered replication |
| `gate3_3d.py` | Gate 5 — 3D contact-environment |
| `RESULTS.md` | Full results, every estimate as point + CI + n |
| `SPEC_bending_and_pairs.md` | Metric & benchmark spec |
| `*.png` | Per-gate figures |

---

*A negative result, reported honestly. The phenomenon is real; local-sequence
predictability of it is not — and locating the missing information in non-local
structure is the actual finding.*
