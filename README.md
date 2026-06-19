# local-bending-probe

**Can protein backbone bending be predicted from *local* sequence — a lightweight,
interpretable alternative to AlphaFold-style models? A gated, falsification-driven
investigation.**

> **Short answer: no — and the *why* is the contribution.** Mutation-induced backbone
> bending is real, but it is not predictable from local sequence (AUC ≈ 0.52, and ≈ 0.59
> in loops even with proper power — every confidence interval touches chance). Adding 3D
> context nudges the result in the theoretically-predicted direction, confirming the
> cause is **tertiary**, not local. This repo is a controlled demonstration of *why*
> AlphaFold needs the non-local information it uses.

---

## Motivation

AlphaFold solved global structure prediction with deep learning + evolutionary
couplings (MSAs) — accurate, but heavy and a black box. I wanted to test, **myself and
from scratch**, whether a deliberately opposite approach could carve out a useful niche
on a *sub-problem* (local backbone bending): **lightweight, interpretable, local, runs
on a laptop.** The hypothesis (an explicit "conditional discriminative pattern mining"
architecture): *given similar baseline geometry, specific amino-acid combinations push a
segment toward bending or not, and those rules are learnable and interpretable.*

Rather than assume it would work, I built a ladder of **falsification gates** and let the
data decide at each step.

## Relation to AlphaFold (read this — it's a framing, not a benchmark)

There is **no head-to-head accuracy number against AlphaFold here, and there shouldn't
be** — it would be a category error. AlphaFold predicts *global 3D structure from
sequence+MSA*; this predicts *local bending change from mutations using local sequence*.
Different inputs, outputs, and task.

The honest comparison is conceptual: **AlphaFold works precisely because it exploits the
non-local / coevolutionary information this project deliberately excludes.** By removing
that information and measuring the result, this project is the *negative-space proof* of
AlphaFold's central design choice — it shows, on real data, that local sequence does not
contain what determines backbone geometry.

## The falsification ladder

| Gate | Question | Result | Verdict |
|---|---|---|---|
| 0 | Is the bending metric sound? | exact on 0–170°; hairpin reads 170° not 10° | ✅ sound |
| 1 | Does the bending signal exist? | 29% of mutations move >2σ, across 38 residues | **GO\*** |
| 2 | Can local sequence predict it? | **AUC 0.524**, Spearman 0.08 (n=72) | ✗ no |
| 2b | What's locally addressable? | residual 7%, mostly artifact | ≈ 0 |
| 3 | Signal in loops (T4L)? | AUC 0.70 but **CI [0.49, 0.89]**, n=10 | underpowered |
| 4 | Does it replicate, powered? | **AUC 0.591, CI [0.49, 0.69]**, n=33 | not distinguishable from chance |
| 5 | Does 3D context rescue it? | 0.52→0.59 overall; cores 0.48→0.58 | mechanism ✓ / predictor ✗ |

\* **GO with a caveat:** 29% is plausibly a *ceiling, not a floor* — T4 lysozyme is among
the most mutation-tolerant proteins known and is core-biased, and sparse-window floors
inflate the fraction. See [RESULTS.md](RESULTS.md).

> Full numbers — each as **point estimate + 90% CI + n** — are in **[RESULTS.md](RESULTS.md)**,
> which also explains why "local-only" appears as 0.524 / 0.516 / 0.477 in different places
> (QC filter + subset definitions, not inconsistency).

## What the evidence actually supports

- **Local sequence → bending: not distinguishable from chance** anywhere tested. Cores
  0.52, loops 0.59 (powered) — every CI includes 0.50.
- **The cause is tertiary, confirmed by the 3D gate.** Adding the contact environment
  lifts the *core* subset (packing-dominated) 0.477 → 0.584 — right direction, right
  place. This is strong **evidence for the mechanism**; it is **not** a validated
  predictor (overall CI [0.52, 0.66] barely clears chance, crude feature).
- These two are kept strictly separate throughout: *why it fails* (established) vs
  *a working predictor* (never reached).

## What this project demonstrates

- Designing **falsification gates** instead of confirming a hypothesis.
- Statistical honesty under self-scrutiny: bootstrap CIs, leakage control (sequence
  culling + held-out families), empirical noise floors, crystallographic-artifact QC,
  and catching/correcting an over-conservative test criterion (the 5.7× differencing
  error-cancellation).
- Reading a **negative result** correctly and not overclaiming — the hardest and most
  valuable discipline in computational science.

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
Structures are fetched from RCSB and cached locally; runs are CPU-only and minutes-scale.

## Limitations (honest, up front)

Small mover counts throughout (10–72 → wide CIs); single-family bias (T4L cores, a
human-lysozyme-dominated loop pool, SNase excluded); a crude 3D feature (contact-AA
composition only, no geometry/energy); single-sequence inputs (no MSA profiles); a
Cα-only bending definition; and observational PDB mutants rather than designed ones.
None of these are hidden — they bound exactly what the conclusions can claim.

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
predictability of it is not. That distinction — and the discipline of holding to it — is
the point.*
