# Results

A gate-by-gate record of the investigation. The discipline of this document:
**every retrieval result is reported as point estimate + 90% bootstrap CI + number
of movers (n), in one breath.** Where a CI includes 0.50, the honest statement is
*"not distinguishable from chance at this sample size"* — never "weak signal."

---

## Read this first: why some numbers differ between runs

You will see the **local-only** model at AUC **0.524**, **0.516**, and **0.477**.
These are not inconsistencies; they are three honestly-different measurements:

| Reported | What it actually is |
|---|---|
| 0.524 (Gate 2) | all 248 T4L single-mutant pairs, **before** the QC filter |
| 0.516 (Gate 5) | the same set **after** QC removed 2 crystallographic artifacts (V111M −32.6°, V111I −19.5°) → 246 pairs / 70 movers |
| 0.477 (Gate 5) | the **core-only** subset (helix+sheet, n=225) — a different, harder population |

So the mover-count wobble (72 vs 70, 248 vs 246) is the QC step removing two artifacts,
**not** noise. Beyond these defined differences, point estimates vary by **±0.02–0.03**
across runs from sampling/seed variation. That wobble is small — but it is *not* small
relative to the effect being chased. The mover counts throughout are 10–72; a method
whose claimed signal is smaller than its own run-to-run noise is not a usable method.
**That the noise and the "signal" are the same size is itself part of the finding.**

---

## Gate 0 — Metric validation (`bending_metric.py`)

**Claim:** the bending metric is correct and free of sign/orientation bugs.
**Numbers:** synthetic windows with known bends recovered exactly — 0°, 30°, 90°,
135°, 160°, 170°. A 170° chain reversal reads **170°, not 10°** (the orientation
guard holds; a sign flip would have made a hairpin look straight).
**Limit:** Cα-only, one geometric definition of "bending." Defensible and intrinsic
(superposition-free), but it is *a* definition, not *the* definition.

## Gate 1 — Does the signal exist? (`feasibility_t4l.py`)

**Claim:** mutation-induced local bending is real, reproducible, and distributed.
**Numbers:** 248 scorable T4L single-mutant windows. Per-window noise floor **0.98°**
(WT crystal-to-crystal), independently corroborated by **0.75°** (same variant
re-crystallized). **29% of mutations exceed 2σ (72/248)**, spread across **38 distinct
residues, only 6% in the mobile hinge** → distributed, not a single-region artifact.
**Honest limit (carry this everywhere):** 29% is plausibly a **ceiling, not a floor.**
T4 lysozyme is among the most mutation-tolerant proteins known and its mutagenesis is
core-biased; and per-window floors estimated on sparsely-sampled windows *inflate* the
above-floor fraction. "GO" means *"there is a signal worth studying,"* not *"29% of
mutations generically move backbones."*

## Gate 2 — Can local sequence predict it? (`gate2_model_feasibility.py`)

**Claim:** a local-sequence model **cannot** predict which mutations move the backbone.
**Numbers:** HistGradientBoosting trained on **136,961 windows / 568 non-redundant
chains** (≤30% identity, X-ray ≤2.0 Å). Absolute-bending RMSE **30.55°** (std 37.97°,
~20% variance reduction). Differenced mover retrieval on 248 T4L pairs:
**AUC 0.524, AP 0.322 (base 0.290), Spearman 0.077, n=72 movers.**
**Methodological highlight (worth pausing on):** the predicted-Δ error is **7.58°, not
the 43.2° a naïve √2×RMSE bound predicts — a 5.7× error cancellation**, because the WT
and mutant inputs are near-identical and their errors cancel in the difference. The
original go/no-go criterion (√2×RMSE) was therefore *over-conservative*; it was caught
and replaced with a direct, leakage-free retrieval measurement. The model fails **not**
from differencing noise but from **insensitivity** — it barely responds to a
single-residue change in a way that tracks reality (Spearman 0.08).
**Limit:** single family (T4L); n=72 movers.

## Gate 2b — What is locally addressable? (`mover_composition.py`)

**Claim:** essentially nothing beyond textbook effects is locally addressable.
**Numbers:** of 72 movers — Pro/Gly involved in **3%** (vs 2% of non-movers; OR 1.23,
p 0.56 → *not* enriched). Non-local contacts present in **93%** of movers (vs 95% of
non-movers; OR 0.64 → ubiquitous, *not* enriched). Residual (neither): **5/72 (7%)**,
of which the two largest were crystallographic artifacts.
**Limit:** "non-local contact" here is a crude binary flag; in a packed protein it is
true almost everywhere, so its non-enrichment is informative but coarse.

## Gate 3 — Do loops rescue it? (T4L) (`loop_gate.py`)

**Claim:** loops *looked* like the one surviving home for local signal — but T4L
cannot certify it.
**Numbers:** movers enriched in loops (**48% vs 25% in helix; OR 2.50, Fisher p 0.041**).
Loop retrieval **AUC 0.700, 90% CI [0.49, 0.89], n=10 movers.**
**Honest limit:** the CI lower bound sits **at chance**. With 10 movers this is **not
distinguishable from chance**; the 0.70 point estimate is a small-sample artifact —
which Gate 4 then confirms. (QC also flagged/removed V111M and V111I here — both in
helix, so they did not touch the loop result.)

## Gate 4 — Powered loop replication (`powered_loop_gate.py`)

**Claim:** with adequate power on *independent* proteins, the loop signal does not hold.
**Numbers:** pooled barnase + human lysozyme + RNase A (SNase fell out — its engineered
variant background yielded only 1 WT crystal, so no noise floor). **n=87 loop pairs,
33 movers. AUC 0.591, 90% CI [0.49, 0.69], AP 0.497 (base 0.38), p@33 0.42.**
**Honest limit:** lower bound at chance → **not distinguishable from chance**; the pool
is human-lysozyme-dominated (27 of 33 movers). T4L's 0.70 regressed to 0.59 under
power — exactly the small-sample-overestimate reading.

## Gate 5 — Does 3D context rescue it? (`gate3_3d.py`)

This gate makes **two distinct claims that must not be blurred.**

**(1) As evidence for the mechanism — suggestive, direction & location correct.**
Adding the 3D contact environment moves mover-retrieval in the predicted direction and
in the predicted place: the **core subset** (packing-dominated, where local-only was
**0.477**) rises to **0.584**. That is exactly where tertiary context *should* matter.
It supports the diagnosis that the cause is tertiary.

**(2) As a validated predictive improvement — null.**
Overall lift is **0.516 → 0.590, 90% CI [0.52, 0.66]** (n=70 movers). The lower bound
barely clears 0.50; this is **not a certified improvement.** Absolute-bending RMSE
improved only 30.5° → 27.9° (9%) — the crude composition feature captures a thin slice
of 3D.

**These are different claims.** The 3D result confirms *why* local prediction fails
(the information is tertiary); it does **not** deliver a working predictor.

**Coincidence, not corroboration:** Gate 5's 0.59 and Gate 4's 0.59 are independent
tests (3D-on-T4L vs powered-loops) that happen to land on the same value. They do not
reinforce each other.

**Limit:** the 3D feature is contact-AA *composition only* — no distances, orientations,
cavity, or energy. A richer 3D model would mean entering established **structure-based
ΔΔG territory** (FoldX / Rosetta / ThermoMPNN), an occupied lane.

---

## Methodological notes worth highlighting

- **Two independent noise-floor estimates agree** (0.98° WT-crystal vs 0.75°
  same-variant) — the detection threshold is empirical, not assumed.
- **5.7× differencing error-cancellation** caught and quantified; the over-conservative
  √2×RMSE gate was corrected to a direct retrieval measurement.
- **Leakage control:** sequence-culled training (≤30% identity); validation families
  fully held out (and explicitly excluded by ID in Gates 4–5).
- **Bootstrap 90% CIs on every AUC**, because the mover counts are small.
- **Crystallographic-artifact QC:** conservative substitutions with implausibly large
  bends (>10°) flagged and removed (e.g. V111M −32.6°).
- **Thresholds are empirical** (per-window floors from redundant crystals), not a priori.

## Honest limitations (what even a complete version cannot claim)

- **Small mover counts throughout (10–72)** → wide CIs; nothing here is high-powered.
- **Single-family bias:** T4L for cores, a human-lysozyme-dominated pool for loops,
  SNase excluded. Not a proteome-wide result.
- **Crude 3D feature** (composition only); a fair test of "does tertiary context help,"
  not of how far full structure-based modeling could go.
- **29% above-floor is a likely ceiling** (mutation-tolerant, core-biased T4L; sparse-
  window floor inflation).
- **Single-sequence inputs** (no evolutionary profiles/MSAs); **Cα-only** bending metric;
  **observational** PDB mutants, not designed.

## Bottom line

> Mutation-induced backbone bending is real (Gate 1) but **not distinguishable-from-chance
> predictable from local sequence** — cores AUC 0.52, loops 0.59 even when powered, both
> with CIs touching 0.50. Adding crude 3D context nudges the core subset in the
> theoretically-predicted direction (0.48 → 0.58), confirming the **cause is tertiary**,
> but reaches no validated predictor. The signal lives in tertiary structure — the domain
> of heavy structure-based methods, not a light interpretable local model.
