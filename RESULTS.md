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
across runs. Every model and bootstrap uses a fixed seed, so this variation comes from
the **input set**: the original scripts ran live RCSB searches, and the PDB changes
between runs. The scripts now pin every search result to `manifests/*.json` on
first run, which removes this source of drift. That wobble is small — but it is *not* small
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

## Statistical revision (applies to Gates 2–5; numbers below not yet re-run)

Two defects in the original analysis were fixed in code after these numbers were produced:

1. **Per-pair bootstrap → residue-cluster bootstrap.** Validation pairs are not
   independent: the 72 T4L movers sit on only 38 residues, and the Gate 4 pool is
   human-lysozyme-dominated. The original CIs resampled pairs as if independent,
   which makes them **too narrow**. `stats_utils.cluster_auc_ci` now resamples whole
   residues (T4L) or (protein, residue) sites (Gate 4). Gate 4 also reports AUC per protein.
2. **Gate 5 lift was never tested as a difference.** The CI [0.52, 0.66] covers the 3D
   model's AUC alone, not the 0.516 → 0.590 **lift**. `stats_utils.paired_delta_auc`
   now scores both models on the same cluster resamples and reports a CI on ΔAUC,
   for both the full set and the core subset. Gate 2's AUC now also gets a CI.

Expected effect: CIs widen. That can only strengthen the local-only null. The Gate 5
lower bound of 0.52 may fall to chance. Until the gates are re-run, read every CI
below as **optimistic**.

## Gate 3 — Do loops rescue it? (T4L) (`loop_gate.py`)

**Claim:** loops *looked* like the one surviving home for local signal — but T4L
cannot certify it.
**Numbers:** movers enriched in loops (**48% of loop mutations are movers vs 25% in helix;
OR 2.50 for loop vs all non-loop, one-sided Fisher p 0.041**). Single protein, not replicated.
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
Overall lift is **0.516 → 0.590** (n=70 movers). The **90% CI [0.52, 0.66] covers only the
3D model's own AUC**, not the lift. The lift itself (paired ΔAUC) was never tested, and
the core-subset 0.477 → 0.584 has no CI at all. This is **not a certified improvement.** Absolute-bending RMSE
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

## v2 — clean labels + delta model (`pairs.py`, `delta_model.py`)

**Why.** Three problems found in the Gate 1–5 pipeline:
1. *Label construction.* There was no per-variant aggregation: T4L site 99 (L99A ligand soaks) supplied 61/248 rows, 4 of them movers, and dropping it raises the mover rate from 29% to 36%. There was no crystal-form matching, no resolution cut on the validation proteins and no ligand check (SPEC §2.3–2.4 were not implemented). The floor came from raw MADs of as few as 3 WT crystals, and the threshold ignored the mutant side's own noise. The label is fragile: 72 movers at 2σ, 41 at 3σ, 30 at |Δ|>3°.
2. *Engine.* Predict-absolute-then-difference has ~30° resolution against a ~2.7° effect, so its null says little about local sequence.
3. *Missing features.* On `feasibility_t4l.csv`, site identity has η² = 0.60 for |Δ| (≈0.28 expected by chance with 69 sites). A leave-one-out same-site diagnostic gets AUC 0.60, against the model's 0.52.

**What v2 does.** It keeps a mutant only when both sides are ≤2.5 Å. It uses a same-crystal-form WT reference and excludes mutant crystals whose ligand state near the window differs from the form's WT. It aggregates to one row per mutation. The floor is σ shrunk (k = 4 pseudo-crystals) toward a B-factor-conditioned prior. Within each crystal form, log σ is regressed on the window's WT B-factor (z within chain), with SE = 1.2533·σ·√(1/m + 1/n), and mover := |Δ|/SE > 2. The model is logistic regression or shallow gradient boosting on site and substitution features. Evaluation is leave-site-out (5-fold × 5 repeats) and leave-protein-out, with residue-cluster bootstrap CIs and a paired ΔAUC of site+subst over site.

**Tests.** Unit and synthetic end-to-end tests pass (`python -m pytest tests/`). They cover the parser against `parse_ca`, form matching, ligand exclusion, aggregation, the resolution cut, AUC equivalence, a planted-signal recovery and a null check.

**Run on real data.** Input sets are pinned in `manifests/`. Re-running `feasibility_t4l.py` on the pinned set reproduces Gate 1 exactly (248 windows, 72 above-floor events, 29%). Clean labels are in `pairs_clean.csv`, model output in `results/delta_model.json` and `delta_model.png`.

*What the filters remove* (`pairs.py` diagnostics):

| protein | structures | dropped: resolution | dropped: ligand mismatch (mutant crystals) | variant-forms with no clean crystal / no WT floor | clean mutations | movers |
|---|---|---|---|---|---|---|
| T4L | 618 | 19 | 112 | 42 / 17 | 113 | 36 |
| human lysozyme | 213 | 3 | 13 | 8 / 38 | 94 | 30 |
| RNase A | 323 | 8 | 10 | 10 / 6 | 13 | 5 |
| barnase | 51 | 7 | 0 | 0 / 5 | 15 | 4 |
| SNase | 290 | 8 | — | — / 92 | 0 | 0 (1 WT crystal) |

T4L falls from 248 window rows to 113 mutations. The drop comes mostly from ligand soaks (the L99A series) and from aggregating to one row per mutation. The T4L mover rate (32%) is close to Gate 1's 29%, so the phenomenon survives clean labelling. Pooled: 235 mutations, 75 movers (|z|>2), 121 sites, 4 proteins.

*Noise-floor prior.* The first real-data run shrank every window's σ toward one pooled σ. Window B-factor correlates with the raw WT σ (Spearman ρ = 0.51), so that pulled flexible windows' σ down, inflated their z, and could manufacture a B-factor → mover association. The default prior is now heteroscedastic (`pairs.PRIOR = "bfactor"`, `fit_sigma_prior`). The fitted slope is positive in every crystal form: 0.30 for T4L and 0.09–0.59 elsewhere, in log σ per 1 SD of B. `--prior pooled` reproduces the old labels. Mover counts barely move: 4 of 235 labels flip at z>2 (two high-B movers lost, two low-B gained), and 6 at z>3.

*Delta model* (B-factor prior, 90% residue-cluster bootstrap CIs):

| features | model | leave-site-out AUC | leave-protein-out AUC |
|---|---|---|---|
| site | logreg | 0.641 [0.55, 0.72] | 0.650 [0.56, 0.72] |
| subst | logreg | 0.522 [0.45, 0.59] | 0.573 [0.50, 0.64] |
| site+subst | logreg | 0.619 [0.53, 0.70] | 0.625 [0.54, 0.70] |
| site | hgb | 0.529 [0.45, 0.60] | 0.519 [0.44, 0.59] |
| subst | hgb | 0.511 [0.44, 0.59] | 0.608 [0.55, 0.67] |
| site+subst | hgb | 0.562 [0.49, 0.63] | 0.565 [0.49, 0.63] |

Paired ΔAUC, site+subst − site (logreg): −0.022 [−0.06, +0.02] site-out and −0.025 [−0.07, +0.02] protein-out. With the pooled prior, site-only hgb scored below chance (0.42) and produced spurious "significant" hgb contrasts. Under the B-factor prior those contrasts are +0.03 and +0.05, with CIs spanning 0. Gradient boosting is still not usable at n = 235, so the tables report logreg.

*Robustness checks* (leave-site-out logreg; scratch analysis, not a committed script):

| noise-floor prior | mover label | movers | site AUC | ΔAUC site+subst − site |
|---|---|---|---|---|
| **B-factor (default)** | z>2 | 75 | 0.645 [0.55, 0.72] | −0.022 [−0.06, +0.02] |
| B-factor | z>3 | 38 | 0.704 [0.63, 0.78] | −0.023 [−0.07, +0.02] |
| B-factor + WT bend | z>2 | 74 | 0.601 [0.50, 0.69] | −0.027 [−0.07, +0.02] |
| B-factor + WT bend | z>3 | 39 | 0.679 [0.60, 0.76] | −0.005 [−0.05, +0.04] |
| pooled (first run) | z>2 | 75 | 0.620 [0.53, 0.70] | +0.001 [−0.04, +0.05] |
| pooled | z>3 | 40 | 0.727 [0.65, 0.79] | −0.009 [−0.05, +0.03] |
| none (raw σ) | z>2 | 81 | 0.547 [0.45, 0.63] | +0.002 [−0.05, +0.06] |
| none (raw σ) | z>3 | 46 | 0.59 [0.51, 0.66] | −0.005 [−0.06, +0.05] |

Within-protein site AUCs under the B-factor prior (z>2) are T4L 0.61, human lysozyme 0.70 and RNase A 0.71. The pooled number is therefore not just between-protein base rates, and leave-protein-out agrees with leave-site-out. The strongest site features are WT window bend (standardized coefficient −0.67; straighter windows move more), window B-factor (+0.38) and non-local contacts (+0.35). The out-of-fold site score also ranks the threshold-free |Δ| in degrees: Spearman 0.30 [0.14, 0.43].

**Reading.**
1. *Substitution identity adds nothing once the site is known.* ΔAUC is ≤ +0.01 with CIs spanning 0 under all eight label definitions. This is the local-sequence question asked with an engine that predicts Δ directly, so it closes the "the engine couldn't see it" loophole in Gates 2–5.
2. *Site context carries a modest signal, and the B-factor confound does not explain it.* Conditioning the prior on B removes movers where B is high, yet site AUC rises slightly (0.62 → 0.64). One residual caveat: WT bend, the strongest feature, also correlates with raw σ (ρ = −0.44). Adding it as a second prior covariate lowers site AUC to 0.60 [0.50, 0.69] at z>2, with the CI touching chance, and to 0.68 [0.60, 0.76] at z>3. Straight windows are both noisier and more mutation-sensitive, and with 4 proteins these cannot be fully separated. **Do not read 0.60 as a lower bound on the site effect.** It corrects only for the two noise covariates we modelled (B-factor and WT bend). Any unmodelled covariate of WT noise that also correlates with the site features would lower it further. Treat it as the upper end of what survives correction so far. Raw σ without shrinkage gives the weakest signal, but it is also the noisiest label (few-crystal MADs), so it does not settle the question.
*Is the WT-bend effect just secondary structure?* No (B-factor prior, z>2; scratch analysis). SS explains about half the variance of WT bend (η² = 0.49), but it does not predict movers:

| check | result |
|---|---|
| mover rate by SS | H 0.30 (n = 141), E 0.35 (n = 34), L 0.33 (n = 60) |
| SS-only model | AUC 0.35 [0.28, 0.44] (no signal; below 0.5 is a grouped-CV artifact) |
| WT bend, residualised on SS | AUC 0.64 [0.54, 0.71] |
| WT bend within helices only | AUC 0.61 [0.50, 0.70] |
| ΔAUC, SS + WT bend − SS only | +0.26 [+0.11, +0.38] |
| ΔAUC, site − (site without WT bend) | +0.10 [+0.03, +0.16] |

The bend signal lives *within* SS classes, so it is not a proxy for SS. What remains open is whether within-class bend marks windows that are noisier in WT crystals or windows that really move on mutation.

3. The honest v2 claim: **"where" weakly predicts movers (AUC ~0.60–0.70 depending on the label), "what" adds nothing.** More proteins with ≥3 WT crystals per form would be needed to separate flexibility-driven noise from flexibility-driven movement.

## Methodological notes worth highlighting

- **Two independent noise-floor estimates agree** (0.98° WT-crystal vs 0.75°
  same-variant) — the detection threshold is empirical, not assumed.
- **5.7× differencing error-cancellation** caught and quantified; the over-conservative
  √2×RMSE gate was corrected to a direct retrieval measurement.
- **Leakage control:** sequence-culled training (≤30% identity); validation families
  fully held out (and explicitly excluded by ID in Gates 4–5).
- **Bootstrap 90% CIs on AUCs**, because the mover counts are small. Originally these
  were per-pair (optimistic) and missing for Gate 2 and the core subset. Both are now
  fixed in code (residue-cluster bootstrap plus paired ΔAUC); see *Statistical revision*.
- **Crystallographic-artifact QC:** conservative substitutions with implausibly large
  bends (>10°) flagged and removed (e.g. V111M −32.6°).
- **Thresholds are empirical** (per-window floors from redundant crystals), not a priori.
  "Mover" = |Δ| > 2× the per-window robust σ, so roughly 5% of non-movers would
  cross it by noise alone. That is the right baseline against which to read the 29%.

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
> theoretically-predicted direction (0.48 → 0.58), **consistent with** a tertiary cause
> (the lift itself is not yet statistically tested), and reaches no validated predictor. The signal lives in tertiary structure — the domain
> of heavy structure-based methods, not a light interpretable local model.
>
> **v2 update (clean labels, delta-targeted model):** substitution identity adds nothing
> over site context (paired ΔAUC −0.02 [−0.06, +0.02]); site context alone reaches
> AUC 0.64 [0.55, 0.72] with a B-factor-conditioned noise floor. That signal survives
> removing the B-factor confound, but only barely clears chance once WT bend also
> enters the prior (0.60 [0.50, 0.69]).
