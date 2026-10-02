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

## v3 — systematic miner (`mine_pairs.py`): 191 proteins, 158 families

> The v3 labels used an uncorrected MAD noise floor. The v4 audit below found that it
> made the z test liberal in forms with few WT crystals, and v4 supersedes the v3
> numbers. The miner and the family definition are unchanged.

**Why.** v2 had 4 usable proteins. That was too few to tell apart two explanations of the site signal: windows that are noisy in WT, or windows that really move on mutation. SPEC §2.1's primary source (RCSB + SIFTS) removes that limit.

**How.**
- *Candidates:* an RCSB facet finds 588 UniProt accessions with ≥ 8 X-ray (≤ 2.5 Å) entities carrying exactly one SIFTS substitution.
- *Pre-screen (metadata only):* keep crystal forms with ≥ 3 WT crystals of one reference construct plus ≥ 1 single mutant, in single-protein-entity entries only. 283 accessions pass. 278 are new; the other 5 are the gate proteins.
- *Download:* 6120 of 6449 entries have PDB-format files.
- *Families:* RCSB sequence search at ≥ 30% identity. A hit links two proteins only if its sequence *is* the other protein's reference construct. Linking through UniProt annotations merged 86 unrelated proteins through fusion constructs (T4L–GPCR, MBP-, GST-, ubiquitin-tagged).

Everything is pinned in `manifests/mined.json`. `pairs.py` applies the same filters as in v2 (B-factor prior, z>2):

| structures | dropped: resolution | dropped: ligand mismatch (mutant crystals) | variant-forms with no clean crystal / no WT floor | clean mutations | movers | proteins | families |
|---|---|---|---|---|---|---|---|
| 7612 | 45 | 1243 | 636 / 274 | 892 | 335 (37.6%) | 191 | 158 |

**Delta model** (90% residue-cluster bootstrap CIs; family-out = 10 folds of whole ≥ 30%-identity families):

| features | model | leave-site-out AUC | leave-family-out AUC |
|---|---|---|---|
| SS only | logreg | 0.473 [0.43, 0.51] | 0.509 [0.47, 0.55] |
| site | logreg | 0.587 [0.55, 0.62] | 0.595 [0.56, 0.63] |
| subst | logreg | 0.564 [0.53, 0.60] | 0.541 [0.50, 0.58] |
| site+subst | logreg | 0.604 [0.57, 0.64] | 0.593 [0.56, 0.63] |
| site | hgb | 0.551 [0.51, 0.59] | 0.546 [0.51, 0.58] |
| site+subst | hgb | 0.559 [0.52, 0.59] | 0.567 [0.53, 0.60] |

| paired contrast (logreg) | leave-site-out | leave-family-out |
|---|---|---|
| site − SS | +0.114 [+0.07, +0.16] | +0.085 [+0.05, +0.12] |
| site+subst − site | +0.017 [−0.01, +0.04] | −0.002 [−0.03, +0.02] |

**Null control** (`pairs.null_rows`). Up to 10 WT crystals per form are held out one at a time. Each is scored as a single-crystal "mutant" against statistics recomputed without it, through the same `label()` code. This gives 2699 pseudo-pairs at 456 sites.
- *z calibration:* 6.9% of pseudo-pairs exceed |z| > 2 (nominal ≈ 4.6%), against 37.6% of real mutations. The z threshold is close to calibrated, and mutation-induced movement is real in aggregate.
- *Where the noise is:* site features partly predict *where* noise exceedances happen. The site model trained on pseudo labels reaches AUC 0.564 [0.52, 0.61] (site-out) and 0.560 [0.51, 0.61] (family-out).

**Decomposition** (scratch analysis, 5 site folds shared between real and null rows):

| score, evaluated on real labels | AUC |
|---|---|
| noise score (site model trained on null labels, other sites) | 0.555 [0.52, 0.59] |
| site model trained on real labels | 0.593 [0.56, 0.63] |
| ΔAUC, (noise score + site features) − noise score | +0.047 [+0.01, +0.08], P(Δ ≤ 0) = 0.02 |
| real-trained site model scored on *null* labels | 0.530 [0.49, 0.57] |

The coefficients differ between the two fits. Non-local contacts carry movement: +0.21 on real labels, −0.14 on null. Window B-factor also flips sign: +0.14 real, −0.12 null. WT bend is negative in both (−0.44 real, −0.24 null), so it is partly a noise marker.

**Reading.**
1. *"What" has no edge.* With 892 mutations and 158 families, substitution identity adds −0.002 [−0.03, +0.02] over the site. This is the project's central question, now answered at 4× the v2 sample and with family hold-out.
2. *"Where" has a small, real edge.* Site features beat SS by +0.085 (family-out), and SS alone is at chance. A noise-only score reaches AUC 0.555 on real labels, which is more than half of the site model's margin above chance (0.593). The movement-specific remainder is ΔAUC ≈ +0.05 [+0.01, +0.08], carried mainly by non-local contacts, consistent with the tertiary-packing reading of Gate 5.
3. Neither is a usable predictor: the best AUC is ≈ 0.60. The v2 single-series numbers (0.64, or 0.60 after WT-bend correction) were within noise of this pooled estimate.

## v4 — audit: label calibration, reliability ceiling, richer features

**Why.** Before reaching for a better model, check two things. Are the labels calibrated? And how much accuracy can these labels support at all?

### Issues found and fixed

| issue | evidence | fix |
|---|---|---|
| Noise floor biased low at small n | σ = 1.4826 × MAD has no finite-sample correction. That gives 0.67σ on average at n = 3 and 0.54σ at the median. On the WT-vs-WT null, 17–20% of pseudo-mutants were "movers" in forms with 3–4 WT crystals (overall 6.9%; nominal 4.6%). | σ = SD / c4(n), which is unbiased. Corrected MAD and Qn remain as options (`SIGMA_EST`). |
| Robust scale vs mixed WT states | With the bias corrected, MAD and Qn still give 5.4–7.3% and 5.9% null false positives, against 3.1–4.0% for SD. The null z has sd 1.4–1.7 but robust sd 0.65–0.8 (heavy tails), a mixture: most WT crystals of a form agree and a few sit in another conformation. A robust scale ignores that minority, and a mutant crystal in the minority state then reads as a mover. | SD by default. Its errors are conservative. |
| SE used the asymptotic median efficiency for every m | 1.2533 overstates the SE of a single crystal (m = 1, 85% of rows). Large-n forms ran at 2% false positives. | Exact small-sample efficiency table, `med_eff(k)`. |
| σ itself is uncertain at small n | Even unbiased, σ from 3–9 crystals leaves z with t-like tails: 5.9–7.8% null false positives in those bins. | δ/SE is mapped through a Student-t with n − 1 + K degrees of freedom onto the normal scale (`Z_CALIB = "t"`). |
| Ligand rule applied only to mutants | WT reference crystals with a ligand next to the window entered the median and σ. | Rule applied to WT too when ≥ 3 clean crystals remain (`WT_LIG_FILTER = "prefer"`). The scaffold is taken from the typical ligand state. |
| Additives counted as ligands | Of 1268 dropped mutant crystals, about 190 failed only on SO4, Cl⁻, glycerol, BME, EDO and similar. Real ligands (CMO/NO in myoglobin, UMP, benzene soaks) are a different, legitimate group. | Common additives are ignored by the rule; metals are not (`ADDITIVES`). The 884 null pseudo-pairs this adds have *fewer* false positives (3.1% vs 3.9%), so additives do not create noise. |

Checked and fine:
- *Within-protein:* the signal is site-level, not between-protein base rates. Within-protein AUC is higher than pooled (0.65 vs 0.60).
- *Bootstrap level:* CIs are the same whether sites, proteins or families are resampled.
- *Fold leakage:* buffering ±4 residues between train and test sites changes leave-site-out AUC by 0.005, so overlapping windows do not leak.

**Calibration after the fixes** (`diagnostics.py`, 3573 WT-vs-WT pseudo-pairs):

| WT crystals behind the floor | 2–4 | 5–9 | 10–19 | 20+ | all |
|---|---|---|---|---|---|
| null false positives, first v3 labels | 17.1% | 12.2% | 4.5% | 2.1% | 6.9% |
| null false positives, v4 labels | 3.5% | 4.9% | 1.3% | 2.0% | **2.4%** |
| real mutations called movers, v4 | 17.9% | 33.0% | 32.6% | 27.3% | 26.5% |

### v4 labels

| structures | dropped: resolution | dropped: ligand mismatch (mutant crystals) | variant-forms with no clean crystal / no WT floor | clean mutations | movers | sites | proteins | families |
|---|---|---|---|---|---|---|---|---|
| 7612 | 45 | 1006 | 480 / 271 | 1045 | 277 (26.5%) | 712 | 207 | 166 |

By secondary structure: helix 23.9% (126/527), strand 34.3% (69/201), loop 25.9% (82/317), χ² p = 0.017. Loops are again not enriched (loop vs rest p = 0.82). 85% of the labels come from a single mutant crystal.

### How much accuracy can these labels support? (split-half test-retest)

Each protein's crystals are split at random into two halves, and labels are rebuilt on each half independently. Each seed yields about 70 mutations measured in both halves (3 seeds).

| corr(z_A, z_B) | all-crystal reliability (Spearman–Brown) | mover agreement | κ | P(mover in B \| mover in A) | direction agreement among shared movers | oracle AUC |
|---|---|---|---|---|---|---|
| 0.75 | ≈ 0.86 | 72% | 0.37 | 62% | 100% | **0.76** |

"Oracle AUC" is the AUC of |z_B| for predicting mover_A: what a predictor that knew each effect as well as half the data could reach. Each half has one mutant crystal and half the WT floor, which is close to a typical single-crystal label. **Against labels like these, accuracy saturates near AUC 0.75–0.8 whatever the model.** Real movers reproduce (same direction every time), but whether a borderline effect crosses |z| > 2 is close to a coin flip.

### Delta model on v4 labels

90% residue-cluster bootstrap CIs; family-out = 10 folds of whole families; within-protein = only mover/non-mover pairs from the same protein.

| features | model | leave-site-out AUC | leave-family-out AUC | within-protein (family-out) |
|---|---|---|---|---|
| SS only | logreg | 0.504 [0.46, 0.55] | 0.523 [0.49, 0.56] | 0.516 |
| site (C-alpha) | logreg | 0.604 [0.57, 0.64] | 0.596 [0.56, 0.63] | 0.648 |
| subst | logreg | 0.572 [0.54, 0.61] | 0.570 [0.53, 0.60] | 0.578 |
| site + subst | logreg | 0.620 [0.58, 0.65] | 0.610 [0.57, 0.65] | 0.653 |
| where (site + full-atom context) | logreg | 0.583 [0.54, 0.62] | 0.574 [0.53, 0.61] | 0.610 |
| where + what (+ subst + interactions) | logreg | 0.609 [0.57, 0.64] | 0.606 [0.57, 0.64] | 0.593 |
| where + what | hgb | 0.613 [0.58, 0.64] | 0.597 [0.56, 0.63] | 0.617 |
| site + ESM-2 site terms | logreg | 0.597 [0.56, 0.63] | 0.594 [0.56, 0.63] | 0.648 |
| site + ESM-2 site + substitution LLR | logreg | 0.602 [0.56, 0.64] | 0.597 [0.56, 0.63] | 0.630 |

| paired contrast (logreg) | leave-site-out | leave-family-out |
|---|---|---|
| site − SS | +0.099 [+0.06, +0.14] | +0.073 [+0.03, +0.11] |
| (site + subst) − site | +0.017 [−0.00, +0.04] | +0.014 [−0.01, +0.04] |
| where − site (full-atom context) | −0.021 [−0.04, −0.00] | −0.023 [−0.04, −0.00] |
| (where + what) − where | +0.026 [+0.00, +0.05] | +0.032 [+0.01, +0.06] |
| ESM site terms − nothing (over site) | −0.006 [−0.01, +0.00] | −0.002 [−0.01, +0.00] |
| ESM substitution LLR − ESM site terms | +0.005 [−0.01, +0.02] | +0.003 [−0.01, +0.01] |

**Null control and decomposition.** The site model trained on the pseudo labels reaches 0.572 leave-site-out and 0.529 family-out on them. As a noise score on the *real* labels it reaches only **0.518 [0.48, 0.55]**. Beyond that noise score, site features add **+0.087 [+0.04, +0.13]** (P = 0.001; v3 found +0.047). Almost all of the site signal is now movement, not noise structure. Standardized coefficients, real vs null labels:

| feature | real movers | null |
|---|---|---|
| WT window bend | −0.47 | −0.05 |
| burial direction (hse_up) | +0.21 | +0.06 |
| window B-factor | +0.22 | −0.03 |
| non-local contacts | +0.12 | +0.02 |
| distance to terminus (log) | +0.03 | +0.32 |

WT bend was partly a noise marker in v3 (−0.24 on the null). With calibrated labels it is movement-specific, which settles the question v2 left open.

### What buys accuracy, what does not (`model_variants.py`, leave-family-out)

| change | AUC | Δ vs logreg baseline |
|---|---|---|
| site, logreg (baseline) | 0.597 | — |
| site, C tuned by nested CV | 0.597 | −0.000 [−0.003, +0.003] |
| site, ridge on continuous \|z\| | 0.606 | +0.009 [−0.001, +0.020] |
| site, regularized boosting | 0.583 | −0.014 [−0.047, +0.020] |
| where + what, logreg (baseline) | 0.607 | — |
| where + what, ridge on continuous \|z\| | 0.626 | +0.019 [+0.004, +0.034] |
| where + what, regularized boosting | 0.629 (within-protein 0.691) | +0.021 [−0.011, +0.054] |
| where + what, ensemble logreg + boosting | **0.630** [0.60, 0.66] | +0.023 [+0.006, +0.040] |

Learning curve (site, logreg): with 25 / 50 / 75 / 100% of training families, family-out AUC is 0.546 / 0.580 / 0.592 / 0.597 (within-protein 0.576 / 0.616 / 0.628 / 0.648). Data still helps but is flattening, about +0.005–0.01 per extra quarter.

**Reading.**
1. *The labels were the weakest link, and fixing them changed the interpretation more than the AUC.* Pooled AUC barely moved (0.595 → 0.596 for site). But the part of it that is noise fell from 0.555 to 0.518, and the movement-specific part nearly doubled (+0.047 → +0.087).
2. *"What" carries a small amount of information, but only in context.* Context-free substitution descriptors add +0.014 (n.s.). Substitution × local-structure terms add +0.03 over the full-atom "where", but that "where" is itself 0.02 below the simple site model, so the net gain over site is about +0.01. ESM-2, which scores the substitution from whole-sequence evolutionary context, adds nothing (+0.003). Whatever decides which substitutions bend the backbone is not in sequence-level substitution plausibility.
3. *Richer site description does not help a linear model.* Full-atom context (φ/ψ, H-bonds, packing, ligand distance, SS position) is largely redundant with the C-alpha site features. Only ligand distance helps on its own (+0.011). The rare physical mechanisms (Gly at positive φ, Pro in a helix) occur in 1% or less of the deposited mutations.
4. *Accuracy ceiling.* Best model ≈ 0.63 family-out (boosting alone: 0.69 within-protein) against an oracle ≈ 0.76 on single-crystal-like labels. Models now capture roughly half of the achievable margin above chance. Two levers remain for the other half: better labels (replicate mutant crystals; 85% of labels rest on one crystal) and more proteins. The learning curve says the second lever is weak.
5. The best variant was picked from about 12 (6 learners × 2 feature sets), so 0.630 is slightly optimistic. The pre-specified logreg numbers are the ones to quote.

## v5 — confidence-weighted training, expanded set, final model

**What changed.**
- *Training weights:* each training row is weighted by label confidence, |z| distance from the threshold (`delta_model.label_weights`). Split-half showed that borderline labels flip between independent crystal sets.
- *Second accuracy metric:* AUC on **confident labels** (|z| > 3 vs |z| < 1), i.e. how well the model separates mutations whose label is not in doubt.
- *Expanded miner set:* the miner threshold is ≥ 4 single-substitution entities (`manifests/mined_ms4.json`; v3/v4 used ≥ 8). That gives 447 mined proteins in 322 families.
- *Family caveat:* one 36-protein family (F003) chains unrelated folds through single-linkage. This only makes family hold-out stricter.

**Data.** 1250 clean mutations, 332 movers (26.6%), 895 sites, 319 proteins, 250 families. WT-vs-WT null: 2.6% false positives (n = 4071), ≤ 4.5% in every n_wt bin. Split-half: corr(z) 0.75, κ 0.35, direction agreement 100%, **oracle AUC 0.75**.

**Delta model** (leave-family-out, confidence-weighted, 90% residue-cluster CIs):

| features | model | AUC, all labels | AUC, confident labels | within-protein |
|---|---|---|---|---|
| SS only | logreg | 0.539 [0.50, 0.57] | 0.550 | 0.527 |
| site | logreg | 0.625 [0.59, 0.66] | 0.684 [0.64, 0.73] | 0.646 |
| site + subst | logreg | **0.639** [0.61, 0.67] | **0.714** [0.67, 0.76] | 0.634 |
| where (site + full-atom context) | hgb | 0.624 [0.59, 0.65] | 0.719 [0.68, 0.76] | 0.654 |
| where + what | logreg | 0.633 [0.60, 0.66] | 0.706 [0.66, 0.75] | 0.593 |
| where + what | hgb | 0.631 [0.60, 0.66] | 0.724 [0.68, 0.76] | 0.686 |
| site + ESM-2 | logreg | 0.625 [0.59, 0.66] | 0.673 | 0.628 |
| ensemble logreg + hgb, where + what (`model_variants.py`) | — | 0.649 [0.62, 0.68] | — | 0.680 |

| paired contrast (logreg, family-out) | ΔAUC |
|---|---|
| site − SS | +0.086 [+0.05, +0.12] |
| (site + subst) − site | +0.014 [−0.00, +0.03] |
| (where + what) − where | +0.018 [−0.00, +0.04] |
| ESM-2 substitution LLR − ESM-2 site terms | +0.002 [−0.01, +0.01] |

Noise decomposition (site): the noise score reaches 0.556 on real labels; site features add **+0.072 [+0.03, +0.11]** beyond it. Learning curve (site, 25 / 50 / 75 / 100% of training families): 0.596 / 0.611 / 0.618 / 0.620. **It is saturated:** 60% more proteins than v4 lifted AUC by about 0.03, and further families now add almost nothing.

**Final model** (`predict.py train`): site + subst logistic regression, confidence-weighted, with Platt calibration on out-of-fold scores. Leave-family-out AUC **0.639 [0.61, 0.67]**, **0.716 [0.67, 0.76]** on confident labels, within-protein 0.635. `predict.py score --pdb ID --chain A --mut L99A,...` ranks mutations on any WT structure.

**Reading (final).**
1. **Mutation-induced local backbone movement is real and calibrated.** 26.6% of clean single mutations move a 5-residue window beyond crystal noise, against 2.6% of WT-vs-WT pseudo-mutants.
2. **"Where" is the information.** Site context (straight WT window, burial, non-local contacts) separates movers from non-movers at AUC ≈ 0.63 on all labels and ≈ 0.70 on confident ones, beyond SS (+0.09), and most of it is movement rather than noise.
3. **"What" (the local sequence change) adds ≈ 0.01–0.02 AUC once "where" is known**, whether it is described physically or by an evolutionary language model. The original hypothesis, that local sequence drives which mutations bend the backbone, is rejected at this resolution.
4. **Ceiling.** Labels support at most AUC ≈ 0.75. The final model reaches 0.64, about 55% of the achievable margin above chance. The learning curve is flat and every model variant is within ±0.02, so the remaining gap is label noise (85% single-crystal labels), not model capacity or data volume.

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
- **v4:** 85% of the clean labels rest on a single mutant crystal, and split-half
  reliability puts the achievable AUC near 0.76. The bending metric sees only a
  5-residue Cα angle; twists and shifts that keep that angle are invisible to it.

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
>
> **v3 update (191 proteins, 158 families, WT-vs-WT null control):** "what" adds
> nothing (−0.002 [−0.03, +0.02]). "Where" predicts movers at AUC 0.60 [0.56, 0.63]
> (family hold-out), beyond SS. More than half of that margin is WT-noise structure; the
> movement-specific part is ΔAUC ≈ +0.05 [+0.01, +0.08], driven by non-local contacts.
>
> **v4 update (audited labels: 2.4% null false positives; 1045 mutations, 207 proteins):**
> the old floor was liberal in few-crystal forms; corrected, the site signal is almost
> all movement (excess over a noise score +0.087 [+0.04, +0.13]). Substitution identity
> adds at most ~0.01–0.03, and only as substitution × structure terms; ESM-2 adds nothing.
> Best model AUC ≈ 0.63 (family hold-out) against a label-reliability ceiling ≈ 0.76.
>
> **v5 (final; 1250 mutations, 319 proteins, 250 families, confidence-weighted):** site +
> substitution logistic regression reaches AUC 0.64 [0.61, 0.67] with whole families held
> out, and 0.71 on confident labels; the substitution's own share is +0.014. The learning
> curve is saturated and the label ceiling is ≈ 0.75. `predict.py` ships the model.
