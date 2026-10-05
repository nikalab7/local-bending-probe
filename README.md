# local-bending-probe

## How much information about protein backbone bending is contained in local sequence?

Modern protein-structure models rely heavily on long-range interactions and evolutionary information. This project asks a simpler question:

**If we deliberately remove all non-local information, how much can local amino-acid sequence alone tell us about mutation-induced backbone bending?**

To answer that, I built a lightweight, interpretable pipeline designed around a single idea:

> Hold local geometry constant, then test whether local sequence patterns can explain which mutations bend the backbone.

The original expectation was that specific sequence motifs would emerge as reliable local drivers of bending. Instead, the project arrived at the opposite conclusion.

On 1,250 clean WT/mutant pairs from 319 proteins, the identity of the substitution adds little or nothing to predicting which mutations move the backbone. *Where* the mutation sits carries a modest signal on the development data: straight, buried windows with non-local contacts move more. **A pre-registered lockbox of 90 unseen proteins did not confirm that signal** (AUC 0.53 [0.40, 0.66]). It was too small to detect an effect of the size seen in development, so every result here is exploratory.

The failure of the local model became the result.

---

## Final result (pre-registered, `PROTOCOL.md`)

The label, features, model, metrics and analysis plan were frozen in `PROTOCOL.md` (commit `fe9c3c7`) before any final evaluation. The lockbox was then scored once (`results/final_lockbox.json`). It holds 90 proteins in 89 families that never entered any decision, none with >= 30% sequence identity to a development protein. All CIs below are 95% family-bootstrap intervals.

| | **lockbox** (n = 111, 30 movers, 89 families) | dev (n = 1250, 447 movers, 250 families; exploratory) |
|---|---|---|
| **model AUC** (site + substitution, logistic regression) | **0.532 [0.403, 0.661]** | 0.662 [0.633, 0.698] |
| label ceiling (split-half oracle, dev) | 0.787 | 0.787 |
| SS-only baseline | 0.530; model − SS **+0.002 [−0.099, +0.117]**, p = 0.98 | 0.576; +0.087 [+0.052, +0.120], p < 0.001 |
| burial-only baseline | 0.449; model − burial **+0.083 [−0.038, +0.204]**, p = 0.18 | 0.550; +0.112 [+0.076, +0.160], p < 0.001 |
| within-protein AUC | 0.375 [0.00, 0.83] (almost no within-protein pairs) | 0.653 [0.635, 0.711] |
| without T4L | 0.538 (trained without T4L) | 0.666 (retrained without T4L) |
| positive control (predict helix) | 0.951 | 0.915 |
| negative control (labels shuffled within families) | 0.538 ± 0.018 (degenerate: 73 of 89 families have one row) | 0.529 ± 0.016 |
| minimum detectable AUC (alpha 0.05, power 0.8) | **0.685** (n_eff 72.9 families) | 0.545 (n_eff 30.4) |
| WT-vs-WT null false-positive rate | 4.7% (median 3 WT crystals per form) | 2.6% (median 7) |
| excluding artifact-suspect rows ("any") | 0.583 [0.382, 0.791] (46 excluded) | 0.701 [0.646, 0.759] (631 excluded) |

**What this means.**
1. **The pre-registered success criterion was not met.** The lockbox AUC CI includes 0.5, and the model does not beat either baseline on the lockbox.
2. **The lockbox is inconclusive rather than negative.** It could only have detected an AUC of about 0.685 or more, which is above the dev estimate of 0.66, and its CI contains 0.66 as well as 0.5. The size rule in `PROTOCOL.md` (>= 30 rows, >= 10 per class, >= 10 families) was met, but it was too lax: the power analysis shows the lockbox is too small for this effect. **All results are therefore exploratory.**
3. **The lockbox labels are noisier.** Lockbox proteins have few WT crystals (median 3 vs 7), so the WT-vs-WT null false-positive rate is 4.7% against 2.6% in dev. About 5 of the 30 lockbox "movers" are expected to be noise.
4. **The dev number contains a between-protein component.** Shuffling labels within families still gives AUC 0.53 on dev, because the model partly learns which families have more movers. The within-protein AUC (0.653) is the cleaner dev estimate of the site-level signal.
5. **Artifacts do not explain the dev signal.** Removing rows with temperature or resolution mismatches, mutant altlocs or lattice contacts at the mutated residue raises dev AUC (0.70). Those cases are noisier, not the source of the signal. They were never model inputs.
6. **"What" vs "where" is unaffected.** The substitution adds about +0.01 on dev, and no added feature (ESM-2, elastic network, full-atom context) helped (`WORKLOG.md`, decision log).

A confirmatory test needs a larger lockbox. The minimum detectable AUC scales as 0.5 + 2.8 x SE, with SE proportional to 1/sqrt(families). Detecting AUC 0.66 needs about 1.4x the current lockbox; detecting AUC 0.60 (plausible with noisier labels), or the +0.09 margin over SS-only, needs about 3–4x. That means 300+ independent families, ideally with more WT crystals per form, e.g. from future PDB depositions.

---

## The question

Protein structure is often discussed as a sequence-to-structure problem, but that framing hides an important distinction.

A protein's backbone can bend for many reasons:

* local amino-acid preferences,
* secondary-structure tendencies,
* packing interactions,
* long-range contacts,
* solvent effects,
* global folding constraints.

The goal of this project was to isolate the first factor.

For a given mutation, two things can be known: **what** it is (the substitution, i.e. the local-sequence change) and **where** it is (the site's structural context). Can "what" predict whether the backbone moves once "where" is known?

If the answer were yes, it would suggest that interpretable local rules explain a meaningful fraction of backbone deformation.

If the answer were no, it would imply that the information lives elsewhere.

---

## Approach

The project was designed as a sequence of falsification tests rather than a search for positive results.

The workflow was:

1. Build a robust bending metric.
2. Verify that mutation-induced bending exists in real structures.
3. Build clean WT/mutant labels: same crystal form, ≤ 2.5 Å, matched ligand state, one row per mutation, and a noise floor that depends on flexibility.
4. Mine the PDB systematically for every protein that can supply such pairs.
5. Predict movers directly from site and substitution features, holding out whole sequence families.
6. Run a null control on WT-vs-WT pseudo-mutants to separate real movement from noise.
7. Audit the labels themselves: check the false-positive rate on the null, and measure test-retest reliability to know how much accuracy is achievable at all.

Every stage had a predefined failure condition.

The objective was not to maximize performance but to determine where the predictive information actually resides.

The first version of the pipeline (Gates 1–5, T4 lysozyme plus four validation proteins) is kept for the record. An audit showed its labels and its prediction engine were too weak to support a negative conclusion. A second audit (v4) found the rebuilt labels' noise floor too liberal when few WT crystals are available, and fixed it. The results below come from the current pipeline (v5: calibrated labels, confidence-weighted training, 319 proteins; details further down and in `RESULTS.md`).

---

## What worked

The phenomenon itself is real.

Across 319 proteins, 26.6% of clean single mutations (332/1250) moved the backbone beyond |z| > 2 of the per-window noise floor. Held-out WT crystals scored the same way as pseudo-mutants cross that threshold only 2.6% of the time, and no more than 4.5% in any group of crystal forms. Mutations do move protein backbones.

The first T4 lysozyme run found 29% (72/248). That run also reported more movers in loops than in helices (48% vs 25%, one-sided Fisher p = 0.041). **That enrichment does not replicate on the clean pooled labels:** loops 28.5% (113/397), all other residues 25.7% (219/853), p = 0.33. A secondary-structure-only model is near chance (AUC 0.54 [0.50, 0.57]).

---

## What failed

The central hypothesis did not survive.

The substitution, the local-sequence part of the question, adds little once the site is known. It was tested three ways (paired ΔAUC, leave-family-out):

* Context-free substitution features (Δvolume, Δhydrophobicity, charge, Pro/Gly, BLOSUM62, cavity × burial): **+0.014 [−0.00, +0.03]**
* Substitution × local-structure terms (Gly at positive φ, Pro strain, overpacking, lost side-chain H-bonds, helix/strand propensity change): +0.018 [−0.00, +0.04] over a full-atom site description
* A protein language model (ESM-2) scoring the substitution in its whole-sequence evolutionary context: **+0.002 [−0.01, +0.01]**

This time the null is informative. The first engine predicted absolute bending (~30° error) and differenced two predictions, so it could not have seen a ~3° effect even if the information were there. The new model predicts movers directly and does find signal in site features. Neither the substitution's physics nor its evolutionary plausibility adds much to that.

---

## The most informative result

Where the mutation sits carries information **on the development data** (v5 numbers below; the frozen v6 pipeline gives 0.662, see *Final result*). The pre-registered lockbox did not confirm it (AUC 0.53 [0.40, 0.66], underpowered), so treat this section as exploratory:

* Site features: AUC **0.63 [0.59, 0.66]** with whole families held out, **0.68** on labels that are clearly mover or clearly not
* Beyond secondary structure: **+0.086 [+0.05, +0.12]**

The null control shows that this is mostly movement, not noise. A site model trained only on WT-vs-WT noise scores the real labels at 0.556. Beyond that noise score, the site features add **+0.072 [+0.03, +0.11]**. (Before the v4 calibration fix, about half the site signal was noise structure.)

The signal is carried by **straight WT windows, burial and non-local contacts**. In the v4 fit the coefficients were, on real movers vs noise: WT bend −0.47 vs −0.05, burial direction +0.21 vs +0.06, non-local contacts +0.12 vs +0.02. The information that local sequence lacks lives, modestly, in the tertiary environment of the site. That is what the first run's core-subset result (0.48 → 0.58, never tested) suggested.

---

## How far can accuracy go?

Two measurements bound what any model can do here:

* **Label reliability.** Splitting each protein's crystals into two halves and rebuilding the labels on each gives two independent measurements of the same mutations. They agree in direction every time, but whether a borderline effect crosses the threshold is noisy (κ = 0.37). Even an oracle that knew each effect as well as half the data would reach only **AUC ≈ 0.75**, and 85% of the labels rest on a single mutant crystal.
* **Learning curve.** With 25 / 50 / 75 / 100% of the training families, AUC is 0.60 / 0.61 / 0.62 / 0.62. Growing the set from 207 to 319 proteins added about 0.03; more families now add almost nothing.

Against that ceiling, the final model (site + substitution logistic regression, trained with label-confidence weights) reaches **0.64 [0.61, 0.67]** family-out and **0.71** on confident labels. The best ensemble reaches 0.65. That is about 55% of the achievable margin above chance. Tuning regularization, a continuous target, full-atom descriptors, boosting and ESM-2 each move AUC by at most ±0.02. The remaining headroom is in the labels (replicate crystals), not in the model or the amount of data.

---

## Using the model

```bash
python predict.py train                                   # once: fits results/mover_model.pkl
python predict.py score --pdb 2LZM --chain A --mut L99A,T26E,V149P,K16E
```

```
mutation    P(mover)  pctile  SS  WT bend  non-local  burial
T26E            0.42     94%  E     20.3          8      20
V149P           0.39     91%  H    113.5          4      21
L99A            0.35     81%  H    109.8          1      31
K16E            0.16     10%  E    101.9          5       4
```

`P(mover)` is the calibrated probability that the mutation moves the 5-residue window around it beyond crystal noise (training base rate 27%). `pctile` ranks it among the 1,250 training mutations. Use it to rank candidate mutations by risk of local backbone change, not as a yes/no call. Caveats: the shipped model is the v5 fit (bend label), and the pre-registered lockbox did not confirm that this kind of model generalizes to unseen proteins (*Final result*).

---

## Why this matters

This project is not an alternative to AlphaFold, and it was never intended to be.

Instead, it explores the negative space around modern structure prediction.

Successful protein models rely on long-range information because proteins themselves rely on long-range interactions.

By deliberately removing that information and measuring what remains, this project provides an empirical demonstration of why local sequence alone is insufficient.

The conclusion is simple:

> Mutation-induced backbone bending is real.
>
> Local sequence does not predict it: once the site is known, the substitution adds about 0.01 AUC, whether it is described physically or by an evolutionary language model.
>
> On the development data the site's structural context carries a modest signal, mostly through tertiary environment: AUC 0.66 (within-protein 0.65), against a label-reliability ceiling of 0.79. A pre-registered lockbox of 90 unseen proteins did not confirm it (AUC 0.53 [0.40, 0.66]); it was too small to detect an effect of this size, so the site signal remains an exploratory finding. At best it ranks mutations by risk; it is not a yes/no predictor.

That result may be less exciting than discovering a new predictor, but it is arguably more informative.

Knowing where the signal is not can be just as valuable as knowing where it is.

---

## Technical highlights

* 1,250 clean WT/mutant pairs from 319 proteins in 250 sequence families, mined systematically from RCSB + SIFTS
* Crystal-form matching, resolution cut, ligand-state matching on both sides (additives ignored), one row per mutation
* Calibrated noise floor: unbiased σ, exact median efficiency, B-factor-conditioned prior, Student-t mapping. WT-vs-WT false positives are 2.4% overall and ≤ 4.9% in every bin.
* WT-vs-WT null control for threshold calibration and noise-vs-movement decomposition
* Split-half test-retest of the labels, giving an explicit accuracy ceiling
* Full-atom site context, substitution × structure terms and ESM-2 language-model features, each tested with a paired ΔAUC
* Leave-site-out and leave-family-out validation (≥ 30% identity families)
* Residue-cluster bootstrap confidence intervals (mutations at the same site are resampled together)
* Paired ΔAUC tests for every model comparison
* All input sets pinned in `manifests/` for reproducibility
* Offline synthetic tests for every label-cleaning rule

The emphasis throughout was on falsification, uncertainty estimation, and honest interpretation rather than benchmark optimization.

---

## Reproducing

```bash
pip install -r requirements.txt
python bending_metric.py          # Gate 0 self-test (no network)
python feasibility_t4l.py         # Gate 1  (downloads T4L PDB entries)
python gate2_model_feasibility.py # Gate 2  (downloads the 30%-culled training set)
python mover_composition.py       # Gate 2b
python loop_gate.py               # Gate 3
python powered_loop_gate.py       # Gate 4  (downloads validation proteins)
python gate3_3d.py                # Gate 5
python summary_figure.py

# v2: clean labels + delta-targeted model (reuses the caches above)
python pairs.py                   # -> pairs_clean.csv + per-protein filter diagnostics
python delta_model.py             # -> results/delta_model.json, delta_model.png
python delta_model.py --z 3       # sensitivity: stricter mover threshold
python delta_model.py --prior pooled  # sensitivity: old single-sigma noise-floor prior

# v3: systematic miner (RCSB + SIFTS), ~6k entries, ~700 MB in mined_pdb/
python mine_pairs.py              # pinned in manifests/mined.json; pairs/delta_model pick it up

# v4: diagnostics, optional language-model features, model variants
python diagnostics.py             # null calibration by n_wt + split-half reliability -> results/diagnostics.json
python plm_features.py            # ESM-2 site log-probs -> results/esm_site_logp.csv (needs torch + fair-esm)
python model_variants.py          # learners, continuous target, learning curve -> results/model_variants.json
python predict.py train           # final model -> results/mover_model.pkl; then predict.py score ...

python -m pytest tests/           # offline synthetic tests (no network)
```

The scripts must run in this order: later gates reuse the PDB caches (`t4l_pdb/`, `cull_pdb/`, `val_pdb/`) that earlier gates download. On the first run, every RCSB search result is pinned to `manifests/*.json` (see `stats_utils.pinned_ids`). Commit those files so that later runs use the same entries, since live searches drift as the PDB grows. To refresh against today's PDB on purpose, delete a manifest.

> **Status of the numbers.** The Gate 1–5 figures quoted in `RESULTS.md` (and the first-run T4L numbers above) come from the original runs, which used a per-pair bootstrap. The scripts now use a residue-cluster bootstrap and a paired ΔAUC test, so CIs are expected to widen somewhat once the gates are re-run. Until then, treat the quoted CIs as optimistic. The headline numbers in this README come from the current (v4) code; the v2 and v3 blocks below are kept as history and used the older, liberal noise floor.

---

## v2: clean labels and a delta-targeted model

An audit of the Gate 1–5 pipeline found problems that bear on the headline conclusion:

* **Labels.** One T4L variant (L99A, around 60 ligand-soak crystals) made up 25% of the validation rows. Mutants were compared to WT crystals of any crystal form. Resolution and ligands were never checked. The mover threshold came from raw per-window MADs of as few as 3 crystals.
* **Engine.** Gate 2 regressed *absolute* bending (RMSE ~30°) and differenced two predictions to find changes with a median of ~2.7°. A tree ensemble returns exactly 0 unless a split touches the mutated position. The engine could not see the effect even if the information were there, so Gates 2–5 do not show that local sequence *lacks* the information.
* **Features.** In the existing T4L labels, site identity explains ~60% of |Δ| variance. A same-site diagnostic (not a valid predictor) reaches AUC 0.60, against 0.52 for the model. Where a mutation sits matters more than what it is, and the model had no site features.

`pairs.py` rebuilds the labels per SPEC §2.2–2.4. It applies a resolution cut (≤2.5 Å), compares each mutant only with WT crystals of the same crystal form, drops mutant crystals whose ligands near the window differ from the form's WT, and aggregates to one row per mutation. The noise floor is shrunk toward a prior that depends on window B-factor (flexible windows are noisier), and the mover threshold is a z-test on the median difference. `delta_model.py` then predicts movers directly from **site** features (SS, B-factor, burial, contacts, WT bend) and **substitution** features (Δvolume, Δhydrophobicity, charge, Pro/Gly, BLOSUM62, cavity × burial). It uses leave-site-out and leave-protein-out CV. The decisive test is the paired ΔAUC of *site+substitution* over *site*: does the substitution identity add anything once the site is known?

Both are covered by offline synthetic tests (`tests/`). Full numbers and robustness checks are in `RESULTS.md` (v2 section).

**v2 results on real data** (235 clean mutations, 75 movers, 121 sites, 4 proteins; SNase has only one WT crystal and drops out):

| model (logistic regression) | leave-site-out AUC | leave-protein-out AUC |
|---|---|---|
| site | 0.64 [0.55, 0.72] | 0.65 [0.56, 0.72] |
| substitution | 0.52 [0.45, 0.59] | 0.57 [0.50, 0.64] |
| site + substitution | 0.62 [0.53, 0.70] | 0.63 [0.54, 0.70] |

* **Substitution identity adds nothing once the site is known.** Paired ΔAUC (site+subst − site) = −0.02 [−0.06, +0.02]. Under all eight label definitions tried (z>2 or z>3; noise-floor prior pooled, B-factor, B-factor + WT bend, or none) it is ≤ +0.01 with CIs spanning 0. This is the local-sequence question asked directly, with an engine that can see a Δ, and the answer is still no.
* **Where the mutation sits carries a modest signal.** The main drivers are WT window bend (straighter windows move more), window B-factor and non-local contacts. A B-factor-conditioned noise floor removes the obvious confound (B correlates with WT noise, ρ = 0.51), and the signal stays. Adding WT bend to the prior as well weakens it to 0.60 [0.50, 0.69] at z>2 (0.68 at z>3). More proteins would be needed to separate flexibility-driven noise from flexibility-driven movement.
* Gradient boosting overfits at this sample size, so the tables report logistic regression. Details and all robustness checks are in `RESULTS.md`.

**v3: 191 proteins from a systematic RCSB + SIFTS miner** (`mine_pairs.py`; 892 clean mutations, 335 movers, 158 sequence families):

| model (logistic regression) | leave-family-out AUC |
|---|---|
| secondary structure only | 0.51 [0.47, 0.55] |
| site | 0.60 [0.56, 0.63] |
| site + substitution | 0.59 [0.56, 0.63] |

* **"What" (the substitution, i.e. local sequence) has no edge:** ΔAUC −0.002 [−0.03, +0.02].
* **"Where" has a small, real edge.** It beats SS by +0.085 [+0.05, +0.12]. A WT-vs-WT null control (held-out WT crystals scored as pseudo-mutants) shows that more than half of its margin above chance is predictability of WT noise (a noise-only score reaches 0.555). The movement-specific remainder is ΔAUC ≈ +0.05 [+0.01, +0.08], carried mainly by non-local contacts.
* The z threshold is close to calibrated: 6.9% pseudo-movers against 37.6% real movers.

**v4: audited labels** (1,045 clean mutations, 277 movers, 207 proteins, 166 families; full tables in `RESULTS.md`):

* **The v3 noise floor was liberal where few WT crystals exist.** Uncorrected MAD underestimates σ by a third at n = 3, so 17–20% of WT-vs-WT pseudo-mutants were "movers" in forms with 3–4 WT crystals. With an unbiased SD-based σ, the exact median efficiency and a Student-t mapping, the false-positive rate is 2.4% overall and ≤ 4.9% in every bin.
* Common crystallization additives no longer veto a crystal (+153 mutations, null unchanged). The ligand rule now also applies to the WT reference.
* With calibrated labels the site signal is movement, not noise structure: excess over a noise-only score +0.087 [+0.04, +0.13], up from +0.047.
* Substitution: +0.014 context-free, ≈ +0.01 net with substitution × structure terms, +0.003 with ESM-2.
* Ceiling: split-half oracle AUC ≈ 0.76. Best model 0.63 family-out.

**v5 (final): confidence-weighted training on an expanded set** (miner threshold ≥ 4 single mutants; 1,250 mutations, 332 movers, 319 proteins, 250 families):

| model (leave-family-out) | AUC, all labels | AUC, confident labels |
|---|---|---|
| secondary structure only | 0.54 [0.50, 0.57] | 0.55 |
| site | 0.63 [0.59, 0.66] | 0.68 |
| site + substitution (final model) | 0.64 [0.61, 0.67] | 0.71 |
| ensemble, all features | 0.65 [0.62, 0.68] | — |

* Substitution over site: +0.014 [−0.00, +0.03]. ESM-2: +0.002.
* Null 2.6%, split-half oracle 0.75, learning curve saturated (0.60 → 0.62).

---

## Final conclusion

The original hypothesis was that local amino-acid patterns drive mutation-induced backbone bending in a predictable way.

After multiple rounds of testing, the evidence does not support that hypothesis.

The signal exists.

The substitution does not predict it, and the site predicts it only weakly.

And that gap turns out to explain something important about protein structure itself.
