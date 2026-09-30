# local-bending-probe

## How much information about protein backbone bending is contained in local sequence?

Modern protein-structure models rely heavily on long-range interactions and evolutionary information. This project asks a simpler question:

**If we deliberately remove all non-local information, how much can local amino-acid sequence alone tell us about mutation-induced backbone bending?**

To answer that, I built a lightweight, interpretable pipeline designed around a single idea:

> Hold local geometry constant, then test whether local sequence patterns can explain which mutations bend the backbone.

The original expectation was that specific sequence motifs would emerge as reliable local drivers of bending. Instead, the project arrived at the opposite conclusion.

On 892 clean WT/mutant pairs from 191 proteins, the identity of the substitution adds nothing to predicting which mutations move the backbone. *Where* the mutation sits carries a small signal, and its movement-specific part comes mainly from non-local contacts.

The failure of the local model became the result.

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

Every stage had a predefined failure condition.

The objective was not to maximize performance but to determine where the predictive information actually resides.

The first version of the pipeline (Gates 1–5, T4 lysozyme plus four validation proteins) is kept for the record. An audit showed its labels and its prediction engine were too weak to support a negative conclusion. The results below come from the rebuilt pipeline (v2/v3; details further down and in `RESULTS.md`).

---

## What worked

The phenomenon itself is real.

Across 191 proteins, 37.6% of clean single mutations (335/892) moved the backbone beyond |z| > 2 of the per-window noise floor. When held-out WT crystals are scored the same way as pseudo-mutants, only 6.9% cross that threshold. Mutations do move protein backbones.

The first T4 lysozyme run found 29% (72/248). That run also reported more movers in loops than in helices (48% vs 25%, one-sided Fisher p = 0.041). **That enrichment does not replicate on the clean pooled labels:** loops 37.3% (107/287), all other residues 37.7% (228/605), p = 0.97. A secondary-structure-only model is at chance (AUC 0.51 [0.47, 0.55]).

---

## What failed

The central hypothesis did not survive.

Substitution features (Δvolume, Δhydrophobicity, charge, Pro/Gly, BLOSUM62, cavity × burial) add nothing once the site is known:

* Paired ΔAUC, site + substitution − site: **−0.002 [−0.03, +0.02]** (leave-family-out)

This time the null is informative. The first engine predicted absolute bending (~30° error) and differenced two predictions, so it could not have seen a ~3° effect even if the information were there. The new model predicts movers directly and does find signal in site features. It finds none in the substitution.

The result holds across label definitions (z > 2 or z > 3; four noise-floor priors), leave-site-out and leave-family-out validation.

---

## The most informative result

Where the mutation sits does carry information:

* Site features: AUC **0.60 [0.56, 0.63]** with whole families held out
* Beyond secondary structure: **+0.085 [+0.05, +0.12]**

The null control shows that much of this is not movement. A site model trained only on WT-vs-WT noise reaches AUC 0.555 on the real labels, which is more than half of the margin above chance. Beyond that noise score, the site features add a movement-specific **+0.05 [+0.01, +0.08]**.

That remainder is carried mainly by **non-local contacts**. Their coefficient is +0.21 on real movers and −0.14 on noise. The information that local sequence lacks lives, weakly, in tertiary contacts, which is what the first run's core-subset result (0.48 → 0.58, never tested) suggested.

---

## Why this matters

This project is not an alternative to AlphaFold, and it was never intended to be.

Instead, it explores the negative space around modern structure prediction.

Successful protein models rely on long-range information because proteins themselves rely on long-range interactions.

By deliberately removing that information and measuring what remains, this project provides an empirical demonstration of why local sequence alone is insufficient.

The conclusion is simple:

> Mutation-induced backbone bending is real.
>
> Local sequence does not predict it: the substitution adds nothing once the site is known.
>
> The site's structural context carries a small real signal, mostly through tertiary contacts. At AUC ≈ 0.60 it is not a usable predictor.

That result may be less exciting than discovering a new predictor, but it is arguably more informative.

Knowing where the signal is not can be just as valuable as knowing where it is.

---

## Technical highlights

* 892 clean WT/mutant pairs from 191 proteins in 158 sequence families, mined systematically from RCSB + SIFTS (7,612 structures screened)
* Crystal-form matching, resolution cut, ligand-state matching, one row per mutation
* Heteroscedastic noise floor (σ prior conditioned on B-factor)
* WT-vs-WT null control for threshold calibration and noise-vs-movement decomposition
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

python -m pytest tests/           # offline synthetic tests (no network)
```

The scripts must run in this order: later gates reuse the PDB caches (`t4l_pdb/`, `cull_pdb/`, `val_pdb/`) that earlier gates download. On the first run, every RCSB search result is pinned to `manifests/*.json` (see `stats_utils.pinned_ids`). Commit those files so that later runs use the same entries, since live searches drift as the PDB grows. To refresh against today's PDB on purpose, delete a manifest.

> **Status of the numbers.** The Gate 1–5 figures quoted in `RESULTS.md` (and the first-run T4L numbers above) come from the original runs, which used a per-pair bootstrap. The scripts now use a residue-cluster bootstrap and a paired ΔAUC test, so CIs are expected to widen somewhat once the gates are re-run. Until then, treat the quoted CIs as optimistic. The v2/v3 numbers in this README come from the current code.

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

---

## Final conclusion

The original hypothesis was that local amino-acid patterns drive mutation-induced backbone bending in a predictable way.

After multiple rounds of testing, the evidence does not support that hypothesis.

The signal exists.

The substitution does not predict it, and the site predicts it only weakly.

And that gap turns out to explain something important about protein structure itself.
