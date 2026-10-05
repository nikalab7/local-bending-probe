# Final-evaluation protocol (frozen)

This file fixes the label, features, model, metrics and analysis plan
**before** the final evaluation. It was committed before `final_eval.py` was
run on either dev or lockbox data. Nothing below changes after that commit;
any later analysis is labelled exploratory. The decisions behind each choice
are in `WORKLOG.md` (decision log).

## 1. Data

| set | source | used for |
|---|---|---|
| **dev** | `pairs.load_proteins()`: T4L, SNase, barnase, human lysozyme, RNase A, and the mined proteins of `manifests/mined_ms4.json` (>= 4 single-substitution entities per UniProt accession) | every decision so far (label choice, thresholds, features, model, error analysis) |
| **lockbox** | `pairs.load_lockbox()`: `manifests/lockbox.json`, built by `python mine_pairs.py --lockbox` | the final evaluation only, once |

Lockbox construction (`mine_pairs.build_lockbox`):
- UniProt accessions with 1-3 single-substitution entities in the PDB (the dev miner required >= 4, so none of them was ever loaded), resolution <= 2.5 A, same pre-screen as dev (a crystal form with >= 3 WT crystals and a same-form single mutant);
- **excluded** if an RCSB sequence search (MMseqs2, >= 30% identity, e-value 1e-3) of its reference sequence hits any entity annotated with a dev accession, or any dev reference sequence; a failed search also excludes;
- lockbox families are >= 30%-identity clusters among lockbox proteins (prefixed `L`), built the same way as dev families.

Before this file was committed, only the lockbox manifest counts were seen (180 proteins, 178 families, 1334 entries, 256 single-substitution entities before label cleaning); no lockbox structure was parsed and no lockbox label was computed.

**Size rule.** The lockbox result is **confirmatory** only if the lockbox has >= 30 labelled mutations, >= 10 movers and >= 10 non-movers, and >= 10 families. Otherwise every result is reported as **exploratory**, and the dev numbers remain exploratory as well (dev data were used for all decisions).

## 2. Label (unchanged from WORKLOG v6)

- WT/mutant crystals of the same crystal form and ligand state; one row per mutation; 5-residue windows covering the mutated residue.
- Per window, four metrics: C-alpha bend, phi, psi, C-alpha torsion b; NCS copies with identical sequence averaged (`pairs.NCS_AVERAGE = True`).
- z per metric: (mutant - WT median) / sigma, sigma = SD/c4 of the WT crystals, B-factor-conditioned log-linear prior (K = 4), Student-t calibration.
- Combined label (`pairs.METRIC = "multi"`): z = the largest |z| of the four, rescaled so that the dev-null-calibrated threshold `MULTI_T = 2.61` maps to `Z_MOVER = 2`. **Mover = |z| > 2.**
- Reported with each run: WT-vs-WT null false-positive rate (pseudo-mutants through the same label code).

## 3. Model inputs (WT-only)

`FEATURES = SITE + SUBST` (`delta_model.py`), all computed from the WT structure and the substitution:
- SITE: SS one-hot (helix/strand/loop from HELIX/SHEET records), site and window B-factor z, C-alpha neighbours within 10 A, half-sphere exposure (up), non-local contacts, WT window bend, log distance to the chain terminus;
- SUBST: volume change (signed, absolute), hydrophobicity change, charge change, to/from Pro, to/from Gly, BLOSUM62, cavity (volume loss x burial), hydrophobicity change x burial.

**Forbidden as model inputs:** anything derived from the mutant structure or mutant experiment (resolution difference, temperature, altlocs, mutant contacts). Crystal contacts, temperature, resolution and altlocs are diagnostics only (section 6). `tests/test_pairs_model.py::test_features_are_wt_only` perturbs mutant coordinates, altlocs, resolution and temperature and asserts the features do not change.

Rejected for the final model (see decision log): full-atom context and substitution x context interactions, ESM-2 site log-probabilities and embeddings, elastic-network mechanics (dev dAUC -0.003 [-0.006, -0.000]), lattice contacts and WT altlocs as inputs, gradient boosting, ensembles.

On dev, mutant-side crystallographic covariates (resolution gap, temperature gap, mutant altlocs) would *raise* AUC by +0.018 [+0.006, +0.029]: part of the label is experimental artifact. They stay forbidden as inputs; section 6 measures how much the result depends on those cases.

## 4. Model

Logistic regression, `C = 0.3`, `class_weight="balanced"`, median imputation, standardization; training weights = label confidence, `clip(| |z| - 2 |, 0.25, 3)` (`delta_model.label_weights`). Baselines use the same model class and weights:
- **SS-only**: helix/strand/loop one-hot;
- **burial-only**: half-sphere exposure (up) + C-alpha neighbours within 10 A.

## 5. Splits

Families = >= 30%-identity sequence clusters (RCSB MMseqs2 search, linked only through exact reference-sequence hits). Dev: 10 folds of whole families (`delta_model.site_folds(family, 10, seed 0)`); `test_no_family_crosses_folds` asserts that no family crosses folds. Lockbox: the model is trained on all dev rows and scored on the lockbox; `final_eval.py` asserts that no family or protein is shared.

## 6. Metrics and analysis plan (`final_eval.py`)

All CIs are 95%, from 2000 bootstrap resamples **of families** (not mutations). Every "A beats B" claim uses a paired family-bootstrap test of dAUC (two-sided p).

**Headline:** lockbox ROC-AUC of the model with its family-bootstrap 95% CI, next to
- the label ceiling: split-half oracle AUC 0.787 (dev),
- the SS-only and burial-only baselines, with paired dAUC (model - baseline).
The model is said to work only if its lockbox AUC CI excludes 0.5 **and** the paired dAUC against both baselines has a CI excluding 0.

**Also reported (pre-specified):**
1. within-protein AUC (mover/non-mover pairs of the same protein only) + CI;
2. with and without T4L: dev — OOF AUC on non-T4L rows, and retrained without T4L; lockbox — model trained on dev with vs without T4L (paired);
3. positive control: same pipeline and splits predicting helix (from the non-SS SITE features + SUBST); expected AUC well above 0.5;
4. negative control: labels (and weights) shuffled within families — dev: retrained 20 times; lockbox: frozen scores, 1000 shuffles; expected AUC ~0.5;
5. power: Kish effective number of families, family-bootstrap SE of the AUC, design effect vs independent rows, minimum detectable AUC and dAUC vs SS-only (alpha 0.05 two-sided, power 0.8);
6. artifact diagnostics (never model inputs): AUC after excluding rows with a WT/mutant temperature mismatch > 50 K, mutant resolution worse by > 0.5 A, mutant altlocs in the window, or the mutated residue in a lattice contact (each, and "any" of these four); separately, any window residue in a lattice contact (broad: ~60% of dev rows); plus the mover rates of the suspect rows. The result "holds" if the AUC after excluding "any" stays within the headline CI and its own CI excludes 0.5;
7. lockbox WT-vs-WT null false-positive rate.

**Secondary:** confident labels (|z| > 3 vs < 1); helix / strand / loop.

**Small groups:** no result is reported for a group with < 30 rows or < 10 per class.

## 7. Order of operations

1. Commit this file, `final_eval.py`, `manifests/lockbox.json`.
2. `python final_eval.py dev` (dev numbers with the frozen pipeline).
3. `python final_eval.py lockbox` — once. The script refuses to run with uncommitted changes to any frozen file, or if `results/final_lockbox.json` already exists.
4. Report: README / RESULTS / WORKLOG updated with the numbers as they come out.
