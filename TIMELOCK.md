# Pre-registration: time-based lockbox

The first lockbox (`PROTOCOL.md`, `results/final_lockbox.json`) was inconclusive. With 111 mutations in 89 families, its minimum detectable AUC was 0.685. This file pre-registers a second, prospective test: the **same frozen pipeline**, applied once to PDB entries released after the cutoff. It changes nothing in the label, features, model or analysis plan of `PROTOCOL.md`. It adds only the data selection, a blinded power gate, and environment pinning.

**Pre-registration commit** = the commit that added this file:
`git log --diff-filter=A --format=%H -- TIMELOCK.md`. `timelock.py eval` refuses to run if any frozen file differs from that commit (see section 5).

## 1. Cutoff and data

- `CUTOFF = 2026-10-05`. Mutant crystals count only if their PDB initial release date is **after** the cutoff. WT reference crystals may be of any date.
- Selection reuses the dev miner unchanged: X-ray, <= 2.5 A, single-substitution entities with a UniProt reference, the same pre-screen (a crystal form with >= 3 WT crystals and >= 1 single mutant), and the same family definition (>= 30% identity, MMseqs2 via RCSB, exact-reference linking).

**Stratum A — new proteins (primary, confirmatory).** Accessions in neither `manifests/mined_ms4.json` (dev) nor `manifests/lockbox.json` (first lockbox). A protein is **excluded** if its reference sequence has a >= 30% sequence-search hit to any entity annotated with a dev or first-lockbox accession, or to any of their reference sequences; a failed search also excludes. Families are prefixed `T`. The question is the same as the first lockbox: does the model generalize to unseen protein families?

**Stratum B — new mutations of known proteins (secondary).** Dev and first-lockbox proteins (same reference construct), using only mutant crystals released after the cutoff. A row is dropped if its mutation (protein, position, WT and mutant residue) was already labelled in dev or the first lockbox; the list is stored in `results/frozen_dev_matrix.npz`. The model was trained on these proteins, so B tests generalization to **new mutations at known proteins**, not to new proteins. It can support only that weaker claim. Its key secondary metric is within-protein AUC.

## 2. Frozen model

Trained on the committed dev matrix `results/frozen_dev_matrix.npz` (1250 rows: z, family, protein, SITE + SUBST features; sha256 `c487c68af78732eafb0802df1dce67cf2b55ffa8566a1f9a6a401c2686421e5e`). It reproduces the first-lockbox score exactly (AUC 0.532099), so no dev rebuild is needed later. Label, features, model, baselines, controls and statistics are those of `PROTOCOL.md` and `final_eval.py`, unchanged.

## 3. Blinded power gate (lesson from the first lockbox)

The first lockbox shrank from 256 single-substitution entities to 111 labelled mutations after QC, and only then did its power turn out to be too low. Power is therefore computed **on the post-QC data, before unblinding**:

1. `python timelock.py build` (selection, download, sha256 of every file → `manifests/timelock_files.json`);
2. `python timelock.py gate --stratum S` builds the labels but reports **only counts**: rows, families, Kish n_eff, rows with >= 5 WT crystals, and the WT-vs-WT null false-positive rate. No mover label, rate or score is printed.
3. Predicted minimum detectable AUC = 0.5 + 2.80 × sqrt(DEFF × (n + 1) / (12 × q(1−q) × n²)), where q = 0.30 is the assumed mover rate and DEFF the design effect (A: 1.1, measured on the first lockbox; B: 0.9, measured on dev). With these values the formula gives SE 0.063 at n = 111; the observed SE was 0.066.
4. **Unblind only if the predicted MDA <= 0.62.** The gate is blinded, so it may be re-run as data accumulate (suggested yearly). It can never be re-run after the evaluation.

## 4. Evaluation (once per stratum)

`python timelock.py eval --stratum S`, run only after a passed gate. It calls `final_eval.run_lockbox` unchanged and writes `results/timelock_S.json`. It refuses to run a second time.

- **Success criterion (stratum A):** identical to `PROTOCOL.md`. The AUC 95% family-bootstrap CI excludes 0.5 **and** the paired dAUC against both the SS-only and the burial-only baselines has a CI excluding 0. Everything else in `PROTOCOL.md` section 6 is reported as well.
- **Null calibration (lesson from the first lockbox):** the dev-calibrated threshold did not transfer. Null FP was 2.6% on dev and 4.7% on the first lockbox, because label reliability depends on the number of WT crystals per form (median 7 vs 3). The threshold stays frozen. If the gate's null FP exceeds 4%, a pre-registered **secondary** analysis restricted to rows with >= 5 WT crystals is reported next to the primary one. The primary result is the all-rows analysis.
- **Stratum B:** the same output, read as "new mutations at known proteins". No claim about new proteins is made from B.

## 5. Environment pinning

- `requirements.lock`: exact versions of the pre-registration environment (Python 3.11.15, x86_64). `eval` refuses if numpy, scipy, scikit-learn or gemmi differ.
- `Dockerfile`: `python:3.11.15-slim` pinned by digest, with the lockfile installed. (The Docker daemon was not running in the pre-registration environment, so the image was not built there.)
- Frozen files, compared with the pre-registration commit by `eval`: `PROTOCOL.md`, `TIMELOCK.md`, `final_eval.py`, `pairs.py`, `delta_model.py`, `structure_features.py`, `plm_features.py`, `stats_utils.py`, `mine_pairs.py`, `timelock.py`, `requirements.lock`, both manifests, and `results/frozen_dev_matrix.npz`.
- Not pinnable: the RCSB search and data APIs (including the server-side MMseqs2 version) and PDB remediation of old entries. `build` records the sha256 of every downloaded file. Any API failure excludes the protein, as in the first lockbox.

## 6. Feasibility, stated plainly (`results/timelock_feasibility.json`)

From the release dates of all 9,662 dev and first-lockbox entries:
- **Rate of newly eligible families:** 18.5 per year (2016–2025, before QC).
- **Post-QC yield:** half of first-lockbox families produced labels, giving about 9 post-QC families per year at 1.25 rows per family.
- **Families needed** (alpha 0.05, power 0.8, q = 0.30, DEFF 1.1): about 107 for a predicted MDA of 0.66, and about 190 for 0.62.
- **Time needed:** about **12 years** for 0.66 and about **21 years** for 0.62.
- **300+ families within ~2 years is not realistic.** Two years would add about 18 post-QC families. Stratum A will very likely not pass its gate within this project's horizon; it is registered so that a future test is possible and cannot be tuned.
- **Stratum B accumulates faster** (125–230 new single-mutant entries per year in these proteins before QC, with an unknown share of new mutations). It may pass its gate within a few years, but it answers the weaker question.

A faster confirmatory route would need a new data source, not more time: for example, purpose-collected WT/mutant crystal series for many unrelated proteins.
