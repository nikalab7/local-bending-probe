# Work log: what was done and what we found

Short summary. Full tables: `RESULTS.md` (v2–v5 and the final pre-registered evaluation), `README.md`, `PROTOCOL.md`.

## The question
Does a single amino-acid substitution bend the protein backbone at that site, and can we predict which mutations will?
- **"What":** which amino acid was replaced by which (the local-sequence change).
- **"Where":** the structural environment of the site (burial, contacts, shape).

A mutation counts as a **mover** when its change clearly exceeds crystal-to-crystal noise. The noise is measured by comparing several WT crystals of the same protein with each other.

## Stages
| stage | what was done | best AUC |
|---|---|---|
| original (main) | 248 windows, essentially one protein (T4L) | 0.52 (CI included chance) |
| v2 | clean labels (same crystal form, ligand state, resolution); site features | 0.64 (4 proteins, wide CI) |
| v3 | systematic miner (RCSB + SIFTS), 191 proteins; WT-vs-WT null control | 0.60 (half the signal turned out to be noise structure) |
| v4 | **fixed the noise floor**: small-n MAD bias, SD estimator, Student-t; additives; WT ligand rule; split-half ceiling; ESM-2 | 0.61 (ensemble 0.63) |
| v5 | confidence-weighted training; 319 proteins; `predict.py` | 0.639; 0.714 on confident labels |
| v6 | NCS-copy averaging; label-metric comparison; **combined label** | 0.662; 0.74 on confident labels |
| final (pre-registered) | `PROTOCOL.md` frozen (commit `fe9c3c7`); fresh lockbox of 90 unseen proteins scored once | **lockbox 0.532 [0.403, 0.661]** (not confirmed, underpowered); dev 0.662 [0.633, 0.698] |

## Main results
**Status: all results are exploratory.** The pre-registered lockbox did not confirm the dev signal (AUC 0.532 [0.403, 0.661], n = 111, 89 families), and it was too small to detect it (minimum detectable AUC 0.685). See *Final result* below.

1. **Bending is real.** 28–36% of clean single mutations move the backbone beyond noise, against 2.6% of WT-vs-WT pseudo-mutants.
2. **"What", bounded (post-hoc, dev):** for the C-alpha bend label, adding substitution features to site features gives +0.004 [−0.015, +0.025]: no detectable gain, upper bound +0.025 AUC. Label noise (kappa 0.37) attenuates effects, so this does not show absence. For the frozen combined label, which includes phi/psi, substitution features add +0.059 [+0.033, +0.089], partly outside Gly/Pro. ESM-2 and substitution x context add nothing beyond them. (Earlier versions said "+0.014, rejected", which was the v5 bend-label number.)
3. **On dev, the information is in "where":** a straight WT window, burial and non-local contacts. That beats secondary structure by +0.09 and burial alone by +0.11, and beats a noise-only score by +0.07. About 0.03 of the pooled dev AUC is between-protein (negative control 0.53); within-protein AUC is 0.65. **Not confirmed on the lockbox.**
4. **The ceiling is set by label noise.** Split-half reliability puts the achievable AUC at about 0.72–0.79. The learning curve is saturated, and every model variant (boosting, ensemble, ESM) lands within ±0.02.

## v6: better labels (null false-positive rate ~2.6–3% for all)
| label | model AUC | confident labels | κ | ceiling (oracle) |
|---|---|---|---|---|
| bend (old) | 0.639 | 0.714 | 0.31 | 0.724 |
| bend + NCS | 0.646 | 0.713 | 0.37 | 0.771 |
| φ + NCS | 0.648 | 0.745 | 0.33 | 0.734 |
| ψ + NCS | 0.611 | 0.635 | 0.35 | 0.770 |
| **bend + φ + ψ + CA torsion (NCS, combined)** | **0.662** | **0.736** | **0.51** | **0.787** |

- **NCS copies** (several copies of the protein in one crystal) clean the label without any new data.
- **The combined label** finds more real movers at the same false-positive rate (36% vs 28%). It is also more reliable and more predictable.

## Error analysis (clustering)
- "Missed" helical movers are not a blind spot. Helices simply move less often (23% vs 33% in strands), and within helices the model ranks well (AUC 0.65).
- Separate per-SS models and SS × feature interactions lower AUC (overfitting).
- One small but real blind spot remains: Gly at positive φ (9 cases, 56% missed).

## Comparison with other work
| | task | result |
|---|---|---|
| Schaefer & Rost 2012 | fragments from **different** proteins, sequence-based | AUC 0.80, but an easier task with no noise control |
| AlphaFold3 (Liu, Calabrese, O'Hern 2026, arXiv 2609.24842) | WT/mutant, full structure prediction | ρ≈0.55; ≈0.2 on large changes; worse outside its training set |
| AlphaFold2 (McBride et al., PRL 2023) | effective strain | correlates on average |
| **this project** | WT/mutant, WT structure only, strict family hold-out | AUC 0.66; 0.74–0.76 on confident labels; ρ≈0.23; ≈0.19 on large changes |

What is different here:
- calibrated labels (crystal form, ligand state, NCS, null control, split-half ceiling);
- prediction from the WT structure alone;
- the "what vs where" decomposition.

An exact comparison needs AlphaFold run on this dataset, which requires a GPU.

## Files
- **labels:** `pairs.py`, `mine_pairs.py`
- **features:** `delta_model.py`, `structure_features.py`, `plm_features.py`
- **diagnostics:** `diagnostics.py`, `model_variants.py`
- **ready-to-use tool:** `predict.py`
  ```bash
  python predict.py score --pdb 2LZM --chain A --mut L99A,T26E
  ```

## Decision log

Every decision that shaped the final pipeline: what was decided, why, on which data, and what was rejected. "Dev" = the data of `pairs.load_proteins()` (T4L, the four validation proteins and the mined proteins); every decision below was made on dev data, so dev numbers are exploratory. The lockbox (PROTOCOL.md) was not used for any of them.

### Labels
| # | decision | why | data | rejected alternatives |
|---|---|---|---|---|
| L1 | Compare a mutant only with WT crystals of the **same crystal form** (space group + cell within tolerance) | different lattices bend windows differently; mixing forms inflated "movers" | T4L + validation proteins (v2) | all WT crystals pooled |
| L2 | **Ligand state** must match near the window (8 A); WT filter "prefer" (applied when >= 3 clean WT remain); common additives ignored | ligands next to the window move it; additives (SO4, Cl, GOL, EDO, BME...) did not add null false positives (3.1% vs 3.9%) | dev, v2-v4 ligand audit | mutant-only rule (asymmetric); "strict" WT filter (loses too many forms); counting additives as ligands (dropped ~190 crystals for nothing) |
| L3 | **One row per mutation** (mutant crystals of a mutation aggregated) | T4L L99A soaks alone gave 61/248 rows in the original set | T4L (v2) | one row per mutant crystal |
| L4 | Resolution <= 2.5 A on both sides | low-resolution backbones add noise | v2 | 3.0 A cut |
| L5 | sigma = **SD / c4** | 1.4826 x MAD is biased low at n = 3-4 (17-20% null FP in those bins); MAD/Qn ignore minority WT states (5-7% null FP vs 3-4% for SD) | dev WT-vs-WT null (v4) | MAD, corrected MAD, Qn |
| L6 | **B-factor-conditioned prior** (log-linear in window B, K = 4 shrinkage) | few-crystal sigma is unstable; flexible windows are noisier | dev (v2-v4) | pooled prior (made site-only hgb score below chance and produced spurious contrasts); no shrinkage; WT bend as a second prior covariate (confounds the strongest feature) |
| L7 | Exact small-sample **median efficiency** for the mutant side, **Student-t** calibration (nu = n - 1 + K) | asymptotic 1.2533 is wrong at m = 1 (85% of rows); sigma uncertainty gives t tails (6-8% null FP) | dev null (v4) | normal calibration |
| L8 | **NCS averaging** of identical-sequence copies (>= 80% coverage) | free replicates: split-half kappa 0.31 -> 0.37, oracle AUC 0.724 -> 0.771 | dev (v6) | first chain only |
| L9 | **Combined label**: max \|z\| over bend, phi, psi, CA torsion b; threshold 2.61 so that the WT-vs-WT null FP equals the single-metric rate (~2.6%) | most reliable (kappa 0.51, oracle 0.787) and most predictable (AUC 0.662) of the five labels compared | dev (v6) | bend only (0.639 / 0.724), bend + NCS (0.646 / 0.771), phi + NCS (0.648 / 0.734), psi + NCS (0.611 / 0.770), CA torsion a; sum or mean of z over metrics; thresholds 2 and 3 on the raw max |
| L10 | Mover = \|z\| > 2 (after rescaling) | conventional; null-calibrated | dev | 3 (too few movers), continuous \|z\| as primary target (kept only as a ridge variant) |

### Data and splits
| # | decision | why | data | rejected alternatives |
|---|---|---|---|---|
| D1 | Systematic miner: UniProt accessions with **>= 4 single-substitution entities**, SIFTS mapping, pre-screen for a usable form | the gate proteins were too few (4 proteins, wide CIs) | RCSB (v3, v5) | hand-picked proteins only; >= 5 entities (v3, 191 proteins; lowered to 4 in v5 for 319 proteins) |
| D2 | **Families** = >= 30%-identity clusters (RCSB MMseqs2 search), linked only through exact reference-sequence hits | linking through any hit merged unrelated proteins via chimeric constructs | dev (v3 fix) | PDB-ID or accession as the hold-out unit; CATH/Pfam (not tried: sequence clustering covers every protein, including constructs without a domain assignment) |
| D3 | **Leave-family-out**, 10 folds of whole families, as the primary CV | site-out leaks family-level information | dev | leave-site-out (reported as secondary until v5), leave-one-protein-out |
| D4 | **Lockbox** = accessions with 1-3 single-substitution entities (never loaded by the dev miner) and no >= 30% hit to any dev protein; one run | every dev family was used in some decision (label comparison, feature choice, error analysis), so none is clean | RCSB (this protocol) | re-using held-out dev folds (all were looked at); a time split (not tried) |

### Features
| # | decision | why | data | rejected alternatives |
|---|---|---|---|---|
| F1 | **SITE** (C-alpha site context) | beats SS by +0.07 to +0.09 (paired); signal survives the noise decomposition (+0.07 beyond a noise score) | dev (v3-v5) | SS only (at chance) |
| F2 | **SUBST** (context-free substitution descriptors) kept | small but consistent gain (+0.014 [-0.00, +0.03]) and needed to state the "what vs where" result. *Post-hoc, frozen combined label: +0.059 [+0.033, +0.089]; bend label +0.004 [−0.015, +0.025]* (`results/posthoc_dev.json`) | dev (v4, v5; post-hoc on v6) | site only |
| F3 | Full-atom CONTEXT and substitution x context INTERACT **not** in the final model | where - site = -0.023 [-0.04, -0.00] (logreg); the best variant (hgb/ensemble on where + what, 0.649) was chosen after the fact among ~12 variants | dev (v4, v5) | where, where + what (logreg and hgb), ensembles |
| F4 | ESM-2 site log-probabilities and LLR **rejected** | +0.002 [-0.01, +0.01] | dev (v4, v5) | site + ESM terms |
| F5 | ESM-2 per-residue **embeddings rejected** | PCA 8/16/32/64 components over site + subst: -0.009 to +0.001, all CIs span 0; embeddings alone 0.591 (-0.057) | dev, bend + NCS label (this round) | any PC count |
| F6 | ENM mechanics (ANM site MSF, window bend response) **rejected** | over site + subst: -0.003 [-0.006, -0.000]; single-feature AUCs 0.45-0.51 | dev, combined label (this round) | ENM features alone or with lattice/altloc (-0.008 [-0.015, -0.002]) |
| F6b | WT altloc fraction and lattice contacts **not inputs** | -0.001 [-0.003, +0.001] and -0.004 [-0.010, +0.001]; lattice is also an experiment property | dev, combined label | as inputs |
| F6c | Mutant-side covariates (resolution gap, temperature gap, mutant altlocs) **forbidden** although they help | +0.018 [+0.006, +0.029]: mutants collected at room temperature against cryo WT move in 57% vs 34% of cases, mutants >0.5 A worse in 55% vs 34% -- part of the label is experiment, not protein; using them would predict the experiment | dev, combined label | as inputs (user rule 3); kept as exclusion diagnostics |
| F7 | Lattice contacts, temperature, resolution gap, altlocs are **diagnostics only**, never inputs | user rule 3: mutant-side information is forbidden; lattice/temperature/altlocs describe the experiment, not the protein | — | using them as covariates (would leak experiment-specific information) |

### Model and training
| # | decision | why | data | rejected alternatives |
|---|---|---|---|---|
| M1 | **Logistic regression, C = 0.3**, balanced classes, median imputation, standardization | nested CV picks the same C (0.000 difference); boosting is not better at this n | dev (v4, model_variants) | nested-C logreg, ridge on \|z\|, HGB, ensembles, per-SS models and SS x feature interactions (both lowered AUC) |
| M2 | **Confidence weights** clip(\| \|z\| - 2 \|, 0.25, 3) | split-half: borderline labels (\|z\| near 2) flip between independent crystal sets, \|z\| < 1 or > 3 rarely do | dev (v4 split-half; adopted in v5) | unweighted; dropping borderline rows |
| M3 | Baselines **SS-only** and **burial-only** (hse_up + n_ca10), same model class | user rule 5 | — | — |

### Evaluation
| # | decision | why | data | rejected alternatives |
|---|---|---|---|---|
| E1 | Family-bootstrap 95% CIs; paired family-bootstrap dAUC for every comparison | rows within a family are correlated; earlier runs used residue-cluster 90% CIs (v4 showed they agree with family resampling) | — | row bootstrap; residue-cluster bootstrap |
| E2 | Headline = lockbox AUC next to the 0.787 ceiling and both baselines; within-protein AUC; with/without T4L; positive (helix) and negative (shuffled within family) controls; power analysis; artifact exclusions as diagnostics | user rules 5-7 | — | — |
| E3 | No result from groups with < 30 rows or < 10 per class | user rule 6 | — | — |
| E4 | Artifact "any" exclusion = temperature gap > 50 K, mutant resolution > 0.5 A worse, mutant altlocs, mutated residue in a lattice contact; window-in-contact reported separately | window contacts cover ~60% of dev rows and are unrelated to mover rate (36% vs 36%), so they would remove most data for nothing | dev (counts and mover rates only, no model scores) | window contact inside "any" |

## Final result (pre-registered)
Rules followed (user's protocol message): freeze before evaluation, a clean lockbox, WT-only inputs, family splits, controls and baselines, family bootstrap with paired tests and power, fixed reporting, and this decision log.

| | lockbox (once) | dev (exploratory) |
|---|---|---|
| model AUC | **0.532 [0.403, 0.661]** | 0.662 [0.633, 0.698] |
| model − SS-only | +0.002 [−0.099, +0.117] | +0.087 [+0.052, +0.120] |
| model − burial-only | +0.083 [−0.038, +0.204] | +0.112 [+0.076, +0.160] |
| within-protein | 0.375 [0.00, 0.83] | 0.653 [0.635, 0.711] |
| positive / negative control | 0.951 / 0.538 | 0.915 / 0.529 |
| minimum detectable AUC | 0.685 | 0.545 |
| null false positives | 4.7% | 2.6% |
| ceiling | 0.787 (dev) | 0.787 |

- **Not confirmed:** the lockbox CI includes 0.5 and the model does not beat the baselines there.
- **Underpowered, not refuted:** the lockbox could detect only AUC >= 0.685; its CI also contains 0.66. Its labels are noisier (median 3 WT crystals per form vs 7).
- **Without T4L:** dev 0.666; lockbox model trained without T4L 0.538. T4L does not drive the result.
- **Artifacts:** excluding temperature/resolution/altloc/lattice-suspect rows raises dev AUC to 0.70. These cases add noise; they are not the signal.

Decisions recorded after the lockbox run (not changes to the frozen analysis):
| # | what | why | data | rejected alternatives |
|---|---|---|---|---|
| P1 | All results labelled exploratory | the protocol's size rule (>= 30 rows, >= 10 per class, >= 10 families) was met, but the power analysis shows the lockbox cannot detect the dev effect, so user rule 2 ("too small -> exploratory") applies | lockbox power analysis | — |
| P2 | The lockbox negative control is reported as degenerate | 73 of 89 lockbox families have one row, so shuffling within families leaves most labels unchanged | lockbox family sizes | — |
| P3 | No re-analysis of the lockbox (no new features, thresholds or subsets) | one run only; anything further would be exploratory and would contaminate the next lockbox | — | re-scoring with other features or thresholds |

| P4 | "What" claim bounded and split by label | post-hoc paired tests: bend label +0.004 [−0.015, +0.025]; combined label +0.059 [+0.033, +0.089]; phi +0.041, psi +0.043, CA torsion +0.027; without Gly/Pro +0.031 [+0.010, +0.054]. The earlier README statement (+0.01) used the v5 bend-label value and was wrong for the frozen label | dev, frozen label (`posthoc_dev.py`) | the unbounded wording "local sequence doesn't determine bending" |
| P5 | Artifact result framed as "quantified in this dataset", not new | cryo-cooling and resolution effects on apparent structural differences are documented (Fraser 2011, Halle 2004, Keedy 2015/2018, Juers & Matthews 2001, Cruickshank 1999); n and family-bootstrap CIs reported (room-temperature mutant vs cryo WT: 0.565 [0.477, 0.684] vs 0.337, n = 115, 20 families; resolution > 0.5 A worse: 0.545 [0.456, 0.646] vs 0.341, n = 101, 33 families) | dev | — |
| P6 | Methods lessons recorded | (a) null FP 2.6% -> 4.7% on the lockbox: reliability depends on WT crystal count (median 7 vs 3); (b) 256 entities -> 111 mutations after QC: power must be computed post-QC, before unblinding; (c) within-family shuffles are degenerate with singleton families; (d) the label definition decides the "what" answer | dev + lockbox counts | — |
| P7 | Time-based lockbox pre-registered (`TIMELOCK.md`): cutoff 2026-10-05; stratum A = new proteins (confirmatory), B = new mutations of known proteins (secondary, weaker question); blinded post-QC gate, unblind only if predicted MDA <= 0.62 | lessons (a)/(b); keeps the frozen pipeline untouched | release dates of 9,662 dev/lockbox entries | re-analysing the opened lockbox; merging lockbox into dev; a lower gate (would repeat the underpowered test) |
| P8 | Environment pinned: `requirements.lock`, digest-pinned `Dockerfile`, pre-registration commit `662994a` (the commit that added TIMELOCK.md), frozen dev matrix `results/frozen_dev_matrix.npz` (reproduces lockbox AUC 0.532099 exactly) | the run in years must use the same pipeline; rebuilding dev labels later could drift with PDB remediation | — | trusting "same code" without version and file checks |
| P9 | Feasibility stated plainly | ~18.5 newly eligible families/yr (2016–2025), 50% post-QC yield -> ~9/yr; ~190 needed for AUC 0.62 (~21 yr), ~107 for 0.66 (~12 yr); 300+ within ~2 years is not realistic | `results/timelock_feasibility.json` | — |

| P10 | Lab / refinement test (`LABTEST.md`, pre-registered `94e9bb7`) | cross-lab WT–WT null 0.068 [0.052, 0.094], not ~20%; within-lab, matched-condition mutational excess +0.259 [+0.207, +0.317]; PDB-REDO paired change +0.034 [−0.052, +0.118] (n = 59); pre-registered decisions: INCONCLUSIVE (primary), INCONCLUSIVE (PDB-REDO) | dev + first lockbox, 8,450 null pseudo-mutants | — |
| P11 | Null bias recorded | holding out all WT crystals gives 5.3% WT–WT FP vs 2.6% with the 10 best-resolution ones; the frozen threshold is somewhat permissive | dev + first lockbox | — |

## Next
- Write-up: `SUMMARY.md` (methods + bounded negative result; exploratory except the lockbox).
- Prospective test: `TIMELOCK.md`. Re-run `python timelock.py build` and `gate` yearly. Run `eval` once per stratum, only after a passed gate. Stratum A will very likely not pass within ~2 years.
- A faster confirmatory route needs a new data source, e.g. purpose-collected WT/mutant crystal series for many unrelated proteins.
- `predict.py` still ships the v5 (bend-label) fit; refit with the frozen label if the tool is used.
