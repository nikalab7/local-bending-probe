# Work log: what was done and what we found

Short summary. Full tables: `RESULTS.md` (v2–v5) and `README.md`.

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
| v6 (current) | NCS-copy averaging; label-metric comparison; **combined label** | **0.662**; 0.74 on confident labels |

## Main results
1. **Bending is real.** 28–36% of clean single mutations move the backbone beyond noise, against 2.6% of WT-vs-WT pseudo-mutants.
2. **"What" adds almost nothing:** +0.014 AUC, and +0.002 even with the ESM-2 language model. The original hypothesis, that local sequence determines bending, is rejected.
3. **The information is in "where":** a straight WT window, burial and non-local contacts. That beats secondary structure by +0.09 and beats a noise-only score by +0.07.
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
| F2 | **SUBST** (context-free substitution descriptors) kept | small but consistent gain (+0.014 [-0.00, +0.03]) and needed to state the "what vs where" result | dev (v4, v5) | site only |
| F3 | Full-atom CONTEXT and substitution x context INTERACT **not** in the final model | where - site = -0.023 [-0.04, -0.00] (logreg); the best variant (hgb/ensemble on where + what, 0.649) was chosen after the fact among ~12 variants | dev (v4, v5) | where, where + what (logreg and hgb), ensembles |
| F4 | ESM-2 site log-probabilities and LLR **rejected** | +0.002 [-0.01, +0.01] | dev (v4, v5) | site + ESM terms |
| F5 | ESM-2 per-residue **embeddings rejected** | PCA 8/16/32/64 components over site + subst: -0.009 to +0.001, all CIs span 0; embeddings alone 0.591 (-0.057) | dev, bend + NCS labels (this round) | any PC count |
| F6 | ENM mechanics (ANM site MSF, window bend response) — ENM_DECISION | ENM_WHY | dev, bend + NCS labels (this round) | ENM_REJECTED |
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

## Running / next
- ESM-2 embedding test (running)
- Artifacts: crystal contacts, data-collection temperature, resolution gap, altlocs (running)
- Mechanical features (elastic network: how much the window bends when the site is pushed) (running)
- Build the combined label into `pairs.py`, re-run the full pipeline, retrain the final model
