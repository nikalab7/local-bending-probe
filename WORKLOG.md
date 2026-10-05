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

## Running / next
- ESM-2 embedding test (running)
- Artifacts: crystal contacts, data-collection temperature, resolution gap, altlocs (running)
- Mechanical features (elastic network: how much the window bends when the site is pushed) (running)
- Build the combined label into `pairs.py`, re-run the full pipeline, retrain the final model
