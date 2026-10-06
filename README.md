# local-bending-probe

## Can local sequence tell which mutations change local backbone geometry?

> **Status.** This is the first version of the project (T4 lysozyme plus four validation
> proteins). It is kept as a record. An earlier version of this README concluded that
> information about mutation-induced bending is absent from local sequence. **The data
> here do not support that conclusion**, and it is withdrawn. The section
> [What can be claimed](#what-can-be-claimed) states what does hold. A rebuilt pipeline
> is in progress on branch `v2-benchmark-rebuild`. It uses per-mutation labels, a
> calibrated noise floor and a replacement metric.

---

## The question

A protein backbone can change shape for many reasons: local residue preferences,
secondary-structure tendencies, packing, long-range contacts and solvent effects. This
project tried to isolate the first one:

**If all non-local information is removed, can local amino-acid sequence predict which
single mutations change the backbone around them?**

The design was a chain of falsification tests ("gates"), each with a failure condition
fixed before it ran. `RESULTS.md` is the gate-by-gate log with the original numbers.

---

## What was done

1. **Metric** (`bending_metric.py`). The angle between principal axes of the two halves of
   a 5-residue Cα window.
2. **Signal** (`feasibility_t4l.py`). WT/mutant crystal pairs of T4 lysozyme. A pair is
   called a *mover* when |Δangle| exceeds 2σ of a per-window floor estimated from WT
   crystals.
3. **Model** (`gate2_model_feasibility.py`). Gradient boosting (HistGradientBoosting),
   trained to predict the **absolute** angle from local sequence on 136,961 windows from
   568 chains culled at ≤ 30% identity. Mutation effects were scored as
   prediction(mutant) − prediction(WT).
4. **Loops and replication** (`loop_gate.py`, `powered_loop_gate.py`). Loop-only retrieval on
   T4L, then on barnase, human lysozyme and RNase A.
5. **3D context** (`gate3_3d.py`). Contact-residue composition added to the features.

---

## What the data show

**The model over-responds.** On the 248 T4L pairs, predicting Δ = 0 for every pair gives
an RMSE of 3.02°. The model's predicted Δ has an RMSE of 7.58° against the observed Δ,
2.5× worse than predicting zero. Its predictions do not track the observed changes
(Spearman 0.08) and are roughly twice their typical size. Mover retrieval is
AUC 0.52. The likely reason: the model learned average residue-to-local-geometry
propensities across proteins, and inside a fixed folded context most of those effects
do not happen.

> `RESULTS.md` (Gate 2) describes this the other way round: a "5.7× error cancellation"
> and a model that "barely responds". Both readings are wrong. The √2×RMSE bound assumes
> independent errors, which two near-identical inputs do not have. Against the
> predict-zero baseline the model is too sensitive, not insensitive.

**The model is weak on absolute geometry, too.** Held-out RMSE is 30.55° against a standard
deviation of 37.97°, so R² ≈ 0.35. (`RESULTS.md` calls the 20% RMSE reduction a
"variance reduction".)

**The metric measures local Cα conformation, not bending.** An ideal, perfectly straight
α-helix reads about 110° and a straight β-strand about 0°. The angle is mostly a
secondary-structure readout, and a "mover" is a pair whose local Cα geometry changed,
not one whose backbone axis bent. The median |Δ| of movers (2.66°) corresponds to Cα
shifts of roughly 0.1–0.2 Å, close to the coordinate error of 1.5–2 Å structures.

**The 248 pairs are not 248 independent observations.** They cover 69 residues. Residue 99
alone contributes 61 pairs (the L99A cavity series with different ligands), and those
pairs scatter like noise. Bootstrap intervals that resample pairs are therefore too
narrow, and the effective number of independent movers is closer to the 38 distinct
residues than to 72.

**Label reliability is unknown.** The rule's false-positive rate on WT-vs-WT pairs was not
measured. Replicates of the same variant often disagree on mover status. The WT floor
(0.75–0.98°) comes from isomorphous WT crystals and probably underestimates noise for
mutants crystallized under other conditions. Without a test-retest ceiling, AUC 0.52
cannot be read as "no signal".

**There is no positive control.** T4L mutagenesis targets stability and cavities, mostly
core hydrophobics. Pro/Gly are involved in only 3% of movers, so the cases where local
sequence is known to matter (X→Pro in a helix, Gly in a turn) are nearly absent.
Nothing in this set could have shown a local effect if one existed.

**The loop enrichment does not survive.** On all pairs, loops have more movers than helices
(Fisher p = 0.041), but the test ignores repeated residues. Gate 4 replicated loop
**AUC**, not the enrichment.

**The 3D lift is untested.** Core-subset AUC moved from 0.48 to 0.58, but no paired test on
the difference was run, and refit noise (±0.02–0.03) is of the same order. It is
suggestive at most.

---

## What can be claimed

> On T4 lysozyme, a gradient-boosting model trained on absolute local geometry across
> proteins does not identify which single mutants change local Cα geometry. Its
> predicted changes do not correlate with the observed ones and are about twice as large.

That is narrow, but it holds. The following do **not** follow from this data:

- that local sequence contains no information about mutation-induced backbone change.
  This is a low-power null from one weak model, with unknown label reliability and no
  positive control;
- that backbone bending is governed primarily by tertiary interactions;
- that ~29% of mutations move the backbone. That figure is 72/248 crystal pairs in one
  protein, with an uncalibrated threshold.

---

## Known gaps in this version

- Analysis unit is the crystal pair, not the variant. `SPEC_bending_and_pairs.md` §2.4
  asks for per-variant aggregation; no script implements it.
- No label-reliability ceiling and no null calibration of the mover rule.
- No positive control.
- Validation proteins are held out by PDB ID, not by sequence identity.
- The confound controls in the SPEC (ligand state, crystal contacts, resolution) are not
  implemented in the T4L pipeline.
- The Gate 2 training set is the first 600 hits of a live RCSB query and changes as the
  PDB grows. IDs are not pinned.

---

## Running it

Python 3 with `numpy`, `scipy`, `scikit-learn` and `matplotlib`. Every gate except
Gate 0 downloads structures from RCSB.

```bash
pip install numpy scipy scikit-learn matplotlib
python bending_metric.py           # Gate 0 self-test (offline)
python feasibility_t4l.py          # Gate 1
python gate2_model_feasibility.py  # Gate 2
python mover_composition.py        # Gate 2b
python loop_gate.py                # Gate 3
python powered_loop_gate.py        # Gate 4
python gate3_3d.py                 # Gate 5
python summary_figure.py
```
