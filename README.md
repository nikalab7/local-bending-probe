# local-bending-probe

## How much information about protein backbone bending is contained in local sequence?

Modern protein-structure models rely heavily on long-range interactions and evolutionary information. This project asks a simpler question:

**If we deliberately remove all non-local information, how much can local amino-acid sequence alone tell us about mutation-induced backbone bending?**

To answer that, I built a lightweight, interpretable pipeline designed around a single idea:

> Hold local geometry constant, then test whether local sequence patterns can explain which mutations bend the backbone.

The original expectation was that specific sequence motifs would emerge as reliable local drivers of bending. Instead, the project arrived at the opposite conclusion.

Local sequence contains enough information to explain some aspects of absolute backbone geometry, but at the sample sizes available here it showed no detectable information about which mutations will change that geometry.

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

Given two protein segments with similar local geometry, can local sequence alone predict which one bends more?

If the answer were yes, it would suggest that interpretable local rules explain a meaningful fraction of backbone deformation.

If the answer were no, it would imply that the information lives elsewhere.

---

## Approach

The project was designed as a sequence of falsification tests rather than a search for positive results.

The workflow was:

1. Build a robust bending metric.
2. Verify that mutation-induced bending exists in real structures.
3. Train a local-sequence model.
4. Test whether prediction survives strict validation.
5. Add structural context and measure what changes.

Every stage had a predefined failure condition.

The objective was not to maximize performance but to determine where the predictive information actually resides.

---

## What worked

The phenomenon itself is real.

In T4 lysozyme, 29% of single mutations (72/248) produced backbone changes larger than twice the measured per-window noise floor (a threshold that pure noise would cross roughly 5% of the time).

Mutations do move protein backbones. The 29% is best read as a ceiling rather than a typical rate: T4 lysozyme is unusually mutation-tolerant, its mutagenesis is core-biased, and noise floors estimated from sparsely sampled windows inflate the above-floor fraction.

The project also confirmed a well-known structural principle:

> Backbone changes are more common in flexible regions.

In T4 lysozyme, mutations in loops were more often movers than mutations in helices (48% of loop mutations were movers vs 25% in helices; OR 2.50 for loop vs all non-loop, one-sided Fisher p = 0.041).

The noise floor was confirmed by two independent estimates (0.98° and 0.75°). The loop enrichment is a single-protein result at p = 0.041 and has not been replicated on other proteins.

---

## What failed

The central hypothesis did not survive.

Models using only local sequence information performed only slightly above chance:

* Overall AUC ≈ 0.52
* Loop-focused replication AUC ≈ 0.59

More importantly, every bootstrapped confidence interval included chance performance. (The overall 0.52 was not bootstrapped in the original run; the scripts now compute a CI for it too.)

The data therefore do not support the claim that local sequence can reliably predict mutation-induced backbone bending.

The result was consistent across multiple validation stages, datasets, and leakage-controlled evaluations.

---

## The most informative result

The strongest evidence came from introducing a small amount of non-local structural information.

When a simple description of the surrounding contact environment was added, the AUC moved in the direction, and in the place, that protein physics predicts.

The shift was largest in protein cores, where packing interactions dominate (0.48 → 0.58).

This is suggestive, not established. The improvement itself (the paired difference between the two models) was not tested in the original run, and the core subset has no confidence interval. The scripts now compute both.

If it holds up, it suggests the missing information is not hidden in more sophisticated local sequence features but lives in tertiary contacts, meaning backbone bending would be governed mainly by tertiary interactions rather than by local residue patterns.

---

## Why this matters

This project is not an alternative to AlphaFold, and it was never intended to be.

Instead, it explores the negative space around modern structure prediction.

Successful protein models rely on long-range information because proteins themselves rely on long-range interactions.

By deliberately removing that information and measuring what remains, this project provides an empirical demonstration of why local sequence alone is insufficient.

The conclusion is simple:

> Mutation-induced backbone bending is real.
>
> Local sequence does not reliably predict it.
>
> Structural context appears to help, consistent with it carrying the information that local sequence lacks, but that lift is not yet statistically established.

That result may be less exciting than discovering a new predictor, but it is arguably more informative.

Knowing where the signal is not can be just as valuable as knowing where it is.

---

## Technical highlights

* 136,961 training windows from 568 non-redundant protein chains
* Family-level holdout evaluation
* Sequence-identity culling
* Residue-cluster bootstrap confidence intervals (mutations at the same site are resampled together)
* Paired ΔAUC test for model comparisons
* Leakage-controlled validation
* Empirical noise-floor estimation
* Statistical enrichment analysis
* Explicit replication stages
* Structural-context ablation testing

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
```

The scripts must run in this order: later gates reuse the PDB caches (`t4l_pdb/`, `cull_pdb/`, `val_pdb/`) that earlier gates download. On the first run, every RCSB search result is pinned to `manifests/*.json` (see `stats_utils.pinned_ids`). Commit those files so that later runs use the same entries, since live searches drift as the PDB grows. To refresh against today's PDB on purpose, delete a manifest.

> **Status of the numbers.** The figures in this README and in `RESULTS.md` come from the original runs, which used a per-pair bootstrap. The scripts now use a residue-cluster bootstrap and a paired ΔAUC test, so CIs are expected to widen somewhat once the gates are re-run. Until then, treat the quoted CIs as optimistic.

---

## Final conclusion

The original hypothesis was that local amino-acid patterns drive mutation-induced backbone bending in a predictable way.

After multiple rounds of testing, the evidence does not support that hypothesis.

The signal exists.

The predictor does not.

And that gap turns out to explain something important about protein structure itself.
