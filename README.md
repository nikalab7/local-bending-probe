# local-bending-probe

> ### ⚠️ Status: v1 below is frozen as a legacy / audit record
>
> An independent mathematical audit ([`AUDIT.md`](AUDIT.md)) found that the v1
> bending metric does not measure what the text below claims, and that the
> quantitative results are not currently supportable. Two defects in the metric
> are established by direct measurement ([`v2/`](v2/)):
>
> * **It is a secondary-structure classifier, not a curvature metric.** On a
>   perfectly *straight* axis it reads 110° for an α-helix, 132° for a π-helix,
>   67° for 3₁₀, 34° for PPII and 0° for a β-strand.
> * **It confuses twist with bend.** In an α-helix it responds to a change in
>   local twist at 1.026° per degree — about 3× its sensitivity to actual axis
>   bending — and its *sign* flips depending on where in the helical turn the
>   window sits, which scrambles any signed analysis (v1's Gate 2 reported
>   Spearman 0.077 and attributed it to the model).
>
> The following claims are **withdrawn pending v2** and should not be cited: the
> 29% above-noise-floor rate, AUC 0.52 / 0.59 / 0.70, the 136,961-window /
> ≤30%-identity leakage-controlled framing, and the conclusion that backbone
> bending is governed primarily by tertiary interactions.
>
> An earlier version of this banner additionally claimed v1's achievable AUC
> ceiling in helices was ~0.50, i.e. that no model could have beaten chance.
> **That claim was wrong** — it came from an averaging bug in the v2 gate, and is
> retracted in [`v2/CORRECTIONS.md`](v2/CORRECTIONS.md). Why v1 measured AUC 0.52
> is still open.
>
> Nothing has been deleted — the research trail is deliberately preserved, including
> the errors, mine as well as v1's. The v2 rebuild is on branch
> `v2-benchmark-rebuild`. The measurement gate **passed on ideal geometry**, which
> establishes the replacement metric is sound; whether the *label* is usable on real
> PDB data is **not yet established** and needs an empirical seed test.
>
> Everything below this banner is the original v1 write-up, unedited.

---

## How much information about protein backbone bending is contained in local sequence?

Modern protein-structure models rely heavily on long-range interactions and evolutionary information. This project asks a simpler question:

**If we deliberately remove all non-local information, how much can local amino-acid sequence alone tell us about mutation-induced backbone bending?**

To answer that, I built a lightweight, interpretable pipeline designed around a single idea:

> Hold local geometry constant, then test whether local sequence patterns can explain which mutations bend the backbone.

The original expectation was that specific sequence motifs would emerge as reliable local drivers of bending. Instead, the project arrived at the opposite conclusion.

Local sequence contains enough information to explain some aspects of absolute backbone geometry, but it contains remarkably little information about which mutations will change that geometry.

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

Across experimental structures, approximately 29% of mutations produced backbone changes larger than the measured structural noise floor.

Mutations do move protein backbones.

The project also confirmed a well-known structural principle:

> Backbone changes are more common in flexible regions.

Mutations were significantly enriched in loops compared with more rigid secondary structures.

These findings survived statistical testing and replication.

---

## What failed

The central hypothesis did not survive.

Models using only local sequence information performed only slightly above chance:

* Overall AUC ≈ 0.52
* Loop-focused replication AUC ≈ 0.59

More importantly, every confidence interval included chance performance.

The data therefore do not support the claim that local sequence can reliably predict mutation-induced backbone bending.

The result was consistent across multiple validation stages, datasets, and leakage-controlled evaluations.

---

## The most informative result

The strongest evidence came from introducing a small amount of non-local structural information.

When a simple description of the surrounding contact environment was added, performance improved in exactly the situations where protein physics predicts it should.

The effect was most visible in protein cores, where packing interactions dominate.

This suggests that the information missing from the local model is not hidden in more sophisticated sequence features.

It is largely absent from local sequence altogether.

Backbone bending appears to be governed primarily by tertiary interactions rather than local residue patterns.

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
> Structural context helps because structural context contains the information that local sequence lacks.

That result may be less exciting than discovering a new predictor, but it is arguably more informative.

Knowing where the signal is not can be just as valuable as knowing where it is.

---

## Technical highlights

* 136,961 training windows from 568 non-redundant protein chains
* Family-level holdout evaluation
* Sequence-identity culling
* Bootstrap confidence intervals
* Leakage-controlled validation
* Empirical noise-floor estimation
* Statistical enrichment analysis
* Explicit replication stages
* Structural-context ablation testing

The emphasis throughout was on falsification, uncertainty estimation, and honest interpretation rather than benchmark optimization.

---

## Final conclusion

The original hypothesis was that local amino-acid patterns drive mutation-induced backbone bending in a predictable way.

After multiple rounds of testing, the evidence does not support that hypothesis.

The signal exists.

The predictor does not.

And that gap turns out to explain something important about protein structure itself.
