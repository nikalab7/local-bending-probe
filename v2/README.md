# v2 — Measurement Gate

**Status: gate PASSED. The label is usable. Proceed to benchmark v2.**

This directory is step 2 of the v2 rebuild: settle the *measurement* before any
model is trained again. No ML, no PDB downloads, no dataset — just geometry,
noise propagation, and an achievable-ceiling calculation.

Everything in `../` (v1) is frozen as the legacy/audit trail. Nothing was deleted.
`../AUDIT.md` records why v1 could not be patched incrementally.

```
python3 measurement_gate.py          # all 7 stages + verdict  (~10 min)
python3 measurement_gate.py 1 4      # selected stages
```

Files: `metrics.py` (candidate metrics + exact ground-truth generator),
`measurement_gate.py` (the gate).

---

## The question

> Can mutation-induced local geometry change be measured at all from this class of
> crystal data?

v1 answered "local sequence does not predict backbone bending." That conclusion is
not separable from "the measurement could not resolve the effect." The gate
separates them, in advance, with no model in the loop.

**Decision rule, fixed before running.** Let `tau*` be the true axis-bend SD a
metric needs for a *perfect* predictor to reach oracle AUC 0.75, at realistic
coordinate error. Reported as a mid-span Cα displacement (sagitta `L·θ/8`) so it
can be judged physically — point mutations move local backbone by ~0.1–0.5 Å.

* `tau*` displacement far above ~0.5 Å for every candidate → **measurement
  bottleneck**; the project's answer is about measurement, not sequence; do not
  resume the v1 ML gates.
* some candidate reaches 0.75–0.85 within ~0.5 Å → **label usable**; promote it and
  build benchmark v2.

## Ground truth

`helix_on_arc(n, r, rise, turn, total_bend_deg)` winds a regular helix (or strand,
or any twist) around an axis that is an exact circular arc turning by a known
angle. That is the calibration target v1 never had: v1's self-test used only
windows whose two halves were exactly collinear (`AUDIT.md` F12), so it could not
have distinguished any of these candidates.

---

## What the gate found

### 1. v1's metric is blind to axis bend in helices — its ceiling was 0.50

The decisive number. Response to *true* axis bend, and the resulting ceiling:

| metric | conformation | slope | noise @0.2 Å | oracle AUC, any τ ≤ 20° |
|---|---|---|---|---|
| `v1_pca5` | **α-helix** | **0.013** | 4.46° | **0.496 – 0.498** |
| `v1_pca5` | 3₁₀-helix | 0.049 | 7.11° | 0.493 – 0.510 |
| `v1_pca5` | PPII | 0.110 | 4.45° | 0.503 – 0.651 |
| `v1_pca5` | β-strand | 0.500 | 3.76° | 0.523 – 0.913 |

v1's metric responds to axis curvature **only in β-strands**. In helical and
turn-like geometry its slope is 5–40× smaller, and the oracle ceiling is ~0.50
*regardless of how large the true effect is* — `tau*` for α-helix is 280° (3.66 Å),
i.e. unreachable.

T4 lysozyme is helix-rich and v1's "core" subset was helix+sheet. So for most of
the windows v1 scored, **no model could have beaten chance.** The reported
AUC 0.477–0.524 is the ceiling, not a model failure. This is the mechanism behind
`AUDIT.md` F8, now localized to a specific geometric cause.

### 2. A straight axis must read zero — only the bisector family does

Worst |metric| on a perfectly **straight** axis (target 0):

| candidate | α-helix | β-strand | π-helix | 3₁₀ | PPII | verdict |
|---|---|---|---|---|---|---|
| `v1_pca5` | 110.4° | 0.0° | 132.4° | 66.8° | 33.7° | FAIL |
| `smooth_chord` w=4 | 11.1° | 0.0° | 7.1° | 11.6° | 5.1° | FAIL |
| `smooth_pca` w=4 | 11.9° | 0.0° | 7.5° | 10.6° | 4.7° | FAIL |
| `circlefit` w=4 | 172.9° | 167.7° | 168.8° | 141.0° | 79.5° | FAIL |
| `bisector` | **0.000** | abstains | **0.000** | **0.000** | **0.000** | PASS |
| **`hybrid`** | **0.000** | **0.000** | **0.000** | **0.000** | **0.000** | **PASS** |

v1's metric is not merely "helix vs strand" — it reads a *different* large offset
for every regular conformation. It is a secondary-structure classifier with a
curvature-shaped name.

`smooth_pca ≈ smooth_chord` throughout, which settles a question the audit left
open: **smoothing was the missing ingredient, the PCA step was incidental.**

### 3. Why smoothing alone cannot fix it, and what does

Sliding-window averaging shrinks the coil by the Dirichlet factor
`|sin(wθ/2) / (w·sin(θ/2))|`, which depends on the twist θ. No single width
collapses every conformation (residual axis wobble, Å):

| w | α-helix | β-strand | π-helix | 3₁₀ | PPII |
|---|---|---|---|---|---|
| 3 | 0.500 | 0.317 | 0.968 | **0.000** | **0.000** |
| 4 | 0.257 | **0.000** | 0.100 | 0.475 | 0.320 |
| 7 | 0.074 | 0.136 | 0.450 | 0.271 | 0.183 |

The bisector construction sidesteps it. For a regular helix the Cα bisector
`b_i = û(Ca[i-1]-Ca[i]) + û(Ca[i+1]-Ca[i])` points radially inward, so two
consecutive bisectors are radial vectors at different heights and
`cross(b_i, b_{i+1})` lies **along the axis for any twist** — exact, not
approximate. Its one failure is honest and reported: at twist ≈180° (β-strand)
bisectors go antiparallel and it **abstains** — precisely where w=4 smoothing is
itself exact.

`HybridAxisBend` dispatches on that degeneracy condition (not on an external SS
annotation, so it is not circular) and divides by a per-branch calibration, so both
sides report the same quantity: **true full-span axis bend in degrees**.

* bisector calibration is essentially universal over helical twists —
  0.6248 / 0.6256 / 0.6307 / 0.6310 for α / π / 3₁₀ / PPII, all **R² = 1.0000**
* hybrid slope vs ground truth: 0.995 (α), 1.004 (3₁₀), 1.005 (PPII), 1.250 (β)
* continuous twist sweep 60°–185°: **no coverage gap**, straight axis reads 0.000
  throughout (worst 1.8° in the handover zone at 175°)

### 4. The real decision axis is differential bias, not noise

A constant offset cancels in a WT→mutant difference, which is the only way v2 uses
the metric. What survives is the offset's *sensitivity to conformation*: if a
mutation nudges the local twist, a twist-coupled offset leaks straight into the
delta and is indistinguishable from real bending.

d(metric)/d(twist) on a straight axis — degrees of **spurious** bend per degree of
twist change:

| candidate | α-helix | β-strand | 3₁₀ | PPII |
|---|---|---|---|---|
| `v1_pca5` | **1.026** | 0.000 | **1.428** | 0.800 |
| `smooth_chord` w=4 | 0.656 | 0.000 | 0.116 | 0.050 |
| `circlefit` w=4 | 0.880 | 0.000 | 6.213 | 7.513 |
| `bisector` / `hybrid` | **0.000** | **0.000** | **0.000** | **0.000** |

Put the two sensitivities side by side and v1's metric stops being ambiguous:

| conformation | d/d(axis bend) | d/d(twist) | **twist : bend** |
|---|---|---|---|
| α-helix | 0.013 | 1.026 | **79 ×** |
| 3₁₀-helix | 0.049 | 1.428 | **29 ×** |
| PPII | 0.110 | 0.800 | 7 × |
| β-strand | 0.500 | 0.000 | 0 × |

**In an α-helix, v1's metric responds ~79× more strongly to a change in local
twist than to actual axis bending.** So whatever it measured in helical windows was
essentially twist change wearing the name "bending" — and re-twisting a helix is
exactly what changing helix propensity by mutation does. This is a confound on
v1's "29% of mutations move backbones" entirely separate from the noise and
threshold problems in `AUDIT.md` F2/F4, and it points the same way.

This is why the noisier metric wins. `smooth_chord` is ~2.5× quieter than the
hybrid (1.2–2.7° vs 3.7–7.8° at σ=0.2 Å), but it carries a 0.656 twist coupling in
α-helices. **Noise averages down with more crystals; differential bias does not.**
So the hybrid is the right trade for a difference measurement.

### 5. The threshold must be solved for, not assumed

v1 fixed the mover cut at `2·σ̂` and never checked its null rate (`AUDIT.md` F2).
The correct multiplier depends on the aggregation and is not 2:

| n_WT | n_mut | SD(Δ) | k for 5% FPR | FPR if k=2 |
|---|---|---|---|---|
| 1 | 1 | 1.41 | 2.77 | 15.7% |
| 5 | 2 | 0.95 | 1.86 | 3.4% |
| 20 | 5 | 0.63 | 1.23 | 0.4% |

(units of the single-measurement σ). v2 computes `k` per variant pair.

### 6. Verdict and its sensitivity

At σ_xyz = 0.20 Å, 5 WT / 2 mutant crystals, 5% null FPR — `tau*` for oracle
AUC 0.75, as Cα displacement:

| candidate | β-strand | α-helix | 3₁₀ | PPII |
|---|---|---|---|---|
| **`hybrid`** | **0.15 Å** | **0.16 Å** | 0.19 Å | 0.29 Å |
| `bisector` | abstains | 0.18 Å | 0.19 Å | 0.29 Å |
| `smooth_chord` | 0.18 Å | 0.40 Å | 0.42 Å | 0.24 Å |
| `v1_pca5` | 0.18 Å | **3.66 Å** | **2.14 Å** | 0.93 Å |

All within the plausible 0.1–0.5 Å band for the hybrid → **LABEL USABLE.**

(These are Monte-Carlo estimates and move by ±0.01–0.02 Å run to run — quoted to
two decimals for comparability, not because the third would be meaningful. The
verdict has a wide margin against the 0.5 Å bound, so the wobble does not touch it.
The one number where it would matter is `v1_pca5` on α-helix, and that one is off
by an order of magnitude, not a rounding.)

And it does not depend on generous aggregation. `tau*` displacement for the
hybrid across the sensitivity grid:

| σ_xyz | 1WT/1mut | 3WT/1mut | 5WT/2mut | 20WT/5mut |
|---|---|---|---|---|
| 0.10 Å | 0.14 Å | 0.12 Å | 0.09 Å | 0.06 Å |
| 0.20 Å | 0.26 Å | 0.22 Å | 0.17 Å | 0.12 Å |
| 0.30 Å | 0.37 Å | 0.30 Å | 0.23 Å | 0.15 Å |

(α-helix; β-strand is similar except 0.58 Å at 0.30 Å / 1WT-1mut.)

Even **one crystal per side** clears the bar at σ ≤ 0.2 Å. That matters, because
redundancy is scarce — from v1's own CSV, median 2 crystals per residue, and only
17% of residues have ≥5:

| crystals per residue | ≥1 | ≥2 | ≥3 | ≥5 | ≥10 |
|---|---|---|---|---|---|
| residues (of 69) | 100% | 62% | 38% | 17% | 3% |

So crystals-per-variant is a resolution *bonus*, not a gate. Worth noting v1 spent
its redundancy the wrong way round — inflating n to 248 rows from 69 residues
(`AUDIT.md` F3) instead of averaging it down into the label.

---

## What this does and does not establish

**Does.** The measurement is not the bottleneck — *provided* the metric is an axis
metric. v1's negative result is explained by its metric being blind to axis bend in
helical geometry (slope 0.013, ceiling 0.50) and by converting twist change into
fake bending at ~1:1. Neither is a fact about protein sequence.

**Does not.** Every number here is exact or perturbed *ideal* geometry. That is the
right scope for a gate — it bounds what the measurement can do independent of any
dataset — but the true effect size τ is empirical. The gate says which τ matters
and how precisely it must be known; it does not say what τ is. Closing that is step
3 of benchmark v2.

Known residuals, not swept under the rug:
* hybrid underestimates bend by ~18% in π-helix (single calibration constant
  across helical twists); tighten with a twist-dependent constant if π-helices
  matter
* the β branch is accurate to ~30° of bend and degrades above it (reads 81.7° for a
  true 60°); fine for the 1–10° regime being chased, not for large bends
* 9-residue span costs coverage vs v1's 5, and is less "local" — both need
  reporting as benchmark coverage
* 1.8° straight-axis residual in the twist handover zone near 175°

## Next

1. Promote `HybridAxisBend` to metric v2 — already calibrated to report true
   full-span axis bend in degrees, with a `branch()` QC field for the output schema.
2. Build benchmark v2: exact RCSB entity/chain + SIFTS mapping → canonical
   sequence → exact variant → median aggregation over redundant structures →
   matched WT/mutant pairs → independent QC → family-level split by *joint*
   sequence clustering (`AUDIT.md` F5 — v1's holdout was PDB-ID only).
   **One biological variant = one observation.** Bootstrap at variant/protein level.
3. Estimate τ empirically and re-read the stage-6 table. That closes the loop.
4. Only then baselines, each read against the ceiling:
   `SS-only → local sequence → direct-Δ → +3D context`.
   If SS-only already explains most of it, the "local model skill" was
   structure-class prediction. **Only a direct-Δ model that fails while the
   ceiling is high is a real negative ML result.**
