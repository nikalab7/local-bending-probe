# v2 — Measurement Gate

**Status: gate passed on the tested ideal conformations outside a newly identified
twist handover band, where the metric now abstains. Whether the LABEL and its
coverage are usable on real PDB data is NOT yet established — that needs the
empirical seed test.** See [`CORRECTIONS.md`](CORRECTIONS.md) for withdrawn claims.

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
`measurement_gate.py` (the gate), `perturbation.py` (controlled deformation and
local Jacobian for **real** coordinates — the tool the seed test needs),
`CORRECTIONS.md` (withdrawn claims).

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

### 1. v1's metric gives the same bend opposite signs depending on helical phase

> **Retraction.** The first version of this section claimed v1's oracle ceiling in
> α-helices was ~0.50, so "no model could have beaten chance". **That was wrong**
> — it came from an averaging bug in my own gate (phase-averaging a signed
> response cancelled it). Full retraction in [`CORRECTIONS.md`](CORRECTIONS.md).
> Why v1 measured AUC 0.52 remains open.

What the corrected measurement shows. Per-phase signed slope vs true axis bend on
an ideal α-helix:

| metric | signed slope range | sign flips? | mean \|slope\| |
|---|---|---|---|
| `v1_pca5` | −0.476 … +0.443 | **yes** | 0.316 |
| `smooth_chord` w=4 | −0.220 … +0.302 | **yes** | 0.177 |
| **`hybrid`** | **+0.862 … +1.120** | **no** | **0.997** |

v1's metric assigns the *same physical bend* a positive or a negative value
depending on where in the helical turn the window sits. That has two distinct
consequences, and conflating them is what produced the retracted claim:

* a **magnitude** score (`|Δ|`, which is what mover AUC uses) is unaffected in
  principle and gets the mean|per-phase| gain — corrected ceiling for v1 on
  α-helix is 0.600 at τ=5° and 0.719 at τ=10°, not 0.50;
* any **signed** analysis is scrambled by a nuisance parameter. v1's Gate 2
  reported `Spearman(pred, obs) = 0.077` and read it as model insensitivity
  ("AF2 disease") — a phase-dependent metric sign produces that number on its
  own, with no model defect required.

The hybrid is phase-consistent, so it is still the better metric — 0.695 vs 0.600
at τ=5° — but that is a quantitative improvement, not a qualitative rescue.

The mean|per-phase| value above is a descriptive summary. Stage 6 keeps the
full phase-slope distribution when it simulates noisy labels; replacing it with
one mean slope would misprice the thresholded AUC.

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

**Handover correction:** the original direct switch was not smooth. A fine twist
sweep found a 3.45° jump on a perfectly straight axis at about 171.4° twist.
The metric now abstains between the admissible bisector and near-180° smoothed
domains; the coverage cost on real structures has not been measured. See
[`CORRECTIONS.md`](CORRECTIONS.md) C4.

* bisector calibration is essentially universal over helical twists —
  0.6248 / 0.6256 / 0.6307 / 0.6310 for α / π / 3₁₀ / PPII, all **R² = 1.0000**
* hybrid slope vs ground truth: 0.995 (α), 1.004 (3₁₀), 1.005 (PPII), 1.250 (β)
* the earlier 10° twist sweep missed a narrow handover failure; the corrected
  metric reports `NaN` in that region instead of a misleading bend

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

Put the two sensitivities side by side (using the **corrected** bend slopes —
mean|per-phase|, per §1 and `CORRECTIONS.md`):

The table uses mean|per-phase| as a descriptive sensitivity. Stage 6 uses the
individual phase slopes instead, because thresholding a noisy delta is nonlinear.

| conformation | d/d(axis bend) | d/d(twist) | **twist : bend** |
|---|---|---|---|
| α-helix | 0.316 | 1.026 | **3.2 ×** |
| 3₁₀-helix | 0.287 | 1.428 | **5.0 ×** |
| PPII | 0.294 | 0.800 | 2.7 × |
| β-strand | 0.500 | 0.000 | 0 × |

**In an α-helix, v1's metric responds ~3× more strongly to a change in local twist
than to actual axis bending** — so twist change is its dominant systematic, and
re-twisting a helix is exactly what changing helix propensity by mutation does.
This is a confound on v1's "29% of mutations move backbones" separate from the
noise and threshold problems in `AUDIT.md` F2/F4.

> An earlier version of this table quoted **79×** and **29×**, using the
> phase-cancelled bend slopes (0.013, 0.049). Those ratios are withdrawn; the
> `d/d(twist)` column is unaffected, because on a straight axis the metric has no
> phase dependence to cancel.

This is why the noisier metric wins. `smooth_chord` is ~2.5× quieter than the
hybrid (1.2–2.7° vs 3.7–7.8° at σ=0.2 Å), but it carries a 0.656 twist coupling in
α-helices, and a twist-coupled systematic is not something more data fixes.

> **Correction.** This section previously asserted flatly that "noise averages down
> with more crystals; differential bias does not". That holds only for **iid**
> noise. Real crystal-to-crystal variation also has components correlated across
> crystals of one variant — packing and space group, refinement protocol and
> deposition era, construct/background, cryo vs RT, ligand state — so the honest
> model is `σ²(n) = σ²_iid/n + σ²_systematic`, which tends to `σ_systematic`, not 0.
> The size of that floor is unknown and only estimable from real redundant
> crystals. See `CORRECTIONS.md` C2.

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

At σ_xyz = 0.20 Å, 5 WT / 2 mutant crystals, 5% null FPR — the coverage-aware
stage-6 simulation gives `tau*` for oracle AUC 0.75, as Cα displacement:

| candidate | β-strand | α-helix | 3₁₀ | PPII |
|---|---|---|---|---|
| **`hybrid`** | **abstains** | **0.17 Å** | 0.20 Å | 0.30 Å |
| `bisector` | abstains | 0.19 Å | 0.20 Å | 0.30 Å |
| `smooth_chord` | 0.18 Å | 0.40 Å | 0.42 Å | 0.24 Å |
| `v1_pca5` | 0.18 Å | **3.66 Å** | **2.14 Å** | 0.93 Å |

The admissible hybrid branches are within the plausible 0.1–0.5 Å band. The
β-strand row is an abstention, so coverage is part of the result rather than a
hidden finite-only filter.

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

On these **synthetic, iid-noise-only** numbers even one crystal per side clears
the bar at σ ≤ 0.2 Å. Do not read that as benchmark readiness: the 1/√n scaling
assumes iid noise, and the systematic floor (`CORRECTIONS.md` C2) is not in this
table because it cannot be estimated without real redundant crystals. Treat the
row as a sensitivity result, not a green light.

Redundancy is scarce anyway — from v1's own CSV, median 2 crystals per residue,
and only 17% of residues have ≥5:

| crystals per residue | ≥1 | ≥2 | ≥3 | ≥5 | ≥10 |
|---|---|---|---|---|---|
| residues (of 69) | 100% | 62% | 38% | 17% | 3% |

Under iid noise alone, crystals-per-variant would be a resolution bonus rather
than a gate. With a systematic floor it may well be a gate — that is one of the
things the seed test has to settle, and it is why the seed test should stratify
redundant crystals by crystal form rather than pooling them. Worth noting v1 spent
its redundancy the wrong way round either way: inflating n to 248 rows from 69
residues (`AUDIT.md` F3) instead of averaging it down into the label.

---

## What this does and does not establish

**Does.** On the ideal conformations tested outside the handover region, the
replacement metric reads zero on a straight axis, is phase-consistent and tracks
controlled bend. And v1's metric has two demonstrable defects — large
SS-dependent offsets, and a phase-dependent sign that scrambles signed analyses
(which is a candidate explanation for its Spearman 0.077). Neither is a fact about
protein sequence.

**Does not.** Three things, and they are load-bearing:

1. **It does not establish that the label is usable on real data.** Every number
   is exact or iid-perturbed *ideal* geometry. Real helices are pre-bent, frayed,
   irregular in rise and twist, with correlated non-Gaussian coordinate error.
   Whether these slopes and noise figures transfer is unmeasured.
2. **It does not pin τ.** The gate says which τ matters and how precisely it must
   be known; it does not say what τ is. That is empirical.
3. **It does not explain v1's AUC 0.52.** The ceiling explanation is retracted
   (§1). That question is reopened, not answered.

It also probed only one deformation mode (a distributed arc). `perturbation.py`
adds a midpoint kink that preserves every Cα-Cα distance exactly; the two modes
turn out to give similar sensitivity for v1 (ratio 0.7–1.2× across conformations),
so this particular worry did not materialise — but which mode real point mutations
produce is open, and it changes what the right metric is.

Known residuals, not swept under the rug:
* hybrid underestimates bend by ~18% in π-helix (single calibration constant
  across helical twists); tighten with a twist-dependent constant if π-helices
  matter
* the β branch is accurate to ~30° of bend and degrades above it (reads 81.7° for a
  true 60°); fine for the 1–10° regime being chased, not for large bends
* 9-residue span costs coverage vs v1's 5, and is less "local" — both need
  reporting as benchmark coverage
* the corrected hybrid abstains near the twist handover; real-data coverage is
  unmeasured

## Seed-test implementation status

`seed_test.py` keeps the estimated crystal-form systematic offset in the null
threshold and simulated AUC; it does not average that offset away as the number
of crystals grows. Thresholds use the actual WT/mutant form labels and counts,
so a shared form is correlated across the pair instead of being treated as two
independent offsets. It checks the mutation identity and the full 9-residue
local sequence against the WT reference before scoring a structure. It accepts
fixed-column `.pdb` and `.ent` files; mmCIF is not parsed.

The estimator measures a real-coordinate Jacobian for every usable WT site and
uses that slope distribution in tau deconvolution and the oracle ceiling. It
reports metric coverage, branch counts, and cross-branch pair rejections, and
uses a site-cluster bootstrap for the tau confidence interval. The harness is
explicitly single-protein; manifests containing multiple `protein_id` values
are rejected until a protein/site hierarchy is implemented.

The automatic verdict is **inconclusive** when fewer than 20 variants are
scored, crystal counts vary by pair, a positive cluster-bootstrap tau lower
bound is absent, or the iid/systematic split cannot be estimated. No real seed
set has been run here.

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
