# Corrections to the v2 measurement gate

Three claims made in the first version of `v2/README.md` and in commit `85844c2`
were overstated or wrong. Two were my own analysis errors; one was an unstated
assumption. They are recorded here rather than edited out of history, because the
whole point of this project's structure is that the trail stays visible.

Prompted by review feedback that the gate was drawing a strong conclusion from
ideal geometry too quickly. That review was right, and chasing it down found a
harder error underneath.

---

## C1 — WITHDRAWN: "v1's α-helix ceiling was 0.50, so AUC 0.52 was the ceiling"

**This was wrong, and the cause was a bug in my own gate.**

The gate computed each metric's response to axis bend by first averaging the
metric over eight helical phases, then fitting a slope against true bend:

```python
vals = [mean over phases of metric(window(bend=b))  for b in BENDS]
slope = lstsq(bends, vals)          # <-- cancellation happens here
```

For v1's metric on an ideal α-helix the per-phase signed slopes are:

| phase | 0° | 45° | 90° | 135° | 180° | 225° | 270° | 315° |
|---|---|---|---|---|---|---|---|---|
| slope | +0.443 | +0.232 | −0.128 | −0.454 | −0.476 | −0.170 | +0.197 | +0.430 |

- mean of signed slopes: **+0.009** ← what the gate reported as "0.013"
- mean of \|signed slopes\|: **+0.316** ← the actual per-window sensitivity

The response **changes sign with helical phase**. Averaging the signed value
before fitting cancelled it to near zero, understating v1's sensitivity by ~34×.
That near-zero slope was then fed into the oracle-ceiling calculation, which
duly returned ~0.50 — a spurious result produced entirely by the averaging order.

Corrected ceilings for α-helix (5 WT / 2 mutant crystals, 5% calibrated FPR):

| | slope | noise | τ=2° | τ=5° | τ=10° | τ=20° |
|---|---|---|---|---|---|---|
| v1_pca5 — as reported (wrong) | 0.013 | 4.46° | 0.496 | 0.506 | 0.495 | 0.505 |
| **v1_pca5 — corrected** | **0.316** | 4.46° | **0.515** | **0.600** | **0.719** | **0.839** |
| hybrid (unaffected) | 0.995 | 7.82° | 0.549 | 0.695 | 0.821 | 0.910 |

So v1's metric was **not** at a 0.50 ceiling in helices. The hybrid is still
better — 0.695 vs 0.600 at τ=5° — but that is a quantitative improvement, not the
qualitative "no model could have beaten chance" I claimed. **Any explanation of
v1's AUC 0.52 as a ceiling effect is withdrawn.**

The hybrid's numbers are unaffected, because the hybrid's response *is* phase
consistent (+0.86 to +1.12 across all phases), so nothing cancelled.

Fixed in code: `phase_slopes()` now returns the per-phase slopes, both means, and
a `sign_consistent` flag; stage 2 prints all of them; stage 6 uses mean|per-phase|.

### What replaces it — a real finding, narrower

The phase dependence is itself worth having, and it survives:

| metric | α-helix signed slope range | sign flips? | mean\|slope\| |
|---|---|---|---|
| `v1_pca5` | −0.476 … +0.443 | **yes** | 0.316 |
| `smooth_chord` w=4 | −0.220 … +0.302 | **yes** | 0.177 |
| **`hybrid`** | **+0.862 … +1.120** | **no** | **0.997** |

v1's metric assigns the *same physical bend* a positive or negative sign
depending on where in the helical turn the window sits. Consequences:

* a **magnitude** score (`|Δ|`, as used for mover AUC) still works, at the
  mean|per-phase| gain — hence C1 above;
* any **signed** analysis is scrambled by a nuisance parameter. v1's Gate 2
  reported `Spearman(pred, obs) = 0.077` and read it as model insensitivity
  ("AF2 disease"). A phase-dependent metric sign would produce that number on
  its own, with no model defect at all.

That is a genuine defect with a specific, bounded consequence — not a ceiling.

---

## C2 — "noise averages down with more crystals; differential bias does not"

True only for **iid** noise, which I did not state and the gate assumed
throughout. Real crystal-to-crystal variation for one variant also contains
components that are *correlated* across crystals and do not shrink with n:

* crystal packing and space group
* refinement protocol, software, and deposition era
* construct, expression background, tags, surface mutations
* cryo vs room temperature
* ligand / ion / buffer state

The honest model is variance components:

```
sigma^2_total(n) = sigma^2_iid / n  +  sigma^2_systematic
```

As n → ∞ this goes to `sigma_systematic`, not 0. So the stage-7 claim that
"more crystals per variant buys resolution as 1/sqrt(n)" holds only until the
systematic floor is reached, and the floor's size is **unknown** — it is
empirically estimable only from real redundant crystals, ideally stratified by
crystal form so the between-form component can be separated from the
within-form one.

Consequence for the verdict: **"clears the bar even with one crystal per side"
is a synthetic sensitivity result, not a statement of benchmark readiness.**
Relabelled as such. `perturbation.repeatability()` is documented as measuring the
iid term only.

---

## C3 — Ideal-geometry PASS does not transfer to real PDB data

The gate ran entirely on `helix_on_arc` output: exact regular geometry, optionally
perturbed by iid Gaussian noise. Real helices are irregular — pre-bent, frayed at
the ends, variable in rise and twist, with non-Gaussian and spatially correlated
coordinate error. Nothing in the gate establishes that its slopes, noise figures,
or `tau*` values hold on actual coordinates.

The verdict is therefore restated as:

> **Gate passed on ideal geometry.** The metric construction is sound and the
> v1 metric's SS-dependent offsets and phase-dependent sign are real defects.
> Whether the *label* is usable is not yet established and requires the
> empirical seed test.

### Also mode-dependent

The gate probed one deformation mode: a smooth circular arc distributed over the
whole span. `perturbation.py` adds a second — a rigid kink at the window midpoint,
which preserves every Cα-Cα distance exactly. Normalised to the same half-axis
turning, the two modes give similar sensitivity for v1 (0.7–1.2× ratio across
conformations), so this particular worry did not materialise. But it had to be
checked, and which mode real point mutations actually produce is an open empirical
question that changes what the right metric is.

---

## What still stands

Unaffected by all three corrections, because each is a direct measurement with no
averaging step and no ideal-vs-real extrapolation in it:

1. **v1's metric has large SS-dependent offsets on a straight axis** — 110.4°
   (α), 132.4° (π), 66.8° (3₁₀), 33.7° (PPII), 0.0° (β). It is an SS classifier.
2. **v1's metric couples twist to apparent bend** at 1.026 °/° in α-helices; the
   hybrid at 0.000. (The *ratio* to bend sensitivity was quoted as 79× using the
   cancelled slope; with the corrected 0.316 it is **3.2×** — still the dominant
   systematic, but not the number I published.)
3. **The bisector construction is exact** — 0.000 on a straight axis for every
   helical twist, R² = 1.0000 against ground truth, calibration constant stable
   to 1% across α/π/3₁₀/PPII, no coverage gap over twists 60–185°.
4. **`smooth_pca ≈ smooth_chord` throughout** — smoothing was the missing
   ingredient, the PCA step was incidental.
5. **The algebraic circle fit is rejected** — residual wobble dominates it
   (141–173° on a straight axis).
6. **v1's `_dihedral` is 180° off standard convention** — harmless for tree
   models, but the feature is not what its docstring names.
7. **Every `AUDIT.md` finding** — those are separate measurements on v1's code and
   committed CSV, and none of them route through the gate's slope calculation.

## Next step unchanged, and now better motivated

The empirical seed test, not the full miner. 20–50 curated variants with multiple
WT and mutant crystals; for **both** v1 and v2 metrics measure within-variant
repeatability (separating iid from systematic where crystal form allows), the
observed WT→mutant Δ distribution, τ, the null mover FPR, and the resulting
ceiling. Plus the real-coordinate Jacobian via `perturbation.jacobian`, to test
whether the ideal-geometry slopes transfer at all.

This environment has no access to any structure database (`files.rcsb.org`,
`www.ebi.ac.uk`, `files.wwpdb.org`, `data.pdbj.org` all unreachable) and no cached
coordinates, so the seed test cannot be run here. `perturbation.py` is validated
against ideal geometry and ready to point at real windows.
