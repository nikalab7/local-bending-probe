# Pre-specification: cavity follow-up to THEORY.md T2

Exploratory. Written and committed **before** any of the quantities below were computed. Only group sizes were counted beforehand: T4L buried cavity rows = 36 at 26 sites; non-T4L = 216 rows in 109 families; T4L WT–WT null rows = 532.

`THEORY.md` T2 found no inward collapse of neighbours around large→small core substitutions: −0.001 A [−0.018, +0.014]. The direction was measured from the WT C-beta of r. Two readings are possible:
- **(a)** the measurement misses the collapse, because the cavity lies beyond C-beta, toward the removed atoms;
- **(b)** removed volume is genuinely not filled.

T4L is the positive control: cavity relaxation there is documented (Eriksson et al. 1992). Script: `cavity.py` → `results/cavity.json`.

## Fixed definitions

- **Cavity group:** buried (hse_up >= the median over all 1361 rows), Δvol <= −25 A^3, no Gly or Pro on either side (as in `THEORY.md`).
- **Superposition, mutant chain choice and multi-crystal median:** as in `theory.py` (C-alpha atoms > 10 A from r; >= 30 atoms).
- **Null rows:**
  - T4L steps: all T4L WT–WT pseudo-mutant rows (dev + lockbox caches);
  - non-T4L step: 1500 non-T4L null rows drawn with seed 0.
  - Null rows carry their real row's wt/mut, so "removed atoms" are defined for them too, and their expected value is 0.
- **CIs:**
  - T4L (one family): bootstrap over residue sites (r), 2000 resamples;
  - non-T4L: family bootstrap, 2000 resamples.
  - All 95%.

**Measurement A (C-beta direction, as in THEORY.md):** unchanged `theory.radial_one`. Negative = toward C-beta (collapse).

**Measurement B (removed-atom direction, new):**
1. Removed atoms = side-chain atoms of WT residue r whose names do not exist in the mutant residue type (e.g. L→A removes CG, CD1, CD2). These are taken from the WT scaffold; rows with none present are excluded.
2. c_rem = centroid of the removed atoms (WT positions).
3. Contact atoms = heavy atoms of residues with |i − r| >= 2, present in both structures, within 6.0 A of any removed atom (WT positions only).
4. For each contact atom, v = unit vector from the atom's WT position toward c_rem. Inward displacement = (x_mut − x_wt) · v; **positive = movement into the cavity.**
5. Row value = mean over >= 3 contact atoms; with several mutant crystals, the median.

## Steps and decision rules

**Step 1: T4L, measurement A.**
- Mean radial displacement of the T4L cavity rows, and the difference vs the T4L null.
- **"T4L shows collapse with A"** if the mean is < 0 with CI < 0 **and** the difference vs null is < 0 with CI < 0.

**Step 2 (only if Step 1 shows no collapse): validate measurement B on T4L.**
- Mean inward displacement of the T4L cavity rows, and the difference vs the T4L null (same measurement on null rows).
- **"B validated"** if the mean is > 0 with CI > 0 **and** the difference vs null is > 0 with CI > 0.
- **If B is validated:** run B **once** on the non-T4L cavity rows (difference vs the non-T4L null, family bootstrap).
  - non-T4L collapse: mean > 0 with CI > 0 and difference vs null with CI > 0;
  - otherwise no collapse outside T4L.
- **If B is not validated on T4L:** stop. Report that the cavity response could not be measured with a validated method, and do not run B on non-T4L proteins.

**Step 3 (only if Step 1 shows collapse in T4L): measurement A on non-T4L.**
- The non-T4L cavity rows with measurement A (difference vs the non-T4L null, family bootstrap).
- If non-T4L shows no collapse, the asymmetry below is the main finding.

**Asymmetry**, reported in every branch where the cavity effect on non-T4L proteins is measured (with A in Step 3, or with B after a validated Step 2):
- **Index** = (overpacking − null) + (cavity − null), with **both terms on the C-beta-based radial scale** (measurement A), non-T4L only, family bootstrap.
- **Reading:**
  - symmetric response (added volume pushed out, removed volume filled equally) → index ≈ 0;
  - "added volume is pushed out, removed volume is not filled" → index > 0.
- The overpacking term uses the `THEORY.md` definition restricted to non-T4L.
- If B was used, the cavity-side inward displacement (B) is reported next to the index, but **not** added into it, since the scales differ.

**What each outcome means**
- **Step 1 shows collapse:** measurement A works and the T2 null is real outside T4L, or T4L is special (Step 3 decides).
- **Step 1 shows no collapse, B validated, non-T4L collapse:** the original measurement missed the collapse. T2 holds with the right direction.
- **Step 1 shows no collapse, B validated, no non-T4L collapse:** asymmetry, now measured with a validated method.
- **B not validated:** inconclusive. The cavity side cannot be measured reliably here.

---

## Results (appended after the run; the pre-specification above is unchanged)

Run: `python cavity.py` → `results/cavity.json`. T4L: 36 cavity rows at 26 sites; 532 T4L null rows (434 with measurement B defined).

| step | result (95% site bootstrap) | decision |
|---|---|---|
| 1. T4L, measurement A (negative = collapse) | mean **+0.030 A [−0.004, +0.060]**; null +0.002; difference +0.028 [−0.005, +0.058] | **no collapse**, so Step 2 runs |
| 2. T4L, measurement B (positive = into the cavity) | mean **−0.008 A [−0.044, +0.031]**; null +0.006; difference −0.014 [−0.050, +0.024] | **B not validated** → **stop** |

- **Outcome, per the pre-registered rule:** inconclusive. The cavity response cannot be measured here with a validated method. Neither measurement detects collapse in T4L, the positive control. B was therefore **not** run on non-T4L proteins, and the asymmetry index was not computed: the rule computes it only in branches where the non-T4L cavity effect is measured.
- **What this means for THEORY.md T2:**
  - The "no measurable collapse" result stands as a statement about these two measurements only.
  - It is **not** evidence that removed volume goes unfilled, because the measurements failed the positive control.
  - The "added volume is pushed out, removed volume is not filled" asymmetry is **not established**.
- **Possible reasons** (not tested):
  - Relaxation into a cavity may be carried by a few atoms moving a lot. A mean over all atoms within 6 A dilutes it.
  - 36 rows at 26 sites give a CI half-width of about 0.035 A, which may be larger than the mean effect.
  - Any further measurement (e.g. the largest displacements per cavity) would have to be pre-registered and validated first against independently published T4L cavity displacements, not against these data.

## Step 3 (filled-fraction via cavity volumes): not run — final stop

- **What was proposed.** A direct filling measurement: filled fraction = 1 − V_real / V_virtual. It would be validated first against published T4L cavity volumes (Xu, Baase, Baldwin & Matthews 1998, *Protein Sci.* 7:158; Eriksson et al. 1992).
- **Why it could not be validated.** The per-mutant volume table is reachable only as a scanned PDF (PMC2143816). From this environment that PDF sits behind bot protection (PMC proof-of-work; Europe PMC Cloudflare challenge), which was not circumvented. Openly reachable sources give only two anchors: L99A ≈ 150 A^3 (41 A^3 pre-existing in WT; Baase et al. 2010, PMC2867005) and I29A ≈ 0 (Xu et al. 1998 abstract). Two points are not enough for a pre-specified agreement criterion.
- **Decision (user):** stop and summarize. **T2 is reported as "not measurable with these data". No further cavity measurements will be made in this project.** Nothing was computed for step 3.
