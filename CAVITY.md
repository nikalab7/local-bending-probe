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
