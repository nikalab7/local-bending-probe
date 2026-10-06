# Pre-specification: physics-derived predictions of local backbone response

Exploratory, after the lockbox run. Written and committed **before** any of the quantities below were computed. Only group sizes were counted beforehand (no outcomes), to choose thresholds that give testable groups. Data: dev + first-lockbox rows (1361), frozen label and structures unchanged. Script: `theory.py` → `results/theory.json`.

## Physical picture

A mutation is a local perturbation of the folded structure's free-energy minimum. To first order, the structural response is the protein's compliance times the perturbing force (linear response: Ikeguchi, Ueno & Kidera 2005, *PRL* 94). Two mutation-specific "forces" are well documented:

- **Steric volume change in the core.**
  - Large→small substitutions leave a cavity that the surrounding structure partly collapses into (Eriksson et al. 1992, *Science* 255, 178).
  - Small→large substitutions introduce strain and displace neighbours (Liu, Baase & Matthews 2000, *J. Mol. Biol.* 295, 127).
- **Backbone torsional preference.** Each residue type has its own Ramachandran distribution (Ting et al. 2010, *PLoS Comput. Biol.* 6). A WT (phi, psi) that is improbable for the new residue should be relieved by a local change.

## Measurements (fixed now)

### Radial displacement around the mutated side chain (T1–T3)

1. **Superposition.**
   - Superpose the mutant crystal on the WT scaffold by Kabsch on C-alpha atoms of residues whose WT C-alpha lies > 10 A from the WT C-alpha of r; require >= 30 such atoms.
   - Mutant chain choice as in `mechanisms.py`: residue r carries the mutant residue type, minimum window C-alpha RMSD.
2. **Contact atoms**, chosen from **WT positions only** so that the selection cannot depend on the response: heavy atoms of residues with |i − r| >= 2, present in both structures, within 7.0 A of the WT C-beta of r. WT Gly rows have no C-beta and are excluded.
3. **Displacement.**
   - u = unit vector from the WT C-beta of r to the contact atom's WT position.
   - Radial displacement of an atom = (x_mut − x_wt) · u, in A.
   - The row value is the mean over its contact atoms (>= 3 atoms required); with several mutant crystals, the median.
4. **Null reference.** The same measurement on 1500 WT-vs-WT pseudo-mutant rows (seed 0): held-out WT crystal vs scaffold, with the residue's own WT C-beta.

**Groups** (all exclude Gly or Pro on either side; burial = hse_up >= dataset median; Δvol = volume(mut) − volume(wt) in A^3, the volume table of `delta_model.py`):
- **overpacking:** buried, Δvol >= +25 (n = 104, 55 families);
- **cavity:** buried, Δvol <= −25 (n = 252, 110 families);
- **neutral:** buried, |Δvol| < 10 (n = 194, 94 families);
- **sensitivity:** burial >= 67th percentile (70 / 189 / 130 rows).

### Ramachandran strain (T4)

- **Reference distributions** from the (phi, psi) of every residue in the WT scaffold structures of all proteins (structure only, no labels). Classes: Gly; Pro; pre-Pro (non-Gly, non-Pro residue followed by Pro); Ile/Val; general (all others).
- **Density estimate:** 10-deg bins, circular Gaussian smoothing (sigma = 1 bin), plus a pseudocount of 1e-4 of the total mass.
- **Strain score S** = log P_class(mut)(phi_r, psi_r) − log P_class(wt)(phi_r, psi_r), at the WT dihedrals of r.
  - If the substitution creates or removes a Pro, add the same difference for residue r−1 (its class changes between pre-Pro and non-pre-Pro).
  - Negative S = the WT backbone is less probable for the new residue (strain).
- **Outcome:** dihedral response at r = max(|z_phi|, |z_psi|) (frozen per-metric z), plus the combined mover label.

## Hypotheses and pre-specified tests

All CIs are 95% family bootstrap (2000 resamples). Groups below 30 rows are descriptive only.

| id | prediction | test | supported if |
|---|---|---|---|
| **T1** (primary) | overpacking pushes neighbours **outward** | mean radial displacement, overpacking group; overpacking − null | mean > 0 with CI > 0 **and** difference vs null with CI > 0 |
| **T2** | cavities pull neighbours **inward** | mean radial displacement, cavity group; cavity − null | mean < 0 with CI < 0 **and** difference vs null with CI < 0 |
| **T3** | the response scales with the volume change (dose-response); neutral substitutions look like the null | Spearman(Δvol, row radial) over all buried non-G/P rows; neutral − null | Spearman > 0 with CI > 0 (the neutral − null CI is expected to include 0) |
| **T4** | backbone strain for the new residue → dihedral change | Spearman(−S, max(\|z_phi\|, \|z_psi\|)) over all rows; mover rate (combined) for S <= −2 vs S > −0.5; paired dAUC (site + subst + S vs site + subst; family-out CV as in `final_eval.py`) | Spearman > 0 with CI > 0 **and** the strained − unstrained mover rate with CI > 0 (the dAUC is reported, with no threshold) |

**Conditional analyses** (run only if at least one of T1–T4 is supported):
- **T5, a physics composite with no fitting.**
  - P = z(max(−S, 0)) + z(|Δvol| × hse_up / 20), with z = standardized over all rows.
  - Compliance C = z(window B-factor of the WT, `b_window`).
  - Score = P + C.
  - Reported: AUC for the combined mover label and Spearman with |z| (all rows), and paired dAUC vs the SS-only baseline, both computed directly from the score (no model).
- **T6, structural-class heterogeneity.** For each supported hypothesis, the same statistic within all-alpha (helix fraction >= 0.4, strand < 0.1), all-beta (strand >= 0.3, helix < 0.1) and alpha/beta (rest) proteins, classified by the WT scaffold's HELIX/SHEET composition. Reported with the difference to the pooled estimate.

**What would count against the physical picture:**
- **T1 or T2 in the opposite direction:** neighbours move toward an added side chain, or away from a cavity.
- **No dose-response (T3 fails):** steric volume change would then not be what drives the local response.

---

## Results (appended after the run; the pre-specification above is unchanged)

Run: `python theory.py` → `results/theory.json`. Radial displacement was computed for 1319 of 1361 rows and 1412 of 1500 null rows.

| id | result (95% family bootstrap) | verdict |
|---|---|---|
| **T1** overpacking pushes outward | overpacking mean radial **+0.072 A [+0.048, +0.092]** (n = 104, 55 families); minus null **+0.068 [+0.042, +0.089]**. Sensitivity (top third buried): +0.061 [+0.037, +0.084] (n = 70) | **supported** |
| **T2** cavities pull inward | cavity −0.001 A [−0.018, +0.014] (n = 252, 110 families); minus null −0.005 [−0.025, +0.012]. Sensitivity: −0.010 [−0.036, +0.012] | **not supported**: no measurable collapse |
| **T3** dose-response | Spearman(Δvol, radial) **+0.302 [+0.228, +0.366]** (n = 673 buried non-G/P rows). Sensitivity: +0.304 [+0.222, +0.401]. Neutral − null +0.013 [−0.004, +0.027] (includes 0, as expected) | **supported** |
| **T4** Ramachandran strain | Spearman(−S, max \|z_phi\|, \|z_psi\|) **+0.116 [+0.030, +0.193]** (n = 1361); with the combined \|z\|: +0.121 [+0.050, +0.185]. Mover rate for S <= −2: 0.786 vs 0.320 for S > −0.5, difference +0.466 [+0.328, +0.640], but the strained group has **n = 28 < 30**, so it is descriptive only. dAUC (site + subst + S vs site + subst) +0.006 [−0.005, +0.016] | **partly supported**: the Spearman criterion is met; the group contrast is below the reporting threshold. `theory.py` printed "supported" without applying the n < 30 rule; this table is the verdict |

**Conditional analyses** (run because T1 and T3 are supported):
- **T5, physics composite (no fitting):** AUC 0.573 [0.538, 0.604], against the frozen model's 0.645 (OOF on the same rows); vs SS-only +0.003. As a predictor of the binary mover label, the unfitted composite is weak.
- **T6, structural classes:**
  - all-alpha (482 rows): overpacking − null +0.098 (n = 41), T3 rho +0.39;
  - alpha/beta (824 rows): +0.050 (n = 56), rho +0.26;
  - all-beta (55 rows, overpacking n = 7): descriptive only.
  - The direction and dose-response effects hold in both classes large enough to test. No claim of heterogeneity is made.

**Reading.**
1. **A physical law, confirmed across many families.** A larger side chain in the core pushes the surrounding atoms outward. The effect is about +0.07 A on average over atoms within 7 A, and it grows with the added volume (rho ≈ 0.3). It is systematic, small, and invisible to a WT-vs-WT comparison.
2. **Cavities: no measurable collapse with this measurement** (T2). *Update (`CAVITY.md`): neither this C-beta-based measurement nor a removed-atom-directed one detects collapse even in T4L, the positive control. So T2's null result says nothing about whether removed volume is filled, and no asymmetry is claimed.* Original reading: This is consistent with cavity-creating mutations often leaving the cavity largely open (Eriksson et al. 1992 report variable, partial relaxation). A limitation of the pre-specified measurement: direction is measured from the WT C-beta, and collapse toward the far end of a removed side chain may be under-captured.
3. **Strain in the backbone torsion predicts dihedral change**, weakly over all rows (rho ≈ 0.12). The few strongly strained cases mostly move (22 of 28), but that group is too small for a claim.
4. **Why this barely shows in the mover AUC.** The binary label asks "beyond crystal noise in one window"; these effects are ~0.07 A shifts, averaged over many atoms and directed. The continuous, direction-aware measure recovers physics that the thresholded label cannot.
