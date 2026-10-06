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
