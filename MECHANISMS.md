# Pre-specification: known-mechanism label check and backbone geometry

Exploratory, after the lockbox run. It was written and committed **before** any of the numbers below were computed, except the descriptive counts already shown in chat (mover rates of X→Pro in helix, Gly at positive phi, X→Gly, Pro→X on dev). Data: dev rows plus the first-lockbox rows. The lockbox has been opened once, so it is no longer a test set and is used here only to add cases. The frozen label and the pipeline are unchanged. Script: `mechanisms.py` → `results/mechanisms.json`.

**Reporting rule (as in PROTOCOL.md):** groups with < 30 rows or < 10 per class are reported descriptively only, with no claim. All CIs are 95% family bootstrap.

## Part 1 — Does the label see effects that physics says must happen?

Groups are defined from the WT scaffold structure and the substitution only. For each group the readout is the mover rate (combined label and bend-only label) vs all other rows, with the difference and its CI. Expected direction:

| id | group | definition | expected |
|---|---|---|---|
| M1 | X→Pro inside a helix | mut = P, WT residue r and r−1, r−2, r−3 all in HELIX records | higher than rest |
| M2 | Gly at positive phi → X | wt = G, WT phi(r) > 0, mut ≠ G | higher |
| M3 | X→Gly | mut = G | higher |
| M4 | Pro→X | wt = P | higher (weaker) |
| M5 | disulfide Cys→X | wt = C and its SG within 2.5 A of another Cys SG in the WT scaffold | higher |
| M6 | cis peptide at r | omega(r−1, r) in the WT scaffold within 30 deg of 0 | higher |
| M7 | negative mechanism control | BLOSUM62 >= 1, no G/P on either side, hse_up below the dataset median (exposed) | **not** higher than rest |

**Diagnosis of misses.** For rows in M1, M2, M5 and M6 that are *not* movers, report:
- the raw change (|delta| in degrees) and the noise floor (sigma, n_wt), to tell "no real change" from "change hidden by a wide floor";
- the largest per-residue phi/psi change in the window, to catch a change the C-alpha bend misses.

## Part 2 — What moves when the label says "moved"? (geometry)

For every row: superpose the mutant crystal's window C-alphas (residues s..s+4) on the WT scaffold's (Kabsch). Then measure:
- the C-alpha RMSD;
- for each of the 4 peptide planes in the window (i → i+1): the carbonyl O displacement, dpsi(i) and dphi(i+1).

Rows with several mutant crystals take the median of each quantity. The same is computed for a WT-vs-WT reference: 1500 null pseudo-mutant rows drawn with seed 0, comparing the scaffold with the held-out WT crystal.

**Peptide flip (primary definition):** a plane with O displacement > 2.0 A and |dpsi(i)| > 60 deg and |dphi(i+1)| > 60 deg. Sensitivity: O displacement > 1.5 A.

| id | hypothesis | expected |
|---|---|---|
| G1 | the flip rate is higher in "dihedral-only movers" (combined-label movers with \|z_bend\| <= 2) than in bend movers, in non-movers, and in the WT-vs-WT null | dihedral-only > bend movers, non-movers, null |
| G2 | in flip planes the motion is compensated (crankshaft): dpsi(i) ≈ −dphi(i+1) | correlation < 0, slope near −1 |
| G3 | flips are enriched in substitutions involving Gly (wt = G or mut = G) | higher flip rate with Gly than without |
| G4 | C-alpha RMSD of dihedral-only movers is close to the null's, i.e. the C-alpha trace stays put | median C-alpha RMSD of dihedral-only movers within the null's interquartile range |

**What would falsify the "flip" reading** of the dihedral finding:
- the flip rate of dihedral-only movers is not above the null's (G1 fails); or
- flip planes are not compensated (G2 fails).

In either case, the substitution signal on phi/psi labels would be something other than peptide flips.

---

## Results (appended after the run; the pre-specification above is unchanged)

Run: `python mechanisms.py` → `results/mechanisms.json`. 1361 rows (dev 1250, first lockbox 111), mover rate 0.35. Geometry was computed for all 1361 rows and 1500 null rows.

### Part 1

| id | n (families) | mover rate, combined | rest | difference [95% CI] | bend label | status |
|---|---|---|---|---|---|---|
| M1 X→Pro inside a helix | 3 (3) | 0.67 | 0.35 | +0.32 [−0.36, +0.68] | 0.33 vs 0.27 | descriptive only |
| M2 Gly at positive phi → X | 16 (7) | 0.62 | 0.35 | +0.28 [+0.01, +0.66] | 0.38 vs 0.27 | descriptive only |
| M3 X→Gly | 60 (35) | **0.60** | 0.34 | **+0.26 [+0.10, +0.38]** | 0.40 vs 0.27 | as expected |
| M4 Pro→X | 16 (10) | 0.56 | 0.35 | +0.21 [−0.05, +0.49] | 0.25 vs 0.27 | descriptive only |
| M5 disulfide Cys→X | 2 (2) | 0.50 | 0.35 | — | — | descriptive only |
| M6 cis peptide at r | 6 (3) | 0.67 | 0.35 | +0.32 [+0.14, +0.67] | 0.33 vs 0.27 | descriptive only |
| M7 conservative, exposed, no G/P (control) | 200 (91) | **0.22** | 0.37 | **−0.15 [−0.22, −0.08]** | 0.17 vs 0.29 | as expected (lower) |

- The label behaves physically sensibly where n allows a test. X→Gly moves more often; conservative exposed substitutions move less often.
- The "physics guarantees it" groups are mostly too small for any claim (3–16 rows).
- **Misses of strong mechanisms (M1, M2, M5, M6; 10 rows).**
  - **8 are real small changes:** the largest per-residue dihedral change is 2–29 deg and the C-alpha RMSD is <= 0.19 A. Examples: Gly→Ala at positive phi in human lysozyme (G105A, G127A: ~10 deg). The physical expectation is not 100%. Ala tolerates the left-handed region at some strain.
  - **2 look hidden by a wide noise floor:** human lysozyme G72A (sigma 12 deg) and O53512 G190P (sigma 5.9 deg, a 165-deg dihedral change, C-alpha RMSD 0.69 A, yet z < 1).
- A 100%-certain effect does not exist in these data. The closest ones occur in a few rows each.

### Part 2

| id | result | vs expectation |
|---|---|---|
| G1 | flip rate: dihedral-only movers 0.031 [0.005, 0.064] (n = 196); bend movers 0.038; non-movers 0.006; null 0.007. Dihedral-only − null = +0.024 [−0.003, +0.054]; − bend movers = −0.007 [−0.041, +0.026] | **not supported**: not above bend movers; above the null only at a CI touching 0 |
| G2 | flip planes of movers (21 planes): corr(dpsi, dphi) −0.41, slope −0.38; null (10 planes): +0.14 | partly compensated, slope far from −1: **not supported** |
| G3 | flip rate with Gly involved 0.091 vs 0.012 without; diff +0.079 [+0.011, +0.139] (n 99 vs 1262) | **supported** |
| G4 | median C-alpha RMSD: dihedral-only movers 0.095 A; bend movers 0.125; non-movers 0.072; null 0.041 (IQR 0.022–0.066) | **not supported**: the dihedral-only movers' C-alpha trace moves more than the null's |

**Reading.**
- **The flip hypothesis is falsified as the main explanation.** Peptide flips are rare (~3% of dihedral-only movers), are not specific to them, and are only partly compensated. Flips happen mainly where Gly is involved.
- **Dihedral-only movers are small, distributed changes.** Their median C-alpha RMSD is 0.1 A, at or below typical coordinate error at ~2 A resolution.
- **New, unplanned observation** (to be tested, not claimed):
  - Non-movers (0.072 A) and conservative exposed substitutions (22% movers) both sit well above the WT-vs-WT null (0.041 A; 2.6% false positives).
  - So mutant-vs-WT comparisons carry an excess difference that WT-vs-WT comparisons do not.
  - Candidates: a real small perturbation from any mutation; or experiment-level differences (different laboratory, refinement program or restraints, deposition era) that the WT-vs-WT null, built mostly from WT series of the same groups, does not sample.
