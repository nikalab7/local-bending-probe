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
