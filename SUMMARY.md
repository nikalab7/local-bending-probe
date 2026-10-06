# Mutation-induced local backbone change in crystal structures: calibrated labels, a pre-registered test, and a bounded negative result

*Status: the pre-registered test was inconclusive (underpowered). Every result below except the lockbox numbers is exploratory, because the development data were used for all decisions.*

## Summary

Do single amino-acid substitutions change the local backbone (a 5-residue window), and can that be predicted from the WT structure and the substitution alone? I built WT/mutant crystal-pair labels with an explicit, calibrated noise floor and trained a simple model. I then froze the whole pipeline (`PROTOCOL.md`) and scored it once on a lockbox of 90 unseen proteins.

- **Pre-registered result.** Lockbox AUC was 0.532 [0.403, 0.661] (95% family bootstrap; n = 111 mutations, 89 families). The model did not beat the secondary-structure-only or burial-only baselines. The lockbox could detect only AUC >= 0.685, so the test is **inconclusive**, not negative.
- **Exploratory (dev, 1,250 mutations, 250 families):**
  - The model reached AUC 0.662 [0.633, 0.698] (within-protein 0.653), against a label-reliability ceiling of 0.787. It beat the SS-only baseline by +0.087 and the burial-only baseline by +0.112.
  - **"What" vs "where" depends on the label.** For C-alpha bending, adding substitution features to site features gave no detectable gain: +0.004 [−0.015, +0.025]. For backbone dihedral changes (phi/psi), substitution features added +0.03 to +0.06.
- **Artifacts, quantified in this dataset.** Pairs with a room-temperature mutant and a cryo-cooled WT were labelled "movers" in 57% of cases, against 34% for the rest. Pairs whose mutant is > 0.5 A worse in resolution: 55% vs 34%. That such mismatches create apparent structural differences is documented in prior work; these numbers measure how much they contaminate a mutant-pair label.

## 1. Data and labels

- **Pairs.** X-ray structures <= 2.5 A from the PDB, mined systematically (RCSB search + SIFTS). A mutant is compared only with WT crystals of the same crystal form and the same ligand state near the window, with one row per mutation.
- **Noise floor.** Per window, from >= 3 WT crystals: sigma = SD/c4, shrunk toward a B-factor-conditioned prior, Student-t calibrated; NCS copies averaged.
- **Label.** Max |z| over C-alpha bend, phi, psi and C-alpha torsion, with the threshold calibrated so that WT-vs-WT pseudo-mutants give ~2.6% false positives.
- **Reliability.** Split-half test-retest gives kappa 0.51 and an oracle AUC of 0.787. This is the ceiling for any predictor.
- **Families.** >= 30% sequence-identity clusters. All evaluation holds out whole families.

## 2. Model and evaluation

- **Model.** Logistic regression on WT-only site features (SS, B-factor, burial, non-local contacts, WT window bend, terminus distance) plus physicochemical substitution descriptors, confidence-weighted.
- **Inputs.** Anything derived from the mutant structure or the mutant experiment is forbidden as input, and a test asserts it. Temperature, resolution, altlocs and lattice contacts are diagnostics only.
- **Statistics.** Family bootstrap; paired tests for every comparison; positive and negative controls; a power analysis.

## 3. Results

### 3.1 Pre-registered lockbox (`results/final_lockbox.json`)

| | lockbox | dev (exploratory) |
|---|---|---|
| model AUC | **0.532 [0.403, 0.661]** | 0.662 [0.633, 0.698] |
| model − SS-only | +0.002 [−0.099, +0.117] | +0.087 [+0.052, +0.120] |
| model − burial-only | +0.083 [−0.038, +0.204] | +0.112 [+0.076, +0.160] |
| within-protein | 0.375 [0.00, 0.83] | 0.653 [0.635, 0.711] |
| positive control (helix) | 0.951 | 0.915 |
| negative control (within-family shuffle) | 0.538 (degenerate, see 4c) | 0.529 |
| minimum detectable AUC | 0.685 | 0.545 |
| WT-vs-WT null false positives | 4.7% | 2.6% |

### 3.2 How much does mutation identity add? (`results/posthoc_dev.json`; dev, after the lockbox run, exploratory)

All values are paired dAUC with 95% family-bootstrap CIs. The upper bound is the largest gain the dev data are compatible with.

| addition | over | dAUC [95% CI] |
|---|---|---|
| physicochemical substitution features, **bend label** | site | **+0.004 [−0.015, +0.025]** |
| physicochemical substitution features, combined label | site | +0.059 [+0.033, +0.089] |
| same, phi label / psi label / C-alpha torsion label | site | +0.041 / +0.043 / +0.027 (all CIs > 0) |
| same, combined label, rows without Gly or Pro on either side | site | +0.031 [+0.010, +0.054] |
| ESM-2 substitution log-likelihood ratio | site + subst | +0.001 [−0.007, +0.008] |
| substitution x local-structure terms | site + subst | −0.000 [−0.012, +0.011] |
| *site descriptors (not mutation identity):* ESM-2 site terms / ESM-2 embeddings / elastic network / full-atom context | site + subst | +0.001 / −0.006 / −0.003 / −0.005 (upper bounds +0.011 / +0.009 / +0.001 / +0.010) |

**Bounded claim.**
- **C-alpha bending:** adding mutation-identity features to site features gave no detectable gain; the upper CI bound is **+0.025 AUC** on dev.
- **Backbone dihedral changes:** the physicochemical substitution descriptors do add signal, including outside Gly/Pro substitutions. The evolutionary (ESM-2) and context-interaction descriptors add nothing detectable beyond them; upper bounds are +0.008 and +0.011.
- **Why this negative result is weaker than it looks.** Label noise (kappa 0.37 for the bend label, 0.51 for the combined label) attenuates every effect toward zero. As a crude scaling, the oracle margin is 0.54 of a perfect predictor's for the bend label. On that scale, the bend upper bound of +0.025 would correspond to roughly +0.05 on noise-free labels. "No detectable gain" is therefore not "no effect".
- **A hypothesis, not tested here.** Substitutions might change local dihedrals without moving the C-alpha trace, as peptide-plane flips do. That would explain why identity matters for phi/psi but not for bending.

### 3.3 Crystallographic mismatch in the label (dev; 95% family-bootstrap CIs)

| group | n (families) | mover rate | rest | difference |
|---|---|---|---|---|
| mutant room temperature (> 250 K), WT cryo (< 150 K) | 115 (20) | 0.565 [0.477, 0.684] | 0.337 [0.298, 0.373] | +0.229 [+0.139, +0.347] |
| temperature mismatch > 50 K | 240 (45) | 0.467 [0.391, 0.585] | 0.332 [0.290, 0.372] | +0.135 [+0.053, +0.253] |
| mutant resolution worse by > 0.5 A | 101 (33) | 0.545 [0.456, 0.646] | 0.341 [0.302, 0.379] | +0.203 [+0.106, +0.311] |
| mutant resolution worse by > 0.3 A | 268 (70) | 0.459 [0.399, 0.546] | 0.330 [0.285, 0.370] | +0.129 [+0.062, +0.221] |

**Not new.** These effects are expected from prior work:
- Cryo-cooling remodels conformational distributions (Fraser et al. 2011: > 35% of side chains across 30 proteins; Halle 2004; Keedy et al. 2015, 2018 with multi-temperature series).
- Cooling repacks the lattice (Juers & Matthews 2001).
- Coordinate error grows with resolution (Cruickshank 1999).

**What this dataset adds** is the size of the effect on a WT/mutant "mover" label: mismatched pairs are labelled movers 1.6–1.7x as often.

**Caveats.**
- The room-temperature group comes from only 20 families, so this is a between-group association, not a within-protein estimate.
- Excluding all mismatch-suspect rows *raises* dev AUC (0.662 → 0.701). These rows add label noise; they are not the source of the site signal.

## 4. Methods lessons

a. **The dev-calibrated null did not transfer.**
   - The WT-vs-WT false-positive rate was 2.6% on dev and 4.7% on the lockbox, because label reliability depends on the number of WT crystals per form (median 7 vs 3).
   - A threshold calibrated on well-studied proteins is too permissive for sparsely studied ones.
   - Calibrate, or at least report, the null within strata of WT-crystal count for every new dataset.

b. **Compute power on the post-QC sample before unblinding.**
   - The lockbox shrank from 256 single-substitution entities to 111 labelled mutations after QC, and only then did it turn out to be underpowered.
   - The pre-registered size rule (>= 30 rows, >= 10 per class, >= 10 families) was met but was far too weak.
   - The follow-up (`TIMELOCK.md`) gates unblinding on a blinded, post-QC power calculation.

c. **Negative controls must match the data's structure.**
   - Shuffling labels within families is degenerate when most families have one row (73 of 89 in the lockbox).
   - On dev the same control gave 0.53, not 0.50, because the model partly ranks families by base rate. Within-protein AUC is the cleaner estimate of the site-level signal.

d. **The label definition decides the "what vs where" answer.** Bend-only and dihedral-including labels support different conclusions about substitution identity. Pre-specify the label and report the component labels.

## 5. Limitations

- One crystal-pair data source; 85% of dev labels (95% of lockbox labels) rest on a single mutant crystal.
- Label noise caps achievable AUC at about 0.79.
- Dev decisions were not cross-validated as a procedure (label, feature and model choice), so dev numbers carry selection optimism.
- The lockbox proteins differ systematically from dev: they are less studied and have fewer crystals.
- No comparison with AlphaFold-type predictors on the same pairs (it needs a GPU).

## 6. Pre-registered follow-up

`TIMELOCK.md` registers a prospective test on entries released after 2026-10-05 with the same frozen pipeline, a pinned environment (`requirements.lock`, `Dockerfile`, pre-registration commit) and a blinded power gate.

Feasibility, stated plainly:
- About 9 new post-QC families appear per year.
- About 190 families are needed to detect AUC 0.62, which is roughly two decades.
- A confirmatory test of the site signal on new proteins is therefore **not realistic within ~2 years** from the PDB alone.

## References

- Cruickshank DWJ (1999). Remarks about protein structure precision. *Acta Cryst.* D55, 583–601.
- Fraser JS, van den Bedem H, Samelson AJ, Lang PT, Holton JM, Echols N, Alber T (2011). Accessing protein conformational ensembles using room-temperature X-ray crystallography. *PNAS* 108, 16247–16252.
- Halle B (2004). Biomolecular cryocrystallography: structural changes during flash-cooling. *PNAS* 101, 4793–4798.
- Juers DH, Matthews BW (2001). Reversible lattice repacking illustrates the temperature dependence of macromolecular interactions. *J. Mol. Biol.* 311, 851–862.
- Keedy DA et al. (2015). Mapping the conformational landscape of a dynamic enzyme by multitemperature and XFEL crystallography. *eLife* 4, e07574.
- Keedy DA et al. (2018). An expanded allosteric network in PTP1B by multitemperature crystallography, fragment screening, and covalent tethering. *eLife* 7, e36307.
