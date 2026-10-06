# Pre-specification: is the mutant-vs-WT excess experimental (laboratory / refinement)?

Exploratory, after the lockbox run. Committed **before** any laboratory, refinement-program or deposition metadata was fetched, and before any PDB-REDO model was downloaded for analysis. A single PDB-REDO file (3KGF) was downloaded only to check its format.

**Motivation** (`MECHANISMS.md` results):
- Conservative, exposed substitutions are labelled movers in 22% of cases, and C-alpha RMSDs of non-movers are 0.072 A.
- The WT-vs-WT null gives 2.6% false positives and 0.041 A.
- Question: is this excess caused by the mutation or by the experiment?

Data: dev and first-lockbox proteins (the lockbox is spent and is used only to add cases), with the frozen label (`PROTOCOL.md`), unchanged. Script: `labtest.py` → `results/labtest.json`, `results/labtest_redo.json`.

## 0. Sanity check: circular arithmetic (done before this file)

- **O53512 G190P (165-deg dihedral change, z < 1).** This is not a wrap error. The WT is bimodal at this window: Gly190 has phi ≈ +135 deg in 5 of 11 WT crystals and phi ≈ −65 deg in 6. The Pro mutant sits in the −65 state, which WT already visits. The 165-deg figure in `MECHANISMS.md` compared the mutant with one scaffold crystal (3KGF) that happens to be in the other state. The wide sigma is the WT bimodality, and the label correctly calls this "not beyond WT variation". It is a population shift, which the label by design does not count as movement.
- **Code check.** Real rows and null rows both go through `form_stats` (circular mean reference, `wrap(v − ref)`) → `rel` → `label`, and NCS averaging wraps as well. A unit test (`test_circular_metric_wraps_across_180`) asserts that WT values straddling ±180 deg give a small sigma and that a mutant at −177 deg gives a small delta.

## 1. Definitions (fixed now)

**Entry metadata (RCSB):** `audit_author.name`, `software` entries with `classification = "refinement"`, `rcsb_accession_info.deposit_date`, `pdbx_database_status.status_code_sf`.

- **Author key:** personal names only (`"Last, F.M."` → lowercase `last,f`). Names without a comma, or containing *project, consortium, center/centre, initiative, genomics, institute, laboratory, university*, are dropped, so a consortium label cannot link unrelated groups.
- **Same lab (two entries):** share at least one author key.
- **Lab overlap f** of an entry e with a reference WT set R: the fraction of R's crystals that share a lab with e.
  - **same-lab:** f >= 0.5
  - **cross-lab:** f = 0
  - **mixed:** 0 < f < 0.5, reported descriptively only.
- **Refinement program family:** the first refinement software name, upper-cased and mapped by prefix: REFMAC, PHENIX, BUSTER, CNS, X-PLOR/XPLOR, SHELX, TNT, PROLSQ, RESTRAIN, else OTHER. The relation of e to R is the fraction of R with the same family. **Same program:** >= 0.5; **different:** 0.
- **Deposition-year gap:** |year(e) − median year(R)|, in bins 0–2, 3–9 and >= 10 years.
- **Condition differences:** Δres = res(e) − median res(R); ΔT = |T(e) − median T(R)|, from the PDB headers parsed by the frozen code.
- **Matched conditions:** |Δres| <= 0.3 A and ΔT <= 50 K, both known.
- **The entity e:**
  - for a null pseudo-mutant: the held-out WT crystal; R = the other WT crystals of the form;
  - for a real row: its mutant crystal(s), with the authors as the union over crystals; R = the WT crystals of the form. With several mutant crystals, years, resolution and temperature are their medians.

## 2. Labels for this test

- `pairs.load_proteins` / `pairs.load_lockbox` are re-run with `pairs.NULL_HELDOUT = 40` (all WT crystals of a form, up to the 40-crystal cap, are held out in turn), instead of the 10 best-resolution ones. The aim is enough cross-lab WT–WT cases.
- The reference WT list per (protein, form) is captured by wrapping `pairs.null_rows` at run time. No frozen file is edited.
- **Check:** real rows must be identical to the frozen-label rows (same mutations, same z). If not, the run stops.

## 3. Primary, decisive test: cross-lab vs same-lab WT–WT null

No mutation is involved, so there is no confound with mutation type.

- **Estimate:** the mover rate (frozen label, |z| > 2) of null pseudo-mutants in cross-lab vs same-lab classes.
- **Within-protein stratified difference:** only (protein, form) strata that have both classes. Strata are weighted by n_cross × n_same / (n_cross + n_same), a CMH-type risk difference.
- **CIs:** 95% family bootstrap (2000 resamples), with weights recomputed in each resample.
- **Reported:** pooled class rates, the stratified difference, and the same difference under matched conditions (the adjusted estimate).

**Decision rule** (on the adjusted, matched-condition, within-protein estimate):
- **experimental excess supported:** difference >= +0.05 with CI excluding 0, **and** the cross-lab rate >= 0.10;
- **rejected:** CI upper bound of the difference < +0.05 **and** cross-lab rate < 0.10;
- otherwise **inconclusive**.

For reference: real mutant rows move at 0.35 overall and 0.22 for conservative exposed substitutions.

## 4. Secondary: same-lab vs cross-lab mutant–WT rows, within protein

- Real rows are classified by the lab overlap of their mutant crystal(s) with R.
- **Estimate:** within-protein stratified difference (cross − same), in strata (protein, form) with both classes, with family-bootstrap CI. Raw and matched-condition versions.
- Also, with CIs (descriptive):
  - the same comparison for refinement program (same vs different);
  - mover rates by deposition-year gap bin;
  - the **mutational excess within each lab class**: real mover rate minus null mover rate, same-lab and cross-lab separately.
- Strata with fewer than 30 rows in total contribute to the pooled estimate but are not reported on their own.

## 5. PDB-REDO re-check (uniform re-refinement)

**Subset**
- Manifest proteins (mined dev + first lockbox), excluding T4L and proteins with > 150 manifest entries (cost).
- PDB-REDO must have a model for >= 90% of the protein's manifest entries.
- Random sample of 40 such proteins, seed 0.
- If fewer than 40 qualify, all qualifying proteins are used, and the count is reported.

**Files**
- Hybrid files: the original PDB header (all records except ATOM/HETATM/ANISOU/TER/MODEL/ENDMDL/CONECT/MASTER/END), followed by the PDB-REDO `*_final.pdb` coordinate records.
- The resolution, temperature, HELIX/SHEET and crystal records the frozen parser reads are therefore unchanged; only coordinates and B-factors differ.
- Entries without a PDB-REDO model keep their original coordinates, and the share of such entries is reported.

**Readout**
- Labels are rebuilt with the frozen code (NULL_HELDOUT = 40).
- Rows present in both builds (same protein and mutation) are compared paired: mover rate, original vs PDB-REDO, difference with family-bootstrap CI.
- The null false-positive rate is compared the same way.
- The excess (real − null) is compared before and after.

**Decision rule**
- **Re-refinement supports the lab/refinement explanation** if the paired mover rate drops by >= 0.05 with CI excluding 0 **and** the excess (real − null) shrinks.
- **It does not** if the CI of the paired difference excludes −0.05.
- Otherwise **inconclusive**.

## 6. What each outcome would mean

- **Primary supported:**
  - WT crystals from different labs differ as much as many "mutations" do.
  - A WT-vs-WT null drawn from one lab's series underestimates the noise of real mutant-vs-WT comparisons.
  - The label then needs a cross-lab-aware floor, or same-lab pairs only.
- **Primary rejected:** lab identity does not explain the excess, so the excess is more likely mutational (small real perturbations), or comes from something not measured here.
- **PDB-REDO:** tells whether the lab effect, if any, comes from refinement choices (removable by uniform re-refinement) or from the crystals and data themselves.
