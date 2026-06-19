# proteinX — Bending Metric & WT/Mutant Pair Benchmark (spec v1)

Gating artifact for the "win where AlphaFold is weak" plan. Until this exists and
is consistent, the headline comparison (proteinX vs AF2 on real mutation-induced
geometry change) is undefined.

## 1. The single bending metric

Canonical implementation: [`bending_metric.py`](bending_metric.py). **Import it; never
reimplement.** Self-test verified: straight = 0.00°, right-angle = 90.00°.

- Window = 5 C-alpha atoms, residues `i .. i+4` (the v2 outcome region).
- Split at midpoint `i+2`; fit a principal axis to each overlapping half
  `{i,i+1,i+2}` and `{i+2,i+3,i+4}`; bending = angle between axes, degrees [0,180].
- **Intrinsic** (C-alpha internal geometry only) ⇒ rotation/translation invariant ⇒
  `delta_bending(wt, mut)` needs **no superposition**. This is the key correctness
  property: it kills alignment-induced artifacts and guarantees the metric is
  identical on experimental and AF-predicted coordinates.
- Used identically by all three consumers: training labels, ground-truth deltas,
  AF2 baseline.

**Noise floor (set empirically, not assumed):** the spread of `bending` across
redundant crystal structures of the *same* variant is the detection floor. A
`delta_bending` smaller than that spread is not real. The same observed-Δ
distribution supplies the quartile buckets for the contrastive explanation layer —
this is how the old a-priori >30°/<10° thresholds get replaced by data.

## 2. WT/mutant pair benchmark — procedure

Goal: tuples `(WT structure, mutant structure, mutation list, residue
correspondence, observed delta_bending, flags)` where both structures are
experimental, so we have ground truth for what a mutation *actually did* to local
geometry.

### 2.1 Source candidate pairs
- **Primary (systematic):** RCSB sequence clusters + SIFTS UniProt mapping → within
  a same-protein cluster, find chains differing by 1–N point substitutions.
- **Seed/sanity:** classic mutational series (T4 lysozyme, SNase, barnase, ubiquitin,
  lac repressor) — high-quality, many matched WT/mutant deposits.

### 2.2 "Matched pair" definition (hard rules)
- Same protein/region; sequence differs only by **substitutions**, `1 ≤ N ≤ 4`.
- **No indels** — they break residue registration. (This rule is what lets window
  `i..i+4` map to the same residues in WT and mutant for free.)
- Otherwise identical over the aligned region.

### 2.3 Confound controls (most false pairs die here)
| Confound | Control |
|---|---|
| Resolution mismatch | both ≤ project threshold (≤2.0–2.5 Å); record the gap |
| Crystal packing | flag windows whose outcome residues sit in a crystal contact |
| Ligand/cofactor/ion mismatch near window | require matched local ligand state, else flag |
| pH / temperature / cryo vs RT | read from header; flag large differences |
| Distal (allosteric) mutations | benchmark focuses on windows **containing** the mutation site; distal effects are explicitly out of local scope |
| Long-range tertiary contact in window | flag → inference **abstains** here (don't exclude at eval; report coverage) |

### 2.4 Registration & aggregation
- Residue correspondence comes straight from the (indel-free) sequence alignment.
- Bending is intrinsic ⇒ **no superposition needed** for the metric. Superposition is
  used only for crystal-contact detection and visualization.
- Many proteins have several structures per variant: aggregate to **one bending per
  variant** (median across redundant structures); keep the spread as the per-variant
  noise estimate (feeds the noise floor in §1).
- `delta_bending = bending(mutant_variant) − bending(WT_variant)` per (mutation, window).

### 2.5 Leakage control
- Training is sequence-culled (≤30–40%). **Hold out whole sequence clusters** for the
  benchmark; no family appears in both train and validation.

### 2.6 Output schema (one row per mutation×window)
`pdb_wt, chain_wt, pdb_mut, chain_mut, mutations[(wt_aa, seqpos, mut_aa)],
window_start..window_end, bending_wt, bending_mut, delta_bending,
within_variant_noise, flags{crystal_contact, ligand_mismatch, resolution_gap,
long_range_contact}, coverage_class{local_reliable | abstain}`

## 3. AlphaFold baseline — fairness rule
Predict **both** WT and mutant with identical AF settings; baseline Δ =
`bending(AF_mut) − bending(AF_wt)`, same metric, same window. Do **not** compare
AF-mutant against the crystal WT — that confounds AF's prediction error with the
mutation effect. Expectation: `|Δ_AF| ≈ 0` (documented point-mutation insensitivity),
so even modest correlation from proteinX is a win. Round out baselines with ESMFold,
an I-sites/fragment lookup, and a naive Pro/Gly heuristic (the floor to beat).

## 4. Tooling recommendation
- Structure I/O + assemblies + altloc/occupancy: **gemmi** (fast, robust mmCIF).
  BioPython acceptable but slower/looser on mmCIF.
- Altloc policy: highest occupancy (ties → 'A'); X-ray only for the metric (resolution
  constraint already implies this); model 1 for any NMR seeds.

## Open decisions before building the miner
1. Pair source: RCSB clusters+SIFTS (systematic) vs curated seed first? (Recommend:
   seed set first to validate the whole pipeline end-to-end on ~tens of clean pairs,
   then scale to systematic.)
2. `N` cap on simultaneous mutations (recommend N=1 for the first benchmark — cleanest
   attribution — then relax).
