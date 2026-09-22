# Independent audit — `local-bending-probe`

Scope: all seven scripts, `SPEC_bending_and_pairs.md`, `RESULTS.md`, `README.md`,
and the committed `feasibility_t4l.csv`. Focus is on mathematical and statistical
correctness, not style.

Every quantitative statement below is reproduced by **`audit_checks.py`** (offline,
deterministic, no PDB downloads). Check IDs in brackets — e.g. `[A1]` — name the
function that produces the number.

```
python3 audit_checks.py          # everything
python3 audit_checks.py A1 B2    # one or more checks
```

A note on direction, because it is easy to misread this audit as hostile to the
project's conclusion. It is not. The stated conclusion is negative ("local sequence
does not predict mutation-induced bending"), and the findings below cut against the
evidence on *both* sides:

- **F1, F2, F4 and F8 undercut the project's positive claims** — Gate 1's "the
  signal is real" and Gate 5's "3D context helps in cores". Gate 1 turns out to be
  the weakest link in the chain, not the foundation.
- **F3 makes the negative result *less* certain than reported, not more** — the true
  confidence intervals are wider than published.
- **F8's positive control comes out in the project's favour** — the architecture
  does recover a planted local effect, which rules out the most obvious way the
  negative result could have been an artifact. That check is missing from the repo
  and is worth adding.

The single summary sentence: most gates measure something other than what they are
captioned as measuring, and the one quantity that would settle the whole project —
the achievable AUC ceiling given the measurement noise (F8) — is never computed.

---

## Severity summary

| # | Finding | Affects | Severity |
|---|---|---|---|
| F1 | The metric is a secondary-structure indicator, not a bending metric: a straight ideal α-helix reads **110°**, a straight β-strand reads **0°** | Gate 0, 2, 5; every document | Critical |
| F2 | The mover test's null false-positive rate is **41 %** at the admitted minimum n=3 WT crystals (14 % at n=10) vs. the reported "29 % of mutations are movers" | Gate 1 (GO verdict) | Critical |
| F3 | Pseudo-replication: 248 "pairs" are **69 residues**, one of them 61×; Kish n_eff = **13.3**. `SPEC` §2.4 mandates per-variant aggregation; no script implements it | Gates 1–5, every CI and p-value | Critical |
| F4 | Metric noise from realistic coordinate error is **2–7°**, above both the adopted 0.98° floor and the 2.66° median "signal" | Gates 1–5 | Critical |
| F8 | The achievable AUC ceiling is never established. A pure measurement noise of SD **1.85°** with **zero** true effect reproduces the observed 29 % base rate exactly — and then the ceiling for *any* predictor is 0.50 | Gates 2–5 (the whole negative result) | Critical |
| F5 | Validation-family holdout is by **PDB ID only**; no sequence-identity exclusion, contradicting `SPEC` §2.5 and `RESULTS.md` | Gate 4 (the replication) | High |
| F6 | "Core" means `ss in "HE"`; there is **no burial/SASA computation anywhere** in the repo | Gate 5, README's mechanism | High |
| F7 | The CIs on the project's headline positive result are **computed and discarded**; model comparison uses non-overlap of marginal CIs instead of a paired bootstrap | Gate 5 | High |
| F9 | Outcome-dependent QC removes points by the magnitude of the outcome, asymmetrically by substitution class | Gates 3–5 | Medium |
| F10 | Unmodeled residues masquerade as wild-type, admitting real mutants into the WT noise-floor reference | Gates 1–5 | Medium |
| F11 | Principal-axis degeneracy at pseudo-angle 60°; ill-conditioned (s₁/s₀ = 0.69–0.98) at 61–80°, i.e. in the loops Gates 3–4 target | Gate 0 | Medium |
| F12 | The self-test only uses exactly-collinear halves, so it never exercises the PCA fit it exists to justify | Gate 0 | Medium |
| F13 | `is_continuous` has no lower bound: duplicated Cα and 0.3 Å steps pass | Gate 0 | Low |
| F14 | Window-range off-by-one: Gate 1 uses `max-4+1`, Gates 2–5 use `max-4` — different floor dictionaries | Gates 1–5 | Low |
| F15 | The documented cross-check metric is unimplemented and off by 90–180° where it matters | Gate 0 | Low |
| F16 | Assorted doc/code contradictions, dead constants, hardcoded figure values, mislabeled plot | — | Low |

---

## F1 — The metric does not measure bending  *(Critical)* `[A1] [A2]`

`bending_metric.py` is declared "THE SINGLE SOURCE OF TRUTH", and its self-test
advertises `straight = 0.00°`. On ideal, textbook, **perfectly straight-axis**
secondary structure:

| conformation | Cα pseudo-angles | `bending_angle` |
|---|---|---|
| ideal α-helix (r=2.3 Å, rise=1.5 Å, 100°/res) | 90.4° | **110.39°** |
| ideal β-strand (r=0.95 Å, rise=3.3 Å, 180°/res) | 120.1° | **0.00°** |

Robust across the whole plausible parameter range (radius, rise, twist, phase):
α-helix **109.8 ± 3.6°**, β-strand **0.5 ± 0.4°** `[A1]`.

The implementation is faithful to its own formula. The defect is that the formula
does not measure the quantity every downstream document attributes to it. An
α-helix has a straight axis; the Cα trace coils around it at ~100°/residue, so two
consecutive 3-atom chords point ~110° apart. The metric therefore reports **local
coil geometry**, which is essentially a secondary-structure label, not curvature of
the backbone axis.

Consequences, in order of importance:

1. **Gate 2's absolute-bending "skill" is secondary-structure prediction.** On a
   mixed H/E/L population the label's SD is ~47° and removing only the
   three SS-class means cuts it to ~20° — a **~58 % RMSE reduction with zero
   knowledge of bending** `[A2]`. Gate 2 reports 20 % (30.55° vs 37.97°). A 20 %
   reduction is *within what SS propensity alone supplies*, so it is not evidence
   that the model learned anything about bending. `RESULTS.md`'s "~20 % variance
   reduction" should not be read as partial success.

2. **The "insensitivity" diagnosis is over-read.** Gate 2 concludes the model
   "barely responds to a single-residue change", names it "AF2 disease", and
   `RESULTS.md`/`README.md` then generalize it to "backbone bending is governed
   primarily by tertiary interactions". But the label is dominated by a ~110°
   SS-driven term and a point mutation almost never changes SS class, so a *small*
   absolute response is what a correctly-fitted model on this label should produce
   regardless of the underlying physics. Small is not the same as uninformative:
   `[B5]` shows this exact architecture still ranks mutations well (AUC 0.925) when a
   real 3° local effect is present, with the differenced prediction attenuated by
   only 0.89×. So low response magnitude alone does not license the tertiary
   conclusion — the low *Spearman* (0.08) is the relevant evidence, and that is
   confounded by label noise (F8), not by response scale.

3. **Scale mismatch.** The label's dynamic range is ~110°, the coordinate noise is
   2–7° (F4), and the effects being chased are 1–3°. The target sits ~2 % of the
   way into the label's own range.

The delta metric (`bending(mut) − bending(wt)`, same window) is not invalidated by
this — both terms carry the same SS offset and it cancels. F1 is a finding about
Gate 2/Gate 5's *architecture* ("predict absolute, then difference") and about the
framing in all three documents, not about the deltas in Gate 1.

**Fix:** if the intent is axis curvature, use a quantity that is zero on regular
secondary structure — e.g. the angle between the two halves' *helical axes*, or a
7–9 residue smoothed-trace curvature, or write it as a deviation from the
conformation-matched expectation. If the intent is genuinely "local coil geometry",
rename it everywhere and drop `straight = 0` from the self-test advertisement.
Either way, Gate 2's absolute-then-difference architecture needs re-justifying
against an SS-only baseline that the repo currently never computes.

---

## F2 — The mover test's null rate is never computed, and it is large  *(Critical)* `[B1] [B2] [B3]`

`feasibility_t4l.py:171-200` implements:

```
wt_med[s] = median(v₁..vₙ)                       # n redundant WT crystals
wt_sig[s] = 1.4826 · MAD(v₁..vₙ)
mover     ⇔  |b − wt_med[s]|  >  2 · wt_sig[s]
```

Windows are admitted at `len(vals) >= 3`. Two errors compound, both in the same
direction.

**(i) `1.4826·MAD` is badly biased low at small n** `[B1]`:

| n | E[σ̂]/σ | bias |
|---|---|---|
| 3 | 0.671 | **−32.9 %** |
| 5 | 0.823 | −17.8 % |
| 10 | 0.913 | −8.7 % |
| 20 | 0.960 | −4.0 % |

The 1.4826 constant is the *asymptotic* consistency factor; there is no small-sample
correction anywhere in the repo.

**(ii) The reference-error term is omitted.** The docstring (`feasibility_t4l.py:12-14`)
justifies the 2σ cut as clearing a "differencing noise ~ sqrt(2)·sigma" band. That is
the noise of *two single measurements*. The code differences a single mutant against
a **median of n**, for which

$$\mathrm{SD}(b-\widehat{\mathrm{med}}) \;=\; \sigma\sqrt{1+\tfrac{\pi}{2n}}$$

using the median's asymptotic variance πσ²/2n. So a nominal "2σ" cut is only 1.62
true SDs at n=3 and 1.86 at n=10 `[B2]`. The stated derivation does not describe the
implemented test. (This formula is illustrative; the FPR table below is exact Monte
Carlo over the estimator as written, with no asymptotic approximation.)

**Combined null false-positive rate** — mutation has *zero* geometric effect `[B2]`:

| n_WT | null FPR | expected false "movers" of 248 |
|---|---|---|
| **3** | **0.413** | **102** |
| 5 | 0.260 | 65 |
| 10 | 0.145 | 36 |
| 20 | 0.094 | 23 |
| 50 | 0.064 | 16 |

An exact 2σ test with σ known would give 0.0455 (11 of 248).

`RESULTS.md` Gate 1 reports **29 % of mutations above floor (72/248)** and issues
the "GO" verdict on it. At the admitted minimum of n=3 the null alone produces
~41 % — *more movers than were observed*. The true per-window n is not recorded
anywhere, so the null rate cannot be pinned down from the committed artifacts; but
it is above the nominal 4.6 % at **every** n, and no gate computes it. `RESULTS.md`
notes that sparse-window floors "inflate the above-floor fraction" — the size of
that inflation is the whole result, and it is never quantified.

**The two floors do not corroborate each other** `[B3]`. The second estimate uses
`np.std(vals, ddof=1)` at `len(vals) >= 2` (`feasibility_t4l.py:216`), then takes a
median across groups. `ddof=1` is unbiased for the *variance*, not the SD:
E[s]/σ = √(2/π) = 0.798 at n=2, and the median-across-groups step pushes it to
≈0.69–0.71. So the two estimators of the same σ are biased low by ≈−9 % and ≈−30 %
respectively — predicting a ratio near 0.77. The reported ratio is
**0.75/0.98 = 0.77**. `RESULTS.md` presents this agreement as "two independent
noise-floor estimates agree → the detection threshold is empirical, not assumed".
The agreement is what two differently-biased estimators of the same quantity are
*expected* to produce; it is not evidence that either is right.

**Fix:** compute the null directly — permute mutant/WT assignment within each
window and report the mover rate under permutation alongside the observed one. Use
a small-sample-corrected scale estimate (or pooled-across-windows σ with a
hierarchical shrink), and put the reference-median error into the threshold:
`2σ̂·√(1+π/2n)`. Record per-window n in the CSV.

---

## F3 — Pseudo-replication invalidates every n, CI and p-value  *(Critical)* `[C1] [C2]`

From the committed `feasibility_t4l.csv` `[C1]`:

- 248 rows, described throughout as "248 T4L single-mutant **pairs**/**windows**"
- **69 distinct mutated residues**
- residue **99 appears 61 times** (the L99A cavity series), residue 44 15×, residue 26 9×
- **Kish effective n = (Σm)² / Σm² = 13.3**, against a nominal 248 → variance
  inflation ≈ **18.7×**, CI widths ≈ **4.3×** too narrow

Rows sharing a residue share the *same* `wt_med[s]` and the *same* `floor`, so they
are not independent draws even before counting repeat crystals of one variant.

`SPEC_bending_and_pairs.md` §2.4 already specifies the fix and is explicit:

> "Many proteins have several structures per variant: aggregate to **one bending per
> variant** (median across redundant structures); keep the spread as the per-variant
> noise estimate"

**No script implements this.** `feasibility_t4l.py` computes redundant-variant groups
(`var_groups`, line 205) for the independent σ, then discards them and reports one
row per PDB entry. Applying §2.4 to the committed CSV gives a per-residue mover
fraction of **23.2 % (16/69)** versus the reported 29.0 % `[C1]`.

**This corrupts the retrieval metrics harder than the counts.** The differenced
prediction depends only on `(scaffold entry features, wt_seq, mut_seq)` —
`gate2_model_feasibility.py:246-249`. For two crystals of the same variant all three
are *identical*, so `dpred` is **bit-identical**. Every repeat crystal is an exact
tie in the ROC, and the ROC's effective sample size is the number of distinct
variants, not rows.

The bootstraps resample **rows** (`loop_gate.py:260-265`,
`powered_loop_gate.py:198-203`, `gate3_3d.py:148-153`). Under a true null with the
real duplication pattern, a nominal 90 % interval should exclude 0.50 in 10 % of
runs `[C2]`:

| bootstrap | excludes 0.50 under the null |
|---|---|
| row (as implemented) | **62–67 %** |
| cluster (by residue) | 10–11 % |

(Range over three runs at different Monte-Carlo budgets; the cluster bootstrap is
correctly calibrated, the row bootstrap is not.)

On the committed CSV a cluster bootstrap on the mover fraction is **2.1× wider**
than the row bootstrap `[C1]`.

Note the direction: this cuts **against** the project's conclusions as well as for
them. The true intervals are wider than published, so `RESULTS.md`'s
"not distinguishable from chance" is if anything *understated*; but Gate 5's
"the lower bound barely clears 0.50" (from `[0.52, 0.66]`) is not supportable at
all, and Gate 3's loop enrichment (`OR 2.50, Fisher p = 0.041`) is a Fisher exact
test on 248 non-independent rows — with n_eff ≈ 13 that p-value carries no weight.
`RESULTS.md` attributes the 0.524/0.516/0.477 spread entirely to the QC step and
"±0.02–0.03 sampling/seed variation"; pseudo-replication is a larger and unlisted
term.

**Fix:** aggregate to one row per variant per §2.4 before any statistic; keep the
within-variant spread as the per-variant noise. Where rows must stay unaggregated,
cluster-bootstrap by variant and report n_variants next to n_rows.

---

## F4 — The claimed signal is below the metric's own noise  *(Critical)* `[A5]`

Propagating iid Cα coordinate error through `bending_angle` (SD of the metric, in
degrees) `[A5]`:

| conformation | 0.05 Å | 0.10 Å | 0.20 Å | 0.30 Å |
|---|---|---|---|---|
| ideal α-helix | ~1.1 | ~2.3 | ~4.5 | ~6.9 |
| ideal β-strand | ~0.8 | ~1.6 | ~3.0 | ~4.5 |
| collinear | ~0.6 | ~1.2 | ~2.4 | ~3.7 |

A 2.0 Å X-ray structure carries roughly 0.1–0.3 Å positional error on a well-ordered
Cα. Against that:

- adopted per-window floor: **0.98°**
- "independent" corroboration: **0.75°**
- median |Δ| of the 72 movers (`gate2_model_feasibility.ABOVE_FLOOR_SIGNAL`): **2.66°**

So the empirical floor is **2–5× smaller** than coordinate-error propagation
predicts, and the median effect called "signal" is **at or below** the metric's own
noise. Two readings, not mutually exclusive: the MAD estimator is biased low (F2),
and/or the redundant WT crystals are not independent samples of coordinate error
(same crystal form, same refinement protocol, re-refinements of the same data),
so crystal-to-crystal spread underestimates true positional uncertainty. Either way
the floor is too low, in the direction that manufactures movers.

**Fix:** cross-check the empirical floor against a coordinate-error propagation
using each entry's deposited B-factors / DPI, and report both. If they disagree by
2–5×, the crystal-to-crystal floor is measuring reproducibility of a pipeline, not
measurement error.

---

## F5 — The replication's leakage control is by PDB ID, not sequence identity  *(High)*

`SPEC_bending_and_pairs.md` §2.5: "**Hold out whole sequence clusters** for the
benchmark; no family appears in both train and validation." `RESULTS.md`:
"Leakage control: sequence-culled training (≤30 % identity); validation families
fully held out (and explicitly excluded by ID in Gates 4–5)."

What the code does:

- `gate2_model_feasibility.fetch_cull_ids` asks RCSB for `group_by
  sequence_identity, similarity_cutoff: 30, representatives`. That makes the
  **training set internally** non-redundant at 30 %.
- Validation exclusion is `pid.upper() in exclude_ids` — a **PDB-ID set**
  (`powered_loop_gate.py:153`, `gate3_3d.py:53`, `gate2_model_feasibility.py:196`).

A grep for any identity/alignment/homology operation between train and validation
returns nothing. Internal non-redundancy of the training set says nothing about the
identity between a training representative and a *held-out* protein. Nothing stops a
30 %-cluster representative from being 50–70 % identical to a validation protein.

This lands hardest on Gate 4, the project's only independent replication. Its pool is
**human-lysozyme-dominated (27 of 33 movers)**. Human lysozyme (P61626) is a c-type
lysozyme with abundant close homologs among the most-deposited proteins in the PDB;
a c-type lysozyme representative in the culled training set would be far above the
30 % identity the document claims. The same applies to RNase A.

Also: `powered_loop_gate.py:220` prints "excluding {len(val_ids)} val ids **+ T4L**",
but `train_loop_model` excludes only `val_ids`. T4L happens to be absent because
`gate2` filtered it at download time, so the effect is nil — but the guarantee is
printed, not enforced.

**Fix:** cluster train and validation *jointly* at the target identity (mmseqs2 /
`cd-hit`, or RCSB cluster membership) and drop every training chain in any validation
protein's cluster. Report how many chains that removes.

---

## F6 — "Core" is secondary structure, not burial  *(High)*

`gate3_3d.py:195`: `core = np.array([x["ss"] in "HE" for x in rows])`.

`RESULTS.md`: "the **core subset** (packing-dominated, where local-only was 0.477)".
`README.md`: "The effect was most visible in **protein cores, where packing
interactions dominate**."

A grep for `sasa|solvent|accessib|burial|buried|dssp` across the repo returns **no
matches**. There is no burial calculation anywhere. "Core" means "helix or sheet".
A solvent-exposed helix counts as core; a buried loop counts as non-core. The
mechanistic claim — the strongest positive claim in `README.md`, the thing the
project's headline conclusion rests on — is supported by a variable that measures
secondary structure.

The narrative it supports is also circular: the "3D contact environment" feature is
the best available proxy for burial, the "core" split is an SS split, and the label
itself is essentially an SS indicator (F1). That a contact-composition feature helps
more on helix/sheet windows than on loop windows is close to expected under F1
alone, with no packing mechanism required.

**Fix:** compute relative solvent accessibility (Shrake–Rupley on the deposited
structure, or a Cα neighbour-count proxy) and define core by RSA, keeping the SS
split as a separate axis. If the effect survives an RSA split, the mechanism claim
is earned; if it only appears on the SS split, say that.

---

## F7 — The headline positive result is reported with no uncertainty  *(High)*

`gate3_3d.py:196-200`:

```python
if core.sum() >= 20 and 0 < mov[core].sum() < core.sum():
    al, _, _ = auc_ci(mov[core], dl[core])
    a3, _, _ = auc_ci(mov[core], d3[core])
    print(f"  core subset (n=...): local AUC={al:.3f} -> local+3D AUC={a3:.3f}")
```

`auc_ci` returns `(auc, lo, hi)` and computes a 3000-sample bootstrap. The core-subset
call sites **discard both bounds**. So the number `RESULTS.md` builds Gate 5 claim (1)
on — "the core subset … rises from 0.477 to 0.584 … exactly where tertiary context
*should* matter … It supports the diagnosis that the cause is tertiary" — and which
`README.md` elevates to "the most informative result", is the one AUC in the repo
reported **without an interval**, on a subset with a handful of movers, even though
the interval was computed and thrown away.

Separately, the model comparison is not a valid paired test. `gate3_3d.py:205` gates
on `lo3 > hiL` — non-overlap of two *marginal* CIs. The quantity of interest is
AUC₃ − AUC_local on the same rows. Because `auc_ci` reseeds
`np.random.default_rng(0)` on every call, both models are already evaluated on
**identical resamples**, so the paired difference is available for free and is
simply never formed. Non-overlap of marginal intervals is neither necessary nor
sufficient for the paired difference to exclude zero.

`RESULTS.md` also phrases the interval ambiguously: "Overall lift is 0.516 → 0.590,
90 % CI [0.52, 0.66]". `[0.52, 0.66]` is the CI of the local+3D AUC, not of the lift
(0.074). No CI on any lift is computed anywhere.

**Fix:** return and print the core-subset bounds (a one-line change), and add a
paired bootstrap on AUC₃ − AUC_local from shared resample indices.

---

## F8 — The achievable AUC ceiling is never established, and it may be 0.50  *(Critical)* `[B4] [B5]`

This is the finding that ties F2 and F4 to the retrieval results, and it is the one
that decides whether the project's headline conclusion means anything.

Mover labels are not ground truth. They are built from a **noisy observation**
`o = t + e`, where `t` is the true mutation effect and `e` is measurement noise, with
`mover ⇔ |o| > 2σ̂`. Any predictor — including a perfect one that knows `t` exactly —
is scored against those noisy labels. So there is a **ceiling AUC** set by the
noise-to-signal ratio, and it is computable. The repo never computes it.

Constraining to the project's own reported 29 % base rate at its own 1.96° threshold
`[B4]`:

| assumed measurement noise `s` | implied true-effect SD `τ` | **oracle AUC** |
|---|---|---|
| 0.98° (the adopted floor) | 1.56° | **0.841** |
| 1.50° | 1.08° | **0.643** |
| 2.00° | ≈0 | **0.501** |
| 3.00° | ≈0 | **0.502** |
| 4.50° | ≈0 | **0.500** |

And the sharpest form of it: with **no true mutation effect whatsoever** (`τ = 0`),
a measurement noise of SD **1.852°** reproduces the observed 29 % base rate *exactly*
`[B4]`. That is only **1.89×** the project's own adopted floor, and it sits inside
the 2–7° range that coordinate-error propagation gives for this metric (F4).

So the reported AUCs — 0.477, 0.516, 0.524, 0.590, 0.591, 0.700 — admit two readings
that the repo never separates:

- **If the noise really is ~0.98°**, an oracle would score ~0.84, there is real
  headroom, and measuring ~0.52 is a genuine model failure. The project's conclusion
  stands.
- **If the noise is ≥1.85°** — which F4's propagation and F2's estimator bias both
  point to — the mover labels contain essentially no recoverable signal, the ceiling
  *is* ~0.50, and **every reported AUC is already at ceiling.** No model could have
  done better, the negative result is unfalsifiable, and the mechanistic conclusion
  ("the information is tertiary") does not follow from it.

Under the second reading, Gate 3's loop AUC of 0.700 and Gate 5's 0.590 *exceed* the
ceiling, which independently marks them as small-sample fluctuations — the same
conclusion `RESULTS.md` reaches, but for a different and firmer reason.

**A correction to an earlier draft of this audit.** I first suspected the AUC ≈ 0.5
was an estimator artifact: with a label dominated by a ~110° SS term (F1), a
`max_depth=6` tree ensemble might never split on the central one-hot block, making
`dpred ≡ 0` and the AUC exactly 0.5 by construction. **That hypothesis is wrong**, and
the positive control the repo is missing is what refutes it `[B5]`. Rebuilding Gate
2's exact architecture with a genuine 3° local effect planted in the label:

- exact zeros in the differenced prediction: **6.5 %** (not ~100 %)
- distinct `|dpred|` values: **230 of 248** (not degenerate)
- predicted/true delta SD: **0.89×** (little attenuation)
- retrieval of the planted movers: **AUC 0.925**

Predict-absolute-then-difference **passes** this control. The architecture is capable
of recovering exactly the kind of effect being chased, so the observed AUC ≈ 0.5 is
not a model-capacity artifact. That result is in the project's favour and should be
in `RESULTS.md`. It also localizes the problem precisely: the failure is in the
**labels**, not the model — which is what makes the ceiling calculation above the
decisive one.

Still missing and still cheap:

- the fraction of `dpred` values that are exactly zero, and the count of distinct
  `|dpred|` among the 248 (ties are guaranteed by F3 and are credited at 0.5)
- a permutation null: shuffle `movers` against `dpred` to get the achievable AUC
  distribution at this n and tie structure
- the ceiling calculation above, with the noise estimate taken from F4 rather than
  from the biased floor

`README.md` frames the outcome as "The failure of the local model became the result."
That framing requires the experiment to have had the power to succeed. `[B5]` shows
the *model* did; `[B4]` shows the *labels* may not have. Establishing which is the
single highest-value piece of work left in the repo.

---

## F9 — Outcome-dependent QC  *(Medium)*

`loop_gate.py:145`, `powered_loop_gate.py:136`, `gate3_3d.py:132`:

```python
qc_bad = is_conservative(wtaa, mutaa) and abs(delta) > 10
```

Points are removed **on the magnitude of the outcome**. There is no independent
evidence of a crystallographic artifact — the criterion *is* the outcome. It is also
asymmetric: a non-conservative substitution with |Δ| > 10° is kept, a conservative
one is dropped, and conservative-vs-not is a function of the substitution identity,
which is exactly what the model reads as input. This is selection on the dependent
variable, correlated with a predictor.

`RESULTS.md` presents it as "Crystallographic-artifact QC: conservative substitutions
with implausibly large bends (>10°) flagged and removed (e.g. V111M −32.6°)", and
uses it to explain the 0.524 → 0.516 shift. The shift is caused by the filter, and
the filter's direction is set by the outcome.

Note that "implausibly large" is itself doing work that F4 undermines: at 0.2–0.3 Å
coordinate error a helical window's metric has a 4–7° SD, so |Δ| ≈ 10–30° is roughly
2–5σ of pure noise — surprising but not impossible, and not obviously an artifact.

**Fix:** decide exclusion from artifact evidence independent of Δ (occupancy, altloc,
B-factors, crystal-contact flags, resolution gap — several already specified in
`SPEC` §2.3 and none implemented), and report results with and without the filter.

---

## F10 — Unmodeled residues masquerade as wild-type  *(Medium)*

The mutation set of a structure is built only over residues **present in that
structure** (`feasibility_t4l.py:159-163`, and repeated in `gate2`, `loop_gate`,
`powered_loop_gate`, `mover_composition`, `gate3_3d`):

```python
muts = tuple(sorted((rs, consensus[rs], aa)
                    for rs, (aa, _) in res.items()
                    if rs in consensus and aa != consensus[rs]))
```

If a mutation site is disordered, truncated, or otherwise unmodeled, it cannot appear
in `muts`, so the structure is classified `len(m) == 0` → **wild-type** and enters
`wt_ids`, the set that defines `wt_med` and `wt_sig` for every window. A real mutant
contaminates the reference *and* the noise floor. By the same mechanism a double
mutant with one site unmodeled is classified as a single mutant and enters the
benchmark with a mis-attributed cause.

Both failures bias the floor upward or downward unpredictably and violate the
"matched pair" rules in `SPEC` §2.2.

Related: the consensus admission threshold is inconsistent between gates —
`>= 10` observations for T4L (`feasibility_t4l.py:155`) but `>= 5` in
`powered_loop_gate.py:94`, so "wild-type" is a different definition in Gate 1 and
Gate 4.

**Fix:** require full coverage of all consensus positions (or at minimum of the
window and every position where any structure differs) before admitting a structure
as WT; record modeled-coverage per structure and exclude on it explicitly.

---

## F11 — Principal-axis degeneracy  *(Medium)* `[A3]`

For a symmetric Cα triple with apex angle φ and bond length d, the variance along
the chord is (2/3)d²sin²(φ/2) and perpendicular to it (2/9)d²cos²(φ/2). These are
equal when tan(φ/2) = 1/√3, i.e. **φ = 60° exactly**. Verified `[A3]`:

| φ | s₁/s₀ | angle(axis, chord) |
|---|---|---|
| 120° | 0.333 | 0° |
| 91° (α-helix) | 0.567 | 0° |
| 80° | 0.688 | 0° |
| 70° | 0.825 | 0° |
| 61° | 0.980 | 0° |
| **60°** | **1.000** | **82.98°** (tie-break) |
| 59° | 0.980 | **90°** |

Two failures. Below 60° the fitted axis is **perpendicular** to the chain, and
`bending_metric.py:51` orients it by `sign(dot(axis, P[-1]-P[0]))` — whose argument
is exactly 0 there, so a ~180° flip is decided by rounding. That regime is rarely
reached: physical Cα pseudo-angles bottom out near 75–80°.

The second failure is reached constantly. s₁/s₀ = 0.69–0.98 for φ = 61–80° means the
principal axis is **ill-conditioned** precisely in tight turns — i.e. in the loop
population that Gates 3 and 4 are built on. There the fitted direction is maximally
noise-amplifying, which is consistent with the elevated metric noise in F4 and means
the loop gates run the metric in its least stable regime. `RESULTS.md` Gate 0 records
the limitation as "Cα-only, one geometric definition"; the conditioning of that
definition is not mentioned.

**Fix:** guard on `s[1]/s[0]` (not just `s[0] < _EPS`) and either abstain or fall
back to the chord when the axis is near-degenerate; propagate that flag into
`coverage_class` as `SPEC` §2.6 already envisages.

---

## F12 — The self-test cannot detect F1, F11, or F13  *(Medium)* `[A4]`

`bending_metric.py:99-103` builds every test window as two straight arms meeting at
an apex, so **both halves are exactly collinear** (s₁/s₀ = 0 to machine precision
for all six test angles) `[A4]`. The SVD returns the chord trivially. The test
validates the `arccos`, not the principal-axis fit it exists to justify — the
docstring's claim that "all three points contribute, so it is more robust than an
endpoint vector" is untested, and replacing the PCA with `axis := P[2]-P[0]` passes
the entire self-test unchanged `[A4]`.

`RESULTS.md` Gate 0 states: "**Claim:** the bending metric is correct and free of
sign/orientation bugs. **Numbers:** synthetic windows with known bends recovered
exactly." What was verified is narrower than what is claimed. No test covers a
non-collinear half, an ideal α-helix or β-strand, the φ→60° degeneracy, a
near-degenerate axis, or coordinate-noise sensitivity — the four things that turn out
to matter.

**Fix:** add assertions on ideal secondary structure (pinning whatever the intended
values are), on `s[1]/s[0]` conditioning, and on noise sensitivity.

---

## F13 — `is_continuous` has no lower bound  *(Low)* `[D1]`

`bending_metric.py:87-93` checks only `d <= 4.5`. A window with a **duplicated Cα**
(0.0 Å step) or 0.3 Å steps passes and returns a confident `0.00°` `[D1]`. Such
windows arise from duplicated atoms, altloc mishandling, and modelling errors — and
`parse_ca`'s altloc logic (highest occupancy, `best[key]` keyed on `(chain, resSeq)`)
keeps one atom per residue, so this is a guard against upstream data problems rather
than a live bug in this pipeline. A lower bound near 2.5 Å would keep genuine
cis-peptides (~2.9 Å) and reject the rest.

`_CA_CA_IDEAL = 3.8` (line 83) is defined and never referenced anywhere in the repo.

---

## F14 — Window-range off-by-one between the gates  *(Low)* `[D2]`

For residues `lo..hi` the last window that fits is `hi-4`.

| site | expression | last start | |
|---|---|---|---|
| `feasibility_t4l.py:172` | `range(min, max-4+1)` | `hi-4` | correct |
| `gate2_model_feasibility.py:162` | `range(min, max-4)` | `hi-5` | drops one |
| `loop_gate.py:121` | `range(min, max-4)` | `hi-5` | drops one |
| `mover_composition.py:64` | `range(min, max-4)` | `hi-5` | drops one |
| `gate3_3d.py:111` | `range(min, max-4)` | `hi-5` | drops one |
| `powered_loop_gate.py:111` | `range(min, max-4)` | `hi-5` | drops one |

Gate 1 builds `wt_med`/`wt_sig` over one more window than Gates 2–5, and `pooled` is
a median over that dictionary — so the "same" floor is not numerically the same
between gates. `featurize_chain` and `build_training` use `range(lo+3, hi-4)` and
drop the last window too. Small in magnitude, but it is one of the unlisted
contributors to the 0.524/0.516 drift that `RESULTS.md` attributes wholly to QC.

---

## F15 — The documented cross-check metric is wrong and unimplemented  *(Low)* `[D3]`

`bending_metric.py:31-33` describes an optional cross-check — "unsigned sum of the
three C-alpha pseudo-bond exterior angles. Should track this within a few degrees;
large disagreement flags an S-shaped (sign-cancelling) window worth inspecting."

Measured `[D3]`:

| conformation | `bending_angle` | Σ exterior angles | disagreement |
|---|---|---|---|
| ideal α-helix | 110.39° | 268.90° | **158.5°** |
| ideal β-strand | 0.00° | 179.59° | **179.6°** |
| collinear | 0.00° | 0.00° | 0° |
| bent 90° | 90.00° | 90.00° | 0° |

The two agree only on the collinear windows of the self-test; on real secondary
structure they differ by 90–180°. The cross-check is also **never implemented**, so
the stated safeguard against S-shaped windows does not exist — and an S-shaped
window remains a real failure mode of a two-half angle (equal and opposite bends
read as straight).

---

## F16 — Documentation and reproducibility defects  *(Low)*

- **Label region misstated.** `gate2_model_feasibility.entry_features` docstring says
  entry geometry is "strictly upstream of the curvature that defines the label, which
  lives in `i+1..i+4`". The label is `bending_angle` over `i..i+4` and depends on
  point `i`, which is also `entry_features`' last point. The two share one atom but no
  angle, so there is no numeric leakage — but the stated guarantee is not the
  implemented one. `bending_metric.py`'s docstring describes the entry region as
  `i-2, i-1`, while `gate2` uses `i-3..i`; the two "canonical" descriptions disagree.
- **`summary_figure.py` hardcodes results.** Its header calls the values "the
  established outputs of gate2 / loop_gate / powered_loop_gate", but they are literals
  that cannot track the scripts. It labels `0.524` as `"Core\n(global)"` — `0.524` is
  the *global, pre-QC* AUC; `RESULTS.md`'s core-only figure is `0.477`. Two different
  populations under one label.
- **`summary_auc.png` is not committed**, though `.gitignore` says "Keep the figures
  (*.png)". Every other figure is present.
- **No pinned dependencies and no test runner.** `bending_metric.py`'s self-test runs
  only under `__main__`; nothing else in the repo has tests. There is no
  `requirements.txt`, so the numpy/scipy/sklearn versions that produced the committed
  figures are unrecorded.
- **`parse_ss` catches only `ValueError`.** A truncated HELIX/SHEET line raises
  `IndexError` at `line[19]` and aborts the run. (The column offsets themselves are
  correct against the PDB spec — HELIX chain at col 20, `initSeqNum` 22–25,
  `endSeqNum` 34–37; SHEET chain at col 22, `initSeqNum` 23–26 — all verified.)
- **Arbitrary chain tie-break.** `pick_t4l_chain` / `pick_chain` break score ties by
  dict iteration order. For a crystal with two NCS copies of the same chain the choice
  is arbitrary, and the two copies have genuinely different local geometry, so an
  arbitrary few tenths of a degree enters each measurement.
- **`env_at` count/composition mismatch.** `gate3_3d.py:45-46` builds `cnt` from
  `valid = mask & (aai >= 0)` but reports `n_contacts = mask.sum()`. Currently dead
  (`parse_ca` admits only the 20 standard residues, so `aai >= 0` always), but the two
  numbers will silently disagree the moment non-standard residues are admitted.
- **`gate3_3d` 3D feature is train/serve mismatched.** In training, `env_at(g, r)` is
  computed from the *same chain* whose coordinates produce the label; at validation it
  comes from a fixed WT scaffold while the label comes from the mutant structure. The
  contact set within 8 Å is itself a function of the local conformation, so the
  training feature is partly a readout of its own target in a way the validation
  feature cannot be. This inflates the absolute-bending RMSE gain (30.5° → 27.9°) that
  `RESULTS.md` cites as evidence the 3D feature "captures a thin slice of 3D", while
  leaving the delta retrieval unaffected — consistent with what is observed.

---

## What would make the conclusions supportable

In dependency order. The first item is the one that decides whether anything else
matters.

1. **Establish the measurement noise independently, then compute the ceiling**
   (F8, F4, F2). Propagate deposited B-factors / DPI through the metric to get a
   coordinate-error-based noise estimate, and compute the oracle AUC at that noise
   and the observed base rate. If the ceiling is ~0.5, the retrieval gates are
   measuring nothing and no amount of modelling will change that. If it is ~0.84,
   the negative result is real and the rest of the list is what makes it publishable.
   This is a day of work and it determines the status of the entire project.
2. **Decide what the metric should measure** (F1). If axis curvature is the target,
   the current definition needs replacing, and Gate 2's
   predict-absolute-then-difference architecture needs a fresh justification against
   an SS-only baseline that the repo never computes.
3. **Aggregate to one row per variant** (F3, `SPEC` §2.4 — already specified).
   Re-derive every n, CI and p-value; cluster-bootstrap what remains.
4. **Recalibrate the mover test** (F2): permutation null, small-sample-corrected
   scale estimate, reference-median error in the threshold, per-window n in the CSV.
5. **Enforce the leakage guarantee** (F5) by joint clustering, and **measure burial**
   (F6) before making the packing argument.
6. **Report the intervals already being computed** (F7), make the model comparison
   paired, and add the tie-fraction and permutation-null diagnostics (F8).

Items 5, 6 and most of F16 are hours of work. Items 1–4 change the headline numbers.

A closing note on framing. `README.md` presents the negative result as robust and
well-characterized: "The signal exists. The predictor does not." The evidence
supports a different and more uncomfortable statement: **the first clause is the one
in doubt.** Gate 1 — "the signal exists" — is the weakest gate in the chain, not the
strongest, because its 29 % above-floor rate is reproduced by pure noise at a
plausible coordinate error (F2, F4, F8), and everything downstream inherits its
labels. The second clause, "the predictor does not", is the one part of the project
that survives a positive control `[B5]`: the architecture demonstrably recovers a
planted local effect, so its failure on real labels is not a capacity artifact.

That inversion is worth stating plainly, because it changes what the project is. As
written it reads as a well-powered negative result about protein physics. What the
artifacts support is a negative result about **this measurement**: a 5-Cα PCA angle
on 2.0 Å crystal structures does not resolve single-mutation backbone effects,
because its noise is comparable to the effects. That is a real, useful, publishable
finding about metric design — and it is a different paper from the one the README
describes.

The documents' habit of reporting limitations candidly (`RESULTS.md`'s insistence on
CIs touching chance, its refusal to call 0.59 a signal, its "the noise and the
signal are the same size is itself part of the finding") is genuine and unusual, and
it is the reason this audit could be done at all. The gap is not honesty; it is that
the limitations which dominate the result are the ones not on the list.
