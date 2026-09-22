"""
v2 MEASUREMENT GATE -- can mutation-induced local geometry change be measured at
all from this class of crystal data?

No machine learning. No PDB downloads. This gate decides whether the *label* is
worth building a benchmark on, before any model is trained again. It exists
because v1's conclusion ("local sequence cannot predict bending") is not
separable from "the measurement could not resolve the effect" -- see AUDIT.md
F1/F2/F4/F8.

Stages, then a verdict:

  1  CRITERION 1  a straight axis must read ~0     (all regular SS types, all phases)
  2  CRITERION 2  a controlled bend must be tracked (slope, linearity, monotonicity)
  3  CRITERION 3  stability under coordinate noise  (noise gain, deg per Angstrom)
  4  conditioning / degeneracy, and twist coverage across the whole range
  5  per-structure uncertainty: B-factor -> sigma, and a calibrated mover threshold
  6  ORACLE AUC CEILING and the minimum detectable effect size
  7  sensitivity of the verdict to coordinate error and crystals per variant

Decision rule (fixed in advance, per the v2 plan)
------------------------------------------------
Let tau* = the true axis-bend SD (degrees) a candidate needs for a PERFECT
predictor to reach a given oracle AUC, at realistic coordinate error. Stage 6
reports tau* and the Ca displacement it corresponds to.

  * no candidate reaches oracle AUC ~0.75 at a physically plausible tau
    -> MEASUREMENT BOTTLENECK. The project's answer is about measurement, not
       about sequence. Do not resume the v1 ML gates.
  * some candidate reaches ~0.75-0.85 at a plausible tau
    -> the label is usable. Promote it to metric v2 and build benchmark v2.

"Physically plausible" is anchored by reporting tau* as a Ca displacement: point
mutations typically move local backbone by ~0.1-0.5 A.

    python3 measurement_gate.py            # all stages
    python3 measurement_gate.py 1 2 6      # selected stages
"""
from __future__ import annotations

import sys

import numpy as np

from metrics import (CAL_BISECTOR, IDEAL, BisectorAxis, HybridAxisBend,
                     SmoothedChord, default_candidates, dirichlet_shrink,
                     helix_on_arc, sigma_from_bfactor)

SIGMA_XYZ = (0.05, 0.10, 0.20, 0.30, 0.50)
SIGMA_REF = 0.20            # the value the verdict is quoted at
PHASES = (0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0)
BENDS = (0.0, 2.5, 5.0, 10.0, 20.0, 30.0, 45.0, 60.0)
SS_LIST = tuple(IDEAL)
SS_MAIN = ("alpha_helix", "beta_strand", "three_ten", "pp_ii")
STRAIGHT_TOL = 5.0
TARGET_FPR = 0.05
N_WT, N_MUT = 5, 2          # a realistic v2 aggregation


def hdr(n, title):
    print()
    print("=" * 78)
    print(f"  STAGE {n}. {title}")
    print("=" * 78)


def _ca(cand, ss, bend, phase=0.0):
    r, rise, turn = IDEAL[ss]
    return helix_on_arc(cand.span, r, rise, turn, bend, phase)


def _mean(vals):
    v = np.array([x for x in vals], float)
    return float("nan") if np.all(np.isnan(v)) else float(np.nanmean(v))


def phase_slopes(cand, ss, bends=BENDS, phases=PHASES):
    """Per-phase SIGNED slope of the metric vs true arc bend.

    Fitting the phase-AVERAGED signed value -- which an earlier version of this
    gate did -- silently cancels any response whose sign depends on helical
    phase, and reports a slope near zero for a metric that is in fact responding
    strongly. That produced a spurious "ceiling 0.50" for v1 on alpha-helices.
    See CORRECTIONS.md.

    Returns (per_phase_slopes, mean_signed, mean_abs, sign_consistent) where
      mean_signed  : what the broken version measured; near 0 under cancellation
      mean_abs     : mean |per-phase slope| -- the sensitivity a MAGNITUDE score
                     such as |delta| actually gets, and the right input to a
                     mover-retrieval ceiling
      sign_consistent : False if the response changes sign across phase, which
                     scrambles any SIGNED analysis (correlation, direction)
    """
    r, rise, turn = IDEAL[ss]
    A = np.column_stack([bends, np.ones(len(bends))])
    out = []
    for ph in phases:
        vals = np.array([cand(helix_on_arc(cand.span, r, rise, turn, b, ph))
                         for b in bends], float)
        if np.isnan(vals).any():
            return None, float("nan"), float("nan"), None
        out.append(float(np.linalg.lstsq(A, vals, rcond=None)[0][0]))
    s = np.array(out)
    return (s, float(s.mean()), float(np.abs(s).mean()),
            bool(not (s.min() < 0 < s.max())))


# --------------------------------------------------------------------------- #
def stage1():
    hdr(1, "CRITERION 1 -- does a STRAIGHT axis read zero?")
    print("  Perfectly regular secondary structure on a perfectly straight axis: an")
    print("  axis-bend metric MUST read ~0. 'nan' = the candidate declares itself")
    print("  inadmissible here, which is a pass-by-abstention, not a silent error.\n")
    print(f"  {'candidate':>20} {'span':>5} " +
          " ".join(f"{s[:9]:>10}" for s in SS_LIST))
    worst = {}
    for c in default_candidates():
        cells, w = [], 0.0
        for ss in SS_LIST:
            v = np.array([c(_ca(c, ss, 0.0, ph)) for ph in PHASES], float)
            if np.all(np.isnan(v)):
                cells.append(f"{'abstain':>10}")
                continue
            m = float(np.nanmax(np.abs(v)))
            cells.append(f"{m:>10.2f}")
            w = max(w, m)
        worst[c.name] = w
        print(f"  {c.name:>20} {c.span:>5} " + " ".join(cells))
    print("\n  (worst |metric| in deg over the phase sweep; straight axis => target 0)")
    print(f"\n  {'candidate':>20} {'worst':>8}  verdict (tol {STRAIGHT_TOL:.0f} deg)")
    for n, w in worst.items():
        print(f"  {n:>20} {w:>8.2f}  {'PASS' if w < STRAIGHT_TOL else 'FAIL'}")
    print("\n  v1_pca5 does not merely confuse helix with strand: it reads a DIFFERENT")
    print("  large offset for every regular conformation, i.e. it is a secondary-")
    print("  structure classifier with a curvature-shaped name (AUDIT.md F1).")
    print("\n  Why sliding-window smoothing alone cannot fix it -- the Dirichlet")
    print("  shrink factor depends on the twist, so no single width collapses every")
    print("  conformation (residual axis wobble, Angstrom):")
    print(f"\n  {'w':>3} " + " ".join(f"{s[:9]:>10}" for s in SS_LIST))
    for w in (3, 4, 5, 6, 7):
        print(f"  {w:>3} " + " ".join(
            f"{IDEAL[s][0]*dirichlet_shrink(w, IDEAL[s][2]):>10.3f}" for s in SS_LIST))
    print("\n  The bisector construction sidesteps this: consecutive Ca bisectors are")
    print("  radial vectors, so their cross product is along the axis for ANY twist.")
    print("  It is exact (0.000) wherever it is admissible, and abstains at twist")
    print("  ~180 deg where bisectors go antiparallel -- which is precisely where")
    print("  w=4 smoothing is itself exact. The hybrid dispatches on that condition.")
    return worst


def stage2():
    hdr(2, "CRITERION 2 -- is a CONTROLLED bend tracked in the right direction?")
    print("  Ground truth: a regular helix wound around a circular arc whose axis")
    print("  turns by a known angle across the span. Reported slope is raw response")
    print("  per degree of TRUE full-span bend. For the calibrated hybrid the target")
    print("  slope is 1.0 (it reports true bend directly); for raw constructions the")
    print("  slope is a geometric constant and only linearity and sign matter.\n")
    out = {}
    for c in default_candidates():
        print(f"  {c.name}")
        for ss in SS_MAIN:
            vals = np.array([_mean([c(_ca(c, ss, b, ph)) for ph in PHASES])
                             for b in BENDS])
            if np.all(np.isnan(vals)):
                print(f"      {ss:12s} abstains (inadmissible construction here)")
                out[(c.name, ss)] = (float("nan"), float("nan"), False)
                continue
            A = np.column_stack([BENDS, np.ones(len(BENDS))])
            slope, icpt = np.linalg.lstsq(A, vals, rcond=None)[0]
            pred = A @ [slope, icpt]
            sst = ((vals - vals.mean()) ** 2).sum()
            r2 = 1 - ((vals - pred) ** 2).sum() / sst if sst > 1e-12 else float("nan")
            _, m_signed, m_abs, consistent = phase_slopes(c, ss)
            print(f"      {ss:12s} " + " ".join(f"{v:6.1f}" for v in vals))
            print(f"      {'':12s} phase-avg slope {m_signed:+.4f}   "
                  f"mean|per-phase| {m_abs:+.4f}   R2 {r2:7.4f}")
            print(f"      {'':12s} sign consistent across phase: "
                  f"{'YES' if consistent else 'NO  <-- signed analyses scrambled'}")
            out[(c.name, ss)] = (float(m_abs), float(r2), bool(consistent))
        print()
    print("  Two slopes are reported because they answer different questions.")
    print("  'phase-avg slope' fits the phase-averaged signed value; 'mean|per-phase|'")
    print("  averages the magnitude of each phase's own slope. They diverge exactly")
    print("  when the response changes SIGN with helical phase -- the same physical")
    print("  bend reading positive at one point in the turn and negative a third of a")
    print("  turn later. When that happens:")
    print("    * a MAGNITUDE score (|delta|) still works, at mean|per-phase| gain;")
    print("    * any SIGNED analysis -- correlation, direction of change -- is")
    print("      scrambled by a nuisance parameter (where in the turn the site sits).")
    print("  v1's Gate 2 reported Spearman(pred, obs) = 0.077 and read it as model")
    print("  insensitivity. A phase-dependent sign would produce that number on its")
    print("  own. The ceiling in stage 6 uses mean|per-phase|, the magnitude gain.")
    return out


def stage3():
    hdr(3, "CRITERION 3 -- stability under coordinate noise")
    print("  SD of each metric (deg) under iid Gaussian Ca perturbation, at a mild")
    print("  10 deg bend -- a realistic operating point. This is the measurement")
    print("  noise every downstream statistic inherits.\n")
    rng = np.random.default_rng(0)
    N = 2500
    out = {}
    for c in default_candidates():
        print(f"  {c.name}")
        print(f"      {'SS':>12} " + " ".join(f"{s:>8.2f}A" for s in SIGMA_XYZ))
        for ss in SS_MAIN:
            P0 = _ca(c, ss, 10.0)
            row = []
            for sg in SIGMA_XYZ:
                v = np.array([c(P0 + rng.normal(0, sg, P0.shape)) for _ in range(N)])
                v = v[np.isfinite(v)]
                sd = float(v.std()) if len(v) > 50 else float("nan")
                row.append(sd)
                if abs(sg - SIGMA_REF) < 1e-9:
                    out[(c.name, ss)] = sd
            cells = " ".join("     abst" if np.isnan(x) else f"{x:>9.2f}" for x in row)
            print(f"      {ss:>12} {cells}")
        print()
    print("  Noise gain is near-linear in sigma_xyz, so these scale. v1 adopted an")
    print(f"  empirical floor of 0.98 deg; at sigma_xyz={SIGMA_REF} A the v1 metric's")
    print("  own noise is several times that (AUDIT.md F4). Note the smoothed and")
    print("  bisector constructions are NOT automatically quieter -- averaging helps,")
    print("  but the bisector differentiates twice, which costs noise. That trade is")
    print("  exactly what stage 6 prices.")
    return out


def stage4():
    hdr(4, "Degeneracy and twist coverage across the whole range")
    print("  (a) s1/s0 of any PCA step. ->1 means the axis direction is undetermined")
    print("      and a ~180 deg flip is decided by rounding (AUDIT.md F11).\n")
    print(f"  {'candidate':>20} " + " ".join(f"{s[:9]:>10}" for s in SS_LIST))
    for c in default_candidates():
        cells = []
        for ss in SS_LIST:
            vs = [c.conditioning(_ca(c, ss, b, ph))
                  for b in (0.0, 30.0) for ph in PHASES]
            cells.append(f"{max(vs):>10.3f}")
        print(f"  {c.name:>20} " + " ".join(cells))

    print("\n  (b) Twist sweep -- expose the handover band where neither branch")
    print("      is reliable. Straight axis, so target is 0.\n")
    bis, sm, hyb = BisectorAxis(9), SmoothedChord(9, 4), HybridAxisBend(9)
    print(f"  {'twist':>6} {'min|sin|':>9} {'branch':>10} {'bisector':>10} "
          f"{'smooth_w4':>10} {'HYBRID':>10}")
    gaps = []
    for turn in range(60, 186, 10):
        rise = float(np.clip(1.5 + (turn - 100) / 80.0 * 1.8, 1.1, 3.4))
        r = np.sqrt(max(3.8 ** 2 - rise ** 2, 1e-6)) / (2 * np.sin(np.radians(turn) / 2))
        P = helix_on_arc(9, r, rise, turn, 0.0)
        bv, hv = bis(P), hyb(P)
        print(f"  {turn:>6} {bis.degeneracy(P):>9.4f} {hyb.branch(P):>10} "
              f"{('       nan' if bv != bv else f'{bv:10.3f}')} {sm(P):>10.3f} "
              f"{hv:>10.3f}")
        if not np.isfinite(hv) or abs(hv) > STRAIGHT_TOL:
            gaps.append((turn, round(hv, 2)))
    dense_gaps = []
    for turn in np.arange(165.0, 195.01, 0.05):
        rise = float(np.clip(1.5 + (turn - 100) / 80.0 * 1.8, 1.1, 3.4))
        r = np.sqrt(max(3.8 ** 2 - rise ** 2, 1e-6)) / (
            2 * np.sin(np.radians(turn) / 2))
        P = helix_on_arc(9, r, rise, turn, 0.0)
        if hyb.branch(P) == "abstain":
            dense_gaps.append(turn)
    print(f"\n  coarse-grid inadmissible twists: {gaps}")
    if dense_gaps:
        intervals = []
        start = previous = dense_gaps[0]
        for turn in dense_gaps[1:]:
            if turn - previous > 0.051:
                intervals.append((start, previous))
                start = turn
            previous = turn
        intervals.append((start, previous))
        print("  dense-grid abstention bands: " + ", ".join(
            f"{lo:.2f} to {hi:.2f} deg" for lo, hi in intervals))
    print("  A direct switch at ~171.4 deg created a 3.45 deg jump on a perfectly")
    print("  straight axis. The hybrid now abstains in the handover band; report")
    print("  its empirical coverage cost on real structures.")

    print("\n  (c) DIFFERENTIAL BIAS -- the part of criterion-1 bias that does NOT")
    print("      cancel in a WT->mutant difference.\n")
    print("      v2 only ever uses the metric as a difference between matched windows")
    print("      of the same protein at the same position, so a CONSTANT offset")
    print("      cancels. What survives is the offset's sensitivity to conformation:")
    print("      if a mutation nudges the local twist, a twist-coupled offset leaks")
    print("      straight into the delta and is indistinguishable from real bending.")
    print("      d(metric)/d(twist) at a straight axis, in deg of spurious bend per")
    print("      deg of twist change:\n")
    print(f"      {'candidate':>20} " + " ".join(f"{s[:9]:>10}" for s in SS_MAIN))
    for c in default_candidates():
        cells = []
        for ss in SS_MAIN:
            r, rise, turn = IDEAL[ss]
            h = 2.0
            vs = []
            for t in (turn - h, turn + h):
                rr = np.sqrt(max(3.8 ** 2 - rise ** 2, 1e-6)) / (
                    2 * np.sin(np.radians(t) / 2))
                vs.append(_mean([c(helix_on_arc(c.span, rr, rise, t, 0.0, ph))
                                 for ph in PHASES]))
            cells.append("     abst" if np.isnan(vs).any()
                         else f"{abs(vs[1] - vs[0]) / (2 * h):>10.3f}")
        print(f"      {c.name:>20} " + " ".join(cells))
    print("\n      Read against the stage-3 noise. A metric with 0.000 here is")
    print("      immune to conformational drift between the two crystals; a metric")
    print("      with ~0.5 turns a 2 deg twist change into ~1 deg of fake bending.")
    print("      This is the real decision axis: the bisector family trades noise")
    print("      for zero differential bias, the smoothed family the reverse.")
    print("      More crystals per variant shrink the IID part of the noise but not")
    print("      a twist-coupled systematic -- nor the correlated part of real")
    print("      crystal-to-crystal variation (packing, refinement era, construct).")
    print("      See CORRECTIONS.md C2: sigma^2(n) = sigma^2_iid/n + sigma^2_syst.")


def stage5():
    hdr(5, "Per-structure uncertainty, and a CALIBRATED mover threshold")
    print("  (a) B-factor -> per-axis positional spread, sigma_1d = sqrt(B/(8 pi^2)):\n")
    print(f"      {'B (A^2)':>9} {'sigma_1d (A)':>14}")
    for B in (5, 10, 15, 20, 30, 40, 60):
        print(f"      {B:>9} {sigma_from_bfactor(B):>14.3f}")
    print("\n      This is positional SPREAD (thermal + static disorder) -- an upper")
    print("      bound on refinement coordinate error, and the right quantity for")
    print("      'how far apart will this Ca sit in two independent crystals'. The")
    print("      Cruickshank DPI is the lower-bound counterpart but needs reflection")
    print("      counts, absent from a coordinate file. Ca B of 10-30 A^2 gives")
    print(f"      0.11-0.20 A, consistent with the {SIGMA_REF} A used for the verdict.")

    print("\n  (b) v1 fixed the mover threshold at 2*sigma_hat and never checked its")
    print("      null rate; it was 14-41%, not ~5% (AUDIT.md F2). v2 must SOLVE for")
    print("      the multiplier that delivers a target FPR under the aggregation")
    print("      actually used. Simulated, null = mutation has no effect:\n")
    rng = np.random.default_rng(1)
    T = 200_000
    print(f"      {'n_wt':>5} {'n_mut':>6} {'SD(delta)':>11} "
          f"{'k for FPR ' + format(TARGET_FPR, '.0%'):>18} {'FPR if k=2':>12}")
    for n_wt, n_mut in ((1, 1), (3, 1), (5, 1), (5, 2), (10, 3), (20, 5)):
        wt = np.median(rng.standard_normal((T, n_wt)), axis=1)
        mut = np.median(rng.standard_normal((T, n_mut)), axis=1)
        d = np.abs(mut - wt)
        sd = float(np.sqrt(np.mean((mut - wt) ** 2)))
        k = float(np.quantile(d, 1 - TARGET_FPR))
        print(f"      {n_wt:>5} {n_mut:>6} {sd:>11.3f} {k:>18.2f} "
              f"{np.mean(d > 2.0):>12.1%}")
    print("\n      Units are the single-measurement sigma. The correct multiplier is")
    print("      not 2 and depends on how many crystals each side has -- so it must")
    print("      be computed per variant pair, not fixed globally.")


def stage6(noise=None, slopes=None):
    hdr(6, "ORACLE AUC CEILING and minimum detectable effect")
    print("  A perfect predictor of the TRUE effect is still scored against labels")
    print("  built from a NOISY measurement, so there is a ceiling. v1 never computed")
    print("  it (AUDIT.md F8). Here it comes from each candidate's OWN measured slope")
    print(f"  and noise at sigma_xyz = {SIGMA_REF} A.\n")

    cands = default_candidates()
    if noise is None or slopes is None:
        print("  (measuring slope and noise; run stages 2 and 3 for detail)\n")
        slopes, noise = {}, {}
        rng = np.random.default_rng(0)
        for c in cands:
            for ss in SS_MAIN:
                vals = np.array([_mean([c(_ca(c, ss, b, ph)) for ph in PHASES])
                                 for b in BENDS])
                if np.all(np.isnan(vals)):
                    slopes[(c.name, ss)] = float("nan")
                    noise[(c.name, ss)] = float("nan")
                    continue
                slopes[(c.name, ss)] = phase_slopes(c, ss)[2]   # mean|per-phase|
                P0 = _ca(c, ss, 10.0)
                v = np.array([c(P0 + rng.normal(0, SIGMA_REF, P0.shape))
                              for _ in range(2500)])
                v = v[np.isfinite(v)]
                noise[(c.name, ss)] = float(v.std()) if len(v) > 50 else float("nan")

    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(7)
    T = 60_000

    def make_ceiling(slope, sig_m):
        """Build a ceiling(tau) closure, drawing the noise ONCE for speed.

        The aggregation noise and the FPR-calibrated threshold depend only on
        (sig_m, N_WT, N_MUT), so they are computed once per candidate rather
        than once per tau -- the bisection below calls this ~50 times per row.
        """
        if not np.isfinite(slope) or not np.isfinite(sig_m) or slope <= 1e-9:
            return lambda tau: float("nan")
        agg = (np.median(rng.normal(0, sig_m, (T, N_MUT)), axis=1)
               - np.median(rng.normal(0, sig_m, (T, N_WT)), axis=1))
        null = (np.median(rng.normal(0, sig_m, (T, N_MUT)), axis=1)
                - np.median(rng.normal(0, sig_m, (T, N_WT)), axis=1))
        thr = float(np.quantile(np.abs(null), 1 - TARGET_FPR))
        z = rng.standard_normal(T)              # reused, scaled by tau

        def ceiling(tau):
            t = tau * z
            mov = np.abs(slope * t + agg) > thr
            if not (0 < mov.sum() < T):
                return float("nan")
            return float(roc_auc_score(mov, np.abs(t)))
        return ceiling

    TAUS = (1.0, 2.0, 5.0, 10.0, 20.0)
    print(f"  aggregation: {N_WT} WT crystals, {N_MUT} mutant crystals, median each;")
    print(f"  threshold calibrated to a {TARGET_FPR:.0%} null FPR.\n")
    print(f"  {'candidate':>20} {'SS':>12} {'slope':>7} {'noise':>7} "
          + " ".join(f"tau={t:<5.0f}" for t in TAUS))
    rows = []
    for c in cands:
        for ss in SS_MAIN:
            sl, sg = slopes[(c.name, ss)], noise[(c.name, ss)]
            if not np.isfinite(sl) or not np.isfinite(sg):
                print(f"  {c.name:>20} {ss:>12} {'abstains / inadmissible':>50}")
                continue
            cf = make_ceiling(sl, sg)
            cs = [cf(t) for t in TAUS]
            print(f"  {c.name:>20} {ss:>12} {sl:>7.3f} {sg:>7.2f} "
                  + " ".join(f"{v:>9.3f}" for v in cs))
            rows.append((c, ss, sl, sg, cs))
    print("\n  tau = SD of the TRUE mutation-induced full-span axis bend, in degrees.")

    print("\n  Inverted -- the number the real-data step has to beat:")
    print(f"  {'candidate':>20} {'SS':>12} {'tau* (AUC .75)':>15} {'tau* (.85)':>11} "
          f"{'Ca displacement':>17}")
    best = []
    for c, ss, sl, sg, _ in rows:
        cf = make_ceiling(sl, sg)
        got = {}
        for target in (0.75, 0.85):
            lo, hi = 0.02, 500.0
            for _ in range(22):
                mid = 0.5 * (lo + hi)
                v = cf(mid)
                if not np.isfinite(v) or v < target:
                    lo = mid
                else:
                    hi = mid
            got[target] = 0.5 * (lo + hi)
        L = IDEAL[ss][1] * (c.span - 1)          # axis arc length of the span
        disp = L * np.radians(got[0.75]) / 8.0   # sagitta ~ L*theta/8
        reach = got[0.75] < 480
        print(f"  {c.name:>20} {ss:>12} "
              f"{(f'{got[0.75]:.2f}' if reach else '>480'):>15} "
              f"{(f'{got[0.85]:.2f}' if got[0.85] < 480 else '>480'):>11} "
              f"{(f'{disp:.2f} A' if reach else 'unreachable'):>17}")
        best.append((c, ss, sl, sg, got[0.75], disp))
    print("\n  'Ca displacement' converts tau*(AUC 0.75) into the mid-span deflection")
    print("  of the axis (sagitta = L*theta/8) -- the physically interpretable form.")
    print("  Point mutations typically move local backbone by ~0.1-0.5 A, so a tau*")
    print("  whose displacement lands far above ~0.5 A is not reachable by the")
    print("  phenomenon under study.")
    return best


def verdict(best):
    hdr("V", "VERDICT")
    print(f"  coordinate error assumed : {SIGMA_REF} A  (Ca B-factor ~20-30 A^2)")
    print(f"  target null FPR          : {TARGET_FPR:.0%}")
    print(f"  aggregation              : {N_WT} WT / {N_MUT} mutant crystals, median")
    print(f"  plausibility bound       : Ca displacement <= 0.50 A\n")
    ok = [(c, ss, d) for c, ss, sl, sg, t75, d in best
          if c.is_axis_metric and d <= 0.50]
    if not ok:
        print("  MEASUREMENT BOTTLENECK -- no axis-bend candidate reaches oracle AUC")
        print("  0.75 within a physically plausible effect size. The project's answer")
        print("  is about measurement resolution, not about sequence. Do NOT resume")
        print("  the v1 ML gates. Options: restrict to a higher-resolution subset,")
        print("  require more crystals per variant, widen the span, or report the")
        print("  measurement limit itself as the result.")
    else:
        print("  GATE PASSED ON IDEAL GEOMETRY -- oracle AUC 0.75 is reachable within")
        print("  a plausible effect size, on exact/iid-perturbed ideal geometry, for:")
        for c, ss, d in sorted(ok, key=lambda x: x[2]):
            print(f"      {c.name:>20}  {ss:<12}  needs Ca displacement {d:.2f} A")
        print()
        print("  Next, in order:")
        print("   1. promote the winning candidate to metric v2 (it is already")
        print("      calibrated to report true full-span axis bend in degrees);")
        print("   2. build benchmark v2 -- exact RCSB entity/chain + SIFTS mapping,")
        print("      canonical sequence, exact variant, median aggregation over")
        print("      redundant structures, matched WT/mutant pairs, independent QC,")
        print("      family-level split. ONE BIOLOGICAL VARIANT = ONE OBSERVATION;")
        print("      bootstrap at variant/protein level (AUDIT.md F3);")
        print("   3. estimate tau empirically and re-read the stage-6 table -- that")
        print("      closes the loop this gate deliberately leaves open;")
        print("   4. only then baselines, each compared against the ceiling:")
        print("         SS-only  ->  local sequence  ->  direct-delta  ->  +3D context")
        print("      If SS-only already explains most of it, the 'local model skill'")
        print("      was structure-class prediction. Only a direct-delta model that")
        print("      fails while the ceiling is HIGH is a real negative ML result.")
    print("\n  NOT a statement that the label is usable on real data. Every number")
    print("  here is exact or iid-perturbed IDEAL geometry; real helices are pre-bent,")
    print("  frayed and irregular, with correlated non-Gaussian coordinate error, and")
    print("  tau is empirical. What this gate establishes is that the metric")
    print("  CONSTRUCTION is sound. The empirical seed test is what decides the label:")
    print("  20-50 curated variants, within-variant repeatability (iid vs systematic,")
    print("  stratified by crystal form), the observed WT->mutant delta distribution,")
    print("  tau, the null mover FPR, and the real-coordinate Jacobian from")
    print("  perturbation.py. Only then is a full RCSB/SIFTS miner worth building.")


def stage7():
    hdr(7, "SENSITIVITY -- what the verdict actually depends on")
    print("  The stage-6 verdict assumes sigma_xyz = 0.20 A and 5 WT / 2 mutant")
    print("  crystals. Both are optimistic in places, so both are swept here. Entry")
    print("  is tau* for oracle AUC 0.75, as a mid-span Ca displacement in Angstrom;")
    print("  '--' means unreachable. Plausibility bound is ~0.50 A.\n")
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(3)
    T = 40_000
    hyb = HybridAxisBend(9)
    combos = ((1, 1), (2, 1), (3, 1), (5, 2), (10, 3), (20, 5))

    for ss in ("alpha_helix", "beta_strand"):
        # calibrated slope ~1, so measure noise per sigma and invert
        print(f"  {ss}")
        print(f"      {'sigma_xyz':>10} " + " ".join(f"{a}WT/{b}mut".rjust(10)
                                                     for a, b in combos))
        r, rise, turn = IDEAL[ss]
        L = rise * (hyb.span - 1)
        for sg_xyz in (0.10, 0.15, 0.20, 0.30):
            P0 = helix_on_arc(hyb.span, r, rise, turn, 10.0)
            v = np.array([hyb(P0 + rng.normal(0, sg_xyz, P0.shape))
                          for _ in range(1500)])
            sig_m = float(v[np.isfinite(v)].std())
            cells = []
            for n_wt, n_mut in combos:
                agg = (np.median(rng.normal(0, sig_m, (T, n_mut)), axis=1)
                       - np.median(rng.normal(0, sig_m, (T, n_wt)), axis=1))
                null = (np.median(rng.normal(0, sig_m, (T, n_mut)), axis=1)
                        - np.median(rng.normal(0, sig_m, (T, n_wt)), axis=1))
                thr = float(np.quantile(np.abs(null), 1 - TARGET_FPR))
                z = rng.standard_normal(T)
                lo, hi = 0.02, 400.0
                for _ in range(20):
                    mid = 0.5 * (lo + hi)
                    t = mid * z
                    mov = np.abs(t + agg) > thr
                    if not (0 < mov.sum() < T):
                        lo = mid
                        continue
                    if roc_auc_score(mov, np.abs(t)) < 0.75:
                        lo = mid
                    else:
                        hi = mid
                tau = 0.5 * (lo + hi)
                d = L * np.radians(tau) / 8.0
                cells.append(f"{d:10.2f}" if tau < 380 else f"{'--':>10}")
            print(f"      {sg_xyz:>10.2f} " + " ".join(cells))
        print()
    print("  Reading: more crystals per variant buys resolution as 1/sqrt(n) -- but")
    print("  ONLY for the iid component simulated here. Real redundant crystals also")
    print("  share systematic effects (packing, space group, refinement protocol and")
    print("  era, construct, cryo vs RT, ligand state) that do NOT shrink with n, so")
    print("  the true curve flattens onto a floor this table does not contain.")
    print("  Treat every entry as an OPTIMISTIC bound. See CORRECTIONS.md C2.")
    print("  So the binding constraint on benchmark v2 is not the number of variants,")
    print("  it is the number of CRYSTALS PER VARIANT. That is the opposite of how v1")
    print("  was assembled: v1 spent its redundancy inflating n (248 rows from 69")
    print("  residues, AUDIT.md F3) instead of averaging it down into the label.")
    print("\n  Coverage cost, from v1's own committed CSV:")
    try:
        import collections
        import csv as _csv
        rows = list(_csv.DictReader(open("../feasibility_t4l.csv")))
        mult = collections.Counter(int(r["resnum"]) for r in rows)
        m = np.array(sorted(mult.values()))
        print(f"      {len(rows)} rows over {len(mult)} distinct residues; "
              f"median {int(np.median(m))} crystals/residue, mean {m.mean():.1f}")
        for k in (1, 2, 3, 5, 10):
            n = int((m >= k).sum())
            print(f"      residues with >= {k:>2} crystals: {n:>3} / {len(mult)} "
                  f"({n/len(mult):.0%})")
        print("\n      So a 5-WT/2-mutant requirement is met by only a minority of")
        print("      variants. v2 must report this as coverage, and weight or stratify")
        print("      by per-variant crystal count rather than pretending it is uniform.")
    except FileNotFoundError:
        print("      (feasibility_t4l.csv not found; run from the v2/ directory)")


STAGES = {"1": stage1, "2": stage2, "3": stage3, "4": stage4, "5": stage5,
          "7": stage7}

if __name__ == "__main__":
    want = sys.argv[1:] or ["1", "2", "3", "4", "5", "6", "7"]
    noise = slopes = None
    for w in want:
        if w == "6":
            continue
        if w == "7":
            continue
        if w not in STAGES:
            sys.exit(f"unknown stage {w}; available 1-7")
        r = STAGES[w]()
        if w == "2":
            slopes = {k: v[0] for k, v in r.items()}
        if w == "3":
            noise = r
    if "6" in want:
        verdict(stage6(noise, slopes))
    if "7" in want:
        stage7()
    print()
