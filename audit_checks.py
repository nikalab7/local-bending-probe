"""
Verification harness for AUDIT.md.

Every quantitative claim in AUDIT.md is reproduced here so it can be checked
rather than taken on trust. Runs offline (no PDB downloads); the only repo data
it reads is the committed `feasibility_t4l.csv`.

    python3 audit_checks.py            # all checks
    python3 audit_checks.py A1 B3      # selected checks

Requires numpy, scipy, scikit-learn. Most checks are seconds; A2, B4, B5 and C2
are Monte-Carlo and take a few minutes each (C2 is the slowest -- it runs a
nested bootstrap). Seeds are fixed, so figures are reproducible run to run;
the Monte-Carlo ones are quoted in AUDIT.md to the precision they are stable at.

Checks
------
A1  metric reads 110 deg on a straight ideal alpha-helix, 0 deg on a beta-strand
A2  the label is separable by secondary-structure class alone
A3  principal-axis degeneracy / sign flip at Ca pseudo-angle 60 deg
A4  the module self-test never exercises the principal-axis fit
A5  metric noise vs coordinate error, against the adopted 0.98 deg floor
B1  bias of sigma_hat = 1.4826 * MAD at small n
B2  null false-positive rate of the |delta| > 2*sigma_hat mover test
B3  bias of the "independent" np.std(vals, ddof=1) floor at n=2
B4  the achievable AUC ceiling: what a PERFECT predictor could score
B5  positive control -- can the Gate-2 architecture recover a planted local effect?
C1  pseudo-replication in feasibility_t4l.csv; Kish effective n
C2  row bootstrap vs cluster bootstrap coverage under the null
D1  is_continuous accepts physically impossible windows
D2  window-range off-by-one between Gate 1 and Gates 2-5
D3  the documented cross-check metric contradicts the metric it checks
"""
from __future__ import annotations
import sys
import numpy as np

from bending_metric import bending_angle, is_continuous, _principal_axis, _CA_CA_IDEAL

D_CA = 3.8
ADOPTED_FLOOR = 0.98      # deg, RESULTS.md Gate 1
CORROB_FLOOR = 0.75       # deg, RESULTS.md Gate 1 "independent" estimate
MOVER_MEDIAN = 2.66       # deg, gate2_model_feasibility.ABOVE_FLOOR_SIGNAL


def hdr(tag, title):
    print()
    print("=" * 78)
    print(f"  {tag}. {title}")
    print("=" * 78)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def coil(n=5, r=2.3, rise=1.5, turn=100.0, phase=0.0):
    """Ideal helical Ca trace. (r,rise,turn) = (2.3, 1.5, 100) is the alpha-helix."""
    return np.array([[r * np.cos(np.radians(turn * i + phase)),
                      r * np.sin(np.radians(turn * i + phase)),
                      rise * i] for i in range(n)])


def sym_triple(phi_deg, d=D_CA):
    """Three Ca with apex angle phi at the middle atom; symmetric, bond length d."""
    h = np.radians(phi_deg) / 2
    return np.vstack([[-d * np.sin(h), -d * np.cos(h), 0.0],
                      [0.0, 0.0, 0.0],
                      [d * np.sin(h), -d * np.cos(h), 0.0]])


def bent_window(theta_deg, scale=D_CA):
    """The window shape used by bending_metric.py's own self-test."""
    t = np.radians(theta_deg)
    d1 = np.array([1.0, 0.0, 0.0])
    d2 = np.array([np.cos(t), np.sin(t), 0.0])
    return np.vstack([-2 * d1, -d1, np.zeros(3), d2, 2 * d2]) * scale


def pseudo_angles(P):
    out = []
    for k in (1, 2, 3):
        u, v = P[k - 1] - P[k], P[k + 1] - P[k]
        out.append(np.degrees(np.arccos(np.clip(
            np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))
    return np.array(out)


def exterior_angle_sum(P):
    """The 'cross-check metric' described in bending_metric.py's docstring."""
    return float((180.0 - pseudo_angles(P)).sum())


# --------------------------------------------------------------------------- #
# A. the metric
# --------------------------------------------------------------------------- #
def A1():
    hdr("A1", "What does the metric read on ideal, UNBENT secondary structure?")
    for name, P in [("ideal alpha-helix", coil()),
                    ("ideal beta-strand", coil(r=0.95, rise=3.3, turn=180.0))]:
        ca = np.linalg.norm(np.diff(P, axis=0), axis=1)
        print(f"  {name}:  Ca-Ca = {ca.round(2)}   pseudo-angles = {pseudo_angles(P).round(1)}")
        print(f"      bending_angle = {bending_angle(P):7.2f} deg")
    print("\n  Both have a perfectly STRAIGHT axis. The docstring self-test advertises")
    print("  'straight = 0.00 deg'. A straight alpha-helix reads ~110 deg.")

    print("\n  robustness over the plausible parameter range:")
    for nm, rr, ri, tt in [("alpha-helix", (2.2, 2.26, 2.3, 2.4), (1.4, 1.5, 1.6),
                            (96, 99, 100, 102, 104)),
                           ("beta-strand", (0.9, 0.95, 1.0), (3.2, 3.3, 3.4),
                            (170, 175, 180, 185, 190))]:
        v = np.array([bending_angle(coil(5, r, ri_, t, ph))
                      for r in rr for ri_ in ri for t in tt
                      for ph in (0, 37, 111, 203)])
        print(f"      {nm}: n={len(v):4d}  bending = {v.mean():7.2f} +/- {v.std():.2f} deg"
              f"   range [{v.min():.2f}, {v.max():.2f}]")
    print("\n  => the metric separates helix from strand by ~110 deg. It is dominated by")
    print("     secondary-structure TYPE, not by curvature of the backbone axis.")


def A2():
    hdr("A2", "How much of the label is explained by secondary-structure class alone?")
    rng = np.random.default_rng(7)

    def loop_window():
        P = [np.zeros(3), np.array([D_CA, 0.0, 0.0])]
        for _ in range(3):
            prev = P[-1] - P[-2]
            prev = prev / np.linalg.norm(prev)
            a = np.radians(180 - rng.uniform(80, 145))
            az = rng.uniform(0, 2 * np.pi)
            t = np.cross(prev, [0, 0, 1.0])
            if np.linalg.norm(t) < 1e-6:
                t = np.cross(prev, [0, 1.0, 0])
            t /= np.linalg.norm(t)
            u = np.cross(prev, t)
            d = np.cos(a) * prev + np.sin(a) * (np.cos(az) * t + np.sin(az) * u)
            P.append(P[-1] + D_CA * d)
        return np.array(P)

    n, noise = 4000, 0.25
    groups = {}
    groups["H"] = np.array([bending_angle(coil(5, rng.uniform(2.2, 2.4), rng.uniform(1.4, 1.6),
                                               rng.uniform(96, 104), rng.uniform(0, 360))
                                         + rng.normal(0, noise, (5, 3))) for _ in range(n)])
    groups["E"] = np.array([bending_angle(coil(5, rng.uniform(0.9, 1.0), rng.uniform(3.2, 3.4),
                                               rng.uniform(170, 190), rng.uniform(0, 360))
                                         + rng.normal(0, noise, (5, 3))) for _ in range(n)])
    groups["L"] = np.array([bending_angle(loop_window() + rng.normal(0, noise, (5, 3)))
                            for _ in range(n)])
    for k, nm in [("H", "helix"), ("E", "strand"), ("L", "loop")]:
        print(f"  {nm:>6} windows: mean {groups[k].mean():7.2f}  sd {groups[k].std():5.2f} deg")
    y = np.concatenate(list(groups.values()))
    resid = np.concatenate([g - g.mean() for g in groups.values()])
    print(f"\n  total SD of label                = {y.std():6.2f} deg")
    print(f"  SD after removing SS-class mean  = {resid.std():6.2f} deg")
    print(f"  => a 3-class SS-only predictor achieves {100*(1-resid.std()/y.std()):.0f}% RMSE "
          "reduction on this mixture,")
    print("     knowing nothing about bending. Gate 2 reports 20% (30.55 vs 37.97 deg).")
    print("     A 20% reduction is therefore not evidence of learned bending; it is within")
    print("     what secondary-structure propensity alone supplies.")


def A3():
    hdr("A3", "Principal-axis degeneracy of a 3-point half (bending_metric.py:42-53)")
    print("  For a symmetric triple with apex angle phi and bond length d:")
    print("      var along chord    = (2/3) d^2 sin^2(phi/2)")
    print("      var perpendicular  = (2/9) d^2 cos^2(phi/2)")
    print("  equal when tan(phi/2) = 1/sqrt(3)  =>  phi = 60 deg exactly.")
    print(f"\n{'phi':>6} {'s0':>8} {'s1':>8} {'s1/s0':>8} {'angle(axis,chord)':>19} {'dot(axis,p2-p0)':>17}")
    for phi in (180, 130, 120, 110, 91, 80, 75, 70, 65, 61, 60, 59, 55, 45):
        P = sym_triple(phi)
        _, s, _ = np.linalg.svd(P - P.mean(axis=0), full_matrices=False)
        ax = _principal_axis(P)
        chord = (P[-1] - P[0]) / np.linalg.norm(P[-1] - P[0])
        ang = np.degrees(np.arccos(np.clip(abs(np.dot(ax, chord)), -1, 1)))
        print(f"{phi:>6} {s[0]:>8.4f} {s[1]:>8.4f} {s[1]/max(s[0],1e-30):>8.4f} "
              f"{ang:>19.2f} {np.dot(ax, P[-1]-P[0]):>17.4f}")
    print("\n  Two distinct failures:")
    print("   (i) below phi = 60 deg the axis is PERPENDICULAR to the chain, and the")
    print("       orientation test at bending_metric.py:51 (sign of dot(axis, p2-p0))")
    print("       has a vanishing argument -> a ~180 deg flip decided by rounding.")
    print("  (ii) s1/s0 -> 1 as phi -> 60 deg. At phi = 70-80 deg (tight turns, i.e.")
    print("       exactly the loop population Gates 3-4 target) s1/s0 = 0.69-0.82, so")
    print("       the axis direction is poorly conditioned and noise-amplifying.")
    print("  Physical Ca pseudo-angles bottom out near 75-80 deg, so the hard flip is")
    print("  rarely reached -- but the ill-conditioning is reached routinely, in loops.")


def A4():
    hdr("A4", "Does the module self-test exercise the principal-axis fit it justifies?")
    print("  bending_metric.py:99-103 builds each test window from two straight arms.")
    print(f"\n{'theta':>8} {'s1/s0 half 1':>14} {'s1/s0 half 2':>14}")
    for th in (0.0, 30.0, 90.0, 135.0, 160.0, 170.0):
        P = bent_window(th)
        r = []
        for half in (P[0:3], P[2:5]):
            _, s, _ = np.linalg.svd(half - half.mean(axis=0), full_matrices=False)
            r.append(s[1] / s[0])
        print(f"{th:>8.1f} {r[0]:>14.2e} {r[1]:>14.2e}")
    print("\n  Every half is EXACTLY collinear (s1 = 0), so the SVD returns the chord")
    print("  trivially. The docstring justifies the PCA fit as 'more robust than an")
    print("  endpoint vector', but the test only ever probes the arccos. Substituting")
    print("  axis := P[2]-P[0] would pass the entire self-test unchanged:")
    for th in (0.0, 90.0, 170.0):
        P = bent_window(th)
        d1 = P[2] - P[0]; d2 = P[4] - P[2]
        naive = np.degrees(np.arccos(np.clip(np.dot(d1, d2) /
                                             np.linalg.norm(d1) / np.linalg.norm(d2), -1, 1)))
        print(f"      theta={th:6.1f}:  PCA {bending_angle(P):7.3f}   endpoint-chord {naive:7.3f}")


def A5():
    hdr("A5", "Metric noise vs Ca coordinate error, against the adopted noise floor")
    rng = np.random.default_rng(0)
    confs = [("ideal alpha-helix", coil()),
             ("ideal beta-strand", coil(r=0.95, rise=3.3, turn=180.0)),
             ("collinear (self-test 'straight')", bent_window(0.0)),
             ("bent 30 deg", bent_window(30.0)),
             ("bent 90 deg", bent_window(90.0))]
    sigmas = (0.05, 0.10, 0.20, 0.30)
    print(f"{'conformation':>34} " + " ".join(f"{s:>7.2f}A" for s in sigmas))
    for name, P0 in confs:
        row = [np.std([bending_angle(P0 + rng.normal(0, sg, P0.shape)) for _ in range(3000)])
               for sg in sigmas]
        print(f"{name:>34} " + " ".join(f"{v:>8.2f}" for v in row))
    print("\n  Adopted per-window floor (RESULTS.md Gate 1) : "
          f"{ADOPTED_FLOOR:.2f} deg")
    print(f"  'Independent' corroboration                  : {CORROB_FLOOR:.2f} deg")
    print(f"  Median |delta| of the 72 'movers'            : {MOVER_MEDIAN:.2f} deg")
    print("\n  A 2.0 A X-ray structure carries ~0.1-0.3 A positional error on a")
    print("  well-ordered Ca. Propagated through the metric that is ~2-7 deg for a")
    print("  helical window -- larger than the adopted floor AND larger than the")
    print("  median effect the project calls signal.")


# --------------------------------------------------------------------------- #
# B. the noise floor and the mover test
# --------------------------------------------------------------------------- #
def B1(T=200_000):
    hdr("B1", "Bias of sigma_hat = 1.4826 * MAD at small n  (truth: sigma = 1)")
    rng = np.random.default_rng(12345)
    print(f"{'n':>5} {'E[sigma_hat]':>14} {'bias':>10}")
    for n in (3, 4, 5, 6, 8, 10, 15, 20, 50):
        v = rng.standard_normal((T, n))
        sig = 1.4826 * np.median(np.abs(v - np.median(v, axis=1)[:, None]), axis=1)
        print(f"{n:>5} {sig.mean():>14.4f} {sig.mean()-1:>+10.1%}")
    print("\n  feasibility_t4l.py:177 admits a window at len(vals) >= 3, where the")
    print("  estimator is biased low by a third. Every gate reuses this estimator.")
    print("  A smaller floor means more 'movers'; the bias is one-directional.")


def B2(T=200_000):
    hdr("B2", "Null false-positive rate of the mover test")
    from scipy.stats import norm
    rng = np.random.default_rng(999)
    print("  Test as implemented:  |b - median(WT_1..WT_n)| > 2 * 1.4826*MAD(WT_1..WT_n)")
    print("  Null: the mutation has NO geometric effect (b and WT_i iid normal).\n")
    print(f"{'n_WT':>6} {'null FPR':>10} {'expected FP of 248':>20} {'true sigma-multiple':>21}")
    for n in (3, 4, 5, 6, 8, 10, 15, 20, 50):
        v = rng.standard_normal((T, n))
        b = rng.standard_normal(T)
        sig = 1.4826 * np.median(np.abs(v - np.median(v, axis=1)[:, None]), axis=1)
        fpr = float(np.mean(np.abs(b - np.median(v, axis=1)) > 2 * sig))
        print(f"{n:>6} {fpr:>10.4f} {fpr*248:>20.1f} {norm.isf(fpr/2):>21.2f}")
    print(f"\n  An exact 2-sigma two-sided test with sigma KNOWN would give "
          f"{2*norm.sf(2):.4f} ({2*norm.sf(2)*248:.1f} of 248).")
    print("\n  Two compounding errors, both in the same direction:")
    print("   (i) sigma_hat is biased low (B1).")
    print("  (ii) the reference is a median of n, so SD(b - med) = sigma*sqrt(1+pi/2n),")
    print("       not sigma. The 2-sigma cut is fewer than 2 true SDs:")
    for n in (3, 5, 10, 20):
        s_d = np.sqrt(1 + np.pi / (2 * n))
        print(f"          n={n:>2}: SD = {s_d:.3f} sigma -> cut is {2/s_d:.2f} true SDs "
              f"-> FPR {2*norm.sf(2/s_d):.4f}")
    print("\n  feasibility_t4l.py:12-14 justifies the 2-sigma cut as clearing a")
    print("  'differencing noise ~ sqrt(2)*sigma' band. That is the noise of two SINGLE")
    print("  measurements; the code differences a single mutant against a MEDIAN. The")
    print("  stated derivation does not describe the implemented test.")
    print("\n  RESULTS.md reports 29% of mutations above floor. At the admitted minimum")
    print("  of n=3 WT crystals the null alone yields ~41%. No gate computes this null.")


def B3(T=200_000):
    hdr("B3", "The 'independent' floor: np.std(vals, ddof=1)  (feasibility_t4l.py:216)")
    rng = np.random.default_rng(5)
    print("  Admitted at len(vals) >= 2. ddof=1 is unbiased for the VARIANCE, not the SD.")
    for n in (2, 3, 4, 5):
        s = rng.standard_normal((T, n)).std(axis=1, ddof=1)
        print(f"      n={n}: E[s] = {s.mean():.4f} sigma   (bias {s.mean()-1:+.1%})")
    print(f"      n=2 analytic: E[s] = sqrt(2/pi) = {np.sqrt(2/np.pi):.4f}")
    print("\n  and the reported figure is the MEDIAN over such groups, which biases further:")
    for g in (5, 10, 20):
        med = np.median(rng.standard_normal((T // 20, g, 2)).std(axis=2, ddof=1), axis=1)
        print(f"      median of {g:>2} two-crystal estimates: E = {med.mean():.4f} sigma "
              f"(bias {med.mean()-1:+.1%})")
    print(f"\n  RESULTS.md treats {ADOPTED_FLOOR:.2f} deg and {CORROB_FLOOR:.2f} deg as two")
    print("  independent estimates that 'agree', and reads the agreement as validation.")
    print(f"  Their ratio is {CORROB_FLOOR/ADOPTED_FLOOR:.2f}; the ratio of the two")
    print("  estimators' biases predicts roughly that. Two estimators of the same sigma,")
    print("  biased low by different known amounts, are expected to disagree in exactly")
    print("  this direction. Their 'agreement' is not evidence that either is correct.")


def B4(N=400_000):
    hdr("B4", "The achievable AUC ceiling -- what could a PERFECT predictor score?")
    from sklearn.metrics import roc_auc_score
    from scipy.stats import norm
    from scipy.optimize import brentq
    rng = np.random.default_rng(0)
    thr = 2 * ADOPTED_FLOOR
    print(f"  Mover labels come from a NOISY observation:  o = t + e")
    print(f"      t = true mutation effect (SD tau)    e = measurement noise (SD s)")
    print(f"      mover  <=>  |o| > {thr:.2f} deg      (the project's 2-sigma cut)")
    print(f"  An ORACLE knowing t exactly scores |t|. Its AUC bounds ANY predictor.\n")
    print(f"{'tau':>7} {'s':>7} {'base rate':>11} {'ORACLE AUC':>12}")
    for tau in (1.0, 2.0, 3.0, 5.0):
        for s in (0.98, 2.0, 3.0, 4.5):
            t = rng.normal(0, tau, N)
            mov = np.abs(t + rng.normal(0, s, N)) > thr
            if 0 < mov.sum() < N:
                print(f"{tau:>7.1f} {s:>7.2f} {mov.mean():>11.1%} "
                      f"{roc_auc_score(mov, np.abs(t)):>12.3f}")

    print("\n  Now constrain to the project's OWN reported base rate of 29%:")
    for s in (0.98, 1.5, 2.0, 3.0, 4.5):
        best = None
        for tau in np.arange(0.0, 6.01, 0.04):
            br = np.mean(np.abs(rng.normal(0, tau, 40000) +
                                rng.normal(0, s, 40000)) > thr)
            if best is None or abs(br - 0.29) < abs(best[1] - 0.29):
                best = (tau, br)
        tau = best[0]
        t = rng.normal(0, tau, N)
        mov = np.abs(t + rng.normal(0, s, N)) > thr
        auc = roc_auc_score(mov, np.abs(t)) if 0 < mov.sum() < N else float("nan")
        print(f"      noise s={s:>4.2f} -> implied true-effect SD tau={tau:>4.2f} deg "
              f"(base {mov.mean():.0%});  ORACLE AUC = {auc:.3f}")

    s0 = brentq(lambda s: 2 * norm.sf(thr / s) - 0.29, 0.5, 20)
    print(f"\n  And with NO true effect at all (tau = 0), the base rate is 29% when")
    print(f"  the measurement noise alone has SD = {s0:.3f} deg -- only "
          f"{s0/ADOPTED_FLOOR:.2f}x the adopted floor,")
    print(f"  and inside the {2}-{7} deg range that coordinate-error propagation gives (A5):")
    for sg in (0.98, 1.5, 1.85, 2.0, 3.0, 4.5):
        print(f"      pure noise SD {sg:>4.2f} deg, zero true effect -> base rate "
              f"{2*norm.sf(thr/sg):>5.1%}")
    print("\n  Reported AUCs: 0.477 / 0.516 / 0.524 / 0.590 / 0.591 / 0.700.")
    print("  If the noise is ~0.98 deg as claimed, an oracle scores ~0.84 and the")
    print("  measured ~0.52 is a real model failure. If the noise is >=1.85 deg, the")
    print("  labels carry no recoverable signal, the ceiling IS ~0.50, and every")
    print("  reported AUC is at ceiling. The repo never distinguishes these cases.")


def B5():
    hdr("B5", "Positive control: can the Gate-2 architecture recover a planted effect?")
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(0)
    N = 60000
    geo = np.column_stack([rng.uniform(80, 150, N), rng.uniform(80, 150, N),
                           rng.uniform(-1, 1, N), rng.uniform(-1, 1, N)])
    seq = rng.integers(0, 20, (N, 5))
    oh = np.zeros((N, 100))
    for p in range(5):
        oh[np.arange(N), p * 20 + seq[:, p]] = 1.0
    X = np.hstack([geo, oh])
    eff = np.linspace(-1.5, 1.5, 20)          # planted central-residue effect, 3 deg p2p
    y = 0.9 * geo[:, 0] + 0.4 * geo[:, 1] + 8 * geo[:, 2] + eff[seq[:, 2]] \
        + rng.normal(0, 5, N)
    print("  Same features (4 geometry + 5x20 one-hot), same estimator settings,")
    print("  same geometry-dominated label. A LOCAL central-residue effect of")
    print(f"  {np.ptp(eff):.1f} deg peak-to-peak is planted by construction.\n")
    m = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.08, max_depth=6,
                                      random_state=0).fit(X, y)
    K = 248
    idx = rng.choice(N, K, replace=False)
    Xwt = X[idx].copy()
    Xmt = Xwt.copy()
    newc = rng.integers(0, 20, K)
    Xmt[:, 44:64] = 0.0
    Xmt[np.arange(K), 44 + newc] = 1.0
    dpred = m.predict(Xmt) - m.predict(Xwt)
    true_d = eff[newc] - eff[seq[idx, 2]]
    print(f"  exact zeros in predicted delta : {np.mean(dpred == 0):.1%}")
    print(f"  distinct |dpred| values        : "
          f"{len(np.unique(np.round(np.abs(dpred), 12)))} of {K}")
    print(f"  SD predicted / SD true delta   : {dpred.std():.3f} / {true_d.std():.3f} "
          f"= {dpred.std()/true_d.std():.2f}x")
    mov = np.abs(true_d) > np.percentile(np.abs(true_d), 71)   # top 29%, as observed
    print(f"  retrieval of planted movers    : AUC = "
          f"{roc_auc_score(mov, np.abs(dpred)):.3f}")
    print("\n  The architecture PASSES this control: predict-absolute-then-difference")
    print("  does recover a genuine 3 deg local effect, and the differenced prediction")
    print("  is neither degenerate nor heavily tied. So 'the estimator ignored the")
    print("  one-hot block' is NOT the explanation for the observed AUC ~ 0.5 --")
    print("  a check worth running, and it comes out in the project's favour.")
    print("  Caveat: here the mover label is noise-free. With labels as noisy as B4")
    print("  suggests, no predictor could reach this AUC. That is the live confound,")
    print("  and it is a property of the LABELS, not of the model.")


# --------------------------------------------------------------------------- #
# C. pseudo-replication
# --------------------------------------------------------------------------- #
def _load_csv():
    import csv
    rows = list(csv.DictReader(open("feasibility_t4l.csv")))
    return (np.array([float(r["delta_deg"]) for r in rows]),
            np.array([float(r["floor_sigma_deg"]) for r in rows]),
            np.array([int(r["resnum"]) for r in rows]))


def C1():
    hdr("C1", "Pseudo-replication in the committed feasibility_t4l.csv")
    import collections
    d, f, res = _load_csv()
    mov = np.abs(d) > 2 * f
    mult = collections.Counter(res)
    m = np.array(list(mult.values()))
    print(f"  rows (called 'pairs' / 'windows')      : {len(d)}")
    print(f"  distinct mutated residues              : {len(mult)}")
    print(f"  most-replicated residue                : {max(mult, key=mult.get)} "
          f"x{max(mult.values())}")
    print(f"  movers                                 : {mov.sum()} ({mov.mean():.1%})")
    print("  top multiplicities                     : "
          + ", ".join(f"res {int(r)} x{c}" for r, c in mult.most_common(8)))
    neff = m.sum() ** 2 / (m ** 2).sum()
    print(f"\n  Kish effective n = (sum m)^2 / sum m^2 = {neff:.1f}   (nominal {len(d)})")
    print(f"  variance inflation                     ~ {len(d)/neff:.1f}x")
    print(f"  CIs computed on rows are too narrow by ~ {np.sqrt(len(d)/neff):.1f}x")
    print(f"\n  distinct per-window floor values       : {len(set(np.round(f,9)))}")
    print("  -> rows sharing a residue share the SAME wt_med and the SAME floor, so")
    print("     they are not independent draws even before counting repeat crystals.")
    print("\n  SPEC_bending_and_pairs.md section 2.4 requires:")
    print('     "aggregate to ONE bending per variant (median across redundant')
    print('      structures); keep the spread as the per-variant noise estimate"')
    print("  No script implements it. Applying it to this CSV, per residue:")
    uniq = np.array(sorted(set(res)))
    per = [abs(np.median(d[res == r])) > 2 * np.median(f[res == r]) for r in uniq]
    print(f"     per-residue mover fraction = {np.mean(per):.1%} ({sum(per)}/{len(uniq)})"
          f"   vs the reported {mov.mean():.1%}")

    rng = np.random.default_rng(3)
    by = {r: np.where(res == r)[0] for r in uniq}
    naive = [mov[rng.integers(0, len(mov), len(mov))].mean() for _ in range(20000)]
    clus = []
    for _ in range(20000):
        pick = rng.integers(0, len(uniq), len(uniq))
        clus.append(mov[np.concatenate([by[uniq[p]] for p in pick])].mean())
    nlo, nhi = np.percentile(naive, [5, 95])
    clo, chi = np.percentile(clus, [5, 95])
    print(f"\n  90% CI on the mover fraction:")
    print(f"      row bootstrap (as used)  [{nlo:.3f}, {nhi:.3f}]  width {nhi-nlo:.3f}")
    print(f"      cluster bootstrap        [{clo:.3f}, {chi:.3f}]  width {chi-clo:.3f}"
          f"   ({(chi-clo)/(nhi-nlo):.2f}x wider)")


def C2(T=120, B=250):
    hdr("C2", "Null coverage of the row bootstrap under the real duplication pattern")
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(11)
    _, _, res = _load_csv()
    uniq = np.array(sorted(set(res)))
    by = {r: np.where(res == r)[0] for r in uniq}
    n = len(res)
    print("  Note on why the score is shared within a residue: the differenced")
    print("  prediction depends only on (scaffold entry features, wt_seq, mut_seq).")
    print("  For two crystals of the SAME variant all three are identical, so dpred is")
    print("  bit-identical -- every repeat crystal is an exact tie in the ROC.")
    print("\n  Simulating a TRUE NULL (no signal) with that structure; a well-calibrated")
    print("  90% CI should exclude 0.50 in 10% of runs.\n")
    for kind in ("row", "cluster"):
        excl = tot = 0
        for _ in range(T):
            lab = np.zeros(n, bool)
            sc = np.empty(n)
            for r in uniq:
                lab[by[r]] = rng.random() < 0.29
                sc[by[r]] = rng.standard_normal()
            if not (0 < lab.sum() < n):
                continue
            b = []
            for _ in range(B):
                if kind == "row":
                    idx = rng.integers(0, n, n)
                else:
                    pick = rng.integers(0, len(uniq), len(uniq))
                    idx = np.concatenate([by[uniq[p]] for p in pick])
                if 0 < lab[idx].sum() < len(idx):
                    b.append(roc_auc_score(lab[idx], sc[idx]))
            if len(b) < 50:
                continue
            lo, hi = np.percentile(b, [5, 95])
            tot += 1
            excl += (lo > 0.5 or hi < 0.5)
        print(f"      {kind:>8} bootstrap: excluded 0.50 in {100*excl/max(tot,1):5.1f}% "
              f"of {tot} null runs   (nominal 10%)")
    print("\n  loop_gate.py:260-265, powered_loop_gate.py:198-203 and gate3_3d.py:148-153")
    print("  all resample ROWS. Under the null that is strongly anti-conservative, so the")
    print("  published intervals -- [0.49,0.89], [0.49,0.69], [0.52,0.66] -- are narrower")
    print("  than the data support. This cuts AGAINST the project's own conclusions too:")
    print("  the true intervals are wider, so 'not distinguishable from chance' is if")
    print("  anything understated, while Gate 5's 'lower bound barely clears 0.50' is not")
    print("  supportable at all.")


# --------------------------------------------------------------------------- #
# D. localized correctness defects
# --------------------------------------------------------------------------- #
def D1():
    hdr("D1", "is_continuous has no lower bound (bending_metric.py:87-93)")
    cases = [("duplicated Ca (0.0 A step)",
              np.array([[0, 0, 0], [3.8, 0, 0], [3.8, 0, 0], [7.6, 0, 0], [11.4, 0, 0]], float)),
             ("0.3 A steps (physically impossible)",
              np.array([[0, 0, 0], [.3, 0, 0], [.6, 0, 0], [.9, 0, 0], [1.2, 0, 0]], float)),
             ("cis-peptide-like 2.9 A (legitimate)",
              np.array([[0, 0, 0], [2.9, 0, 0], [6.7, 0, 0], [10.5, 0, 0], [14.3, 0, 0]], float))]
    for nm, P in cases:
        try:
            b = f"{bending_angle(P):.2f} deg"
        except Exception as e:
            b = f"raised {type(e).__name__}"
        print(f"  {nm:>38}: is_continuous={str(is_continuous(P)):>5}   bending={b}")
    print("\n  Only an UPPER bound (4.5 A) is checked. A duplicated or badly-placed Ca")
    print("  passes and yields a confident 0.00 deg. A lower bound near 2.5 A would")
    print("  keep genuine cis-peptides (~2.9 A) and reject the rest.")
    print(f"\n  _CA_CA_IDEAL = {_CA_CA_IDEAL} is defined and never used anywhere in the repo.")


def D2():
    hdr("D2", "Window-range off-by-one between Gate 1 and Gates 2-5")
    lo, hi = 1, 50
    print(f"  residues {lo}..{hi}; the last window that fits (i..i+4) starts at i = {hi-4}")
    print(f"\n  feasibility_t4l.py:172            range(min, max-4+1) -> last start "
          f"{list(range(lo, hi-4+1))[-1]}   correct")
    for f_ in ("gate2_model_feasibility.py:162", "loop_gate.py:121",
               "mover_composition.py:64", "gate3_3d.py:111", "powered_loop_gate.py:111"):
        print(f"  {f_:<33} range(min, max-4)   -> last start "
              f"{list(range(lo, hi-4))[-1]}   drops window {hi-4}")
    print("\n  Gate 1 builds wt_med/wt_sig over one more window than Gates 2-5 do. The")
    print("  pooled floor is a median over that dictionary, so the 'same' floor is not")
    print("  numerically the same between gates -- one contributor to the 0.524 / 0.516")
    print("  drift that RESULTS.md attributes wholly to the QC step.")
    print("  featurize_chain / build_training use range(lo+3, hi-4) and drop it too.")


def D3():
    hdr("D3", "The documented cross-check metric contradicts the metric it checks")
    print("  bending_metric.py:31-33: 'Cross-check metric (optional): unsigned sum of the")
    print("  three C-alpha pseudo-bond exterior angles. Should track this within a few")
    print("  degrees; large disagreement flags an S-shaped window.'\n")
    print(f"{'conformation':>34} {'bending':>10} {'ext-angle sum':>15} {'disagreement':>13}")
    for nm, P in [("ideal alpha-helix", coil()),
                  ("ideal beta-strand", coil(r=0.95, rise=3.3, turn=180.0)),
                  ("collinear ('straight')", bent_window(0.0)),
                  ("bent 90 deg", bent_window(90.0)),
                  ("bent 170 deg", bent_window(170.0))]:
        b, e = bending_angle(P), exterior_angle_sum(P)
        print(f"{nm:>34} {b:>10.2f} {e:>15.2f} {abs(b-e):>13.2f}")
    print("\n  The two quantities agree only on the collinear windows of the self-test.")
    print("  On real secondary structure they differ by 90-180 deg. The cross-check is")
    print("  also never implemented in code, so the stated safeguard against S-shaped")
    print("  (sign-cancelling) windows does not exist.")


CHECKS = {"A1": A1, "A2": A2, "A3": A3, "A4": A4, "A5": A5,
          "B1": B1, "B2": B2, "B3": B3, "B4": B4, "B5": B5,
          "C1": C1, "C2": C2, "D1": D1, "D2": D2, "D3": D3}

if __name__ == "__main__":
    want = [a.upper() for a in sys.argv[1:]] or list(CHECKS)
    bad = [w for w in want if w not in CHECKS]
    if bad:
        sys.exit(f"unknown check(s): {bad}; available: {list(CHECKS)}")
    for w in want:
        CHECKS[w]()
    print()
