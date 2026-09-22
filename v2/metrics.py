"""
Candidate local-geometry metrics for the v2 benchmark.

Why this file exists
--------------------
The v1 metric (5-Ca PCA half-axis angle) reads 110 deg on a *straight* ideal
alpha-helix and 0 deg on a *straight* ideal beta-strand -- see AUDIT.md F1. It is
therefore dominated by secondary-structure type rather than by curvature of the
backbone axis. Before any model is trained again, the metric itself has to be
chosen empirically against stated criteria.

All candidates share one signature so the gate can treat them uniformly:

    m = Candidate(...)
    m.span            -> number of consecutive Ca required
    m(ca)             -> float, degrees  (a turning angle)
    m.is_axis_metric  -> True if it claims to measure AXIS bend (so a straight
                         helix/strand must read ~0)

Ground truth
------------
`helix_on_arc` winds a regular helix (or strand, or any (r, rise, turn)) around an
axis that is an exact circular arc of known total turning. That gives an analytic
"true axis bend" to calibrate against, which no synthetic construction in v1 had:
v1's self-test used only windows whose halves were exactly collinear, so it could
not distinguish these candidates at all (AUDIT.md F12).
"""
from __future__ import annotations

import numpy as np

_EPS = 1e-12

# Textbook backbone parameters (radius A, rise per residue A, twist deg/residue).
IDEAL = {
    "alpha_helix": (2.30, 1.50, 100.0),
    "beta_strand": (0.95, 3.30, 180.0),
    "pi_helix":    (2.63, 1.15,  87.0),
    "three_ten":   (1.90, 2.00, 120.0),
    "pp_ii":       (1.28, 3.12, 120.0),   # polyproline II, common in loops
}


# --------------------------------------------------------------------------- #
# ground truth: a regular helix wound around a circular-arc axis
# --------------------------------------------------------------------------- #
def helix_on_arc(n, r=2.30, rise=1.50, turn=100.0, total_bend_deg=0.0, phase=0.0):
    """n Ca wound at radius `r` around an axis that is a circular arc.

    The axis tangent turns by exactly `total_bend_deg` across the whole span, so
    `total_bend_deg` is the ground-truth axis bend. 0 => perfectly straight axis.

    Returns (n, 3). Ca-Ca spacing is set by (r, rise, turn) and is ~3.8 A for the
    IDEAL presets.
    """
    L = rise * (n - 1)                      # axis arc length spanned
    kappa = np.radians(total_bend_deg) / L if abs(total_bend_deg) > _EPS else 0.0

    pts = np.empty((n, 3))
    for i in range(n):
        t = rise * i
        if kappa == 0.0:
            A = np.array([0.0, 0.0, t])
            N = np.array([1.0, 0.0, 0.0])
        else:
            R = 1.0 / kappa
            A = np.array([R * (1 - np.cos(kappa * t)), 0.0, R * np.sin(kappa * t)])
            N = np.array([np.cos(kappa * t), 0.0, -np.sin(kappa * t)])
        B = np.array([0.0, 1.0, 0.0])       # arc lies in x-z, so binormal is +y
        th = np.radians(turn * i + phase)
        pts[i] = A + r * (np.cos(th) * N + np.sin(th) * B)
    return pts


def true_half_chord_bend(total_bend_deg):
    """Axis turning a *half-chord* construction sees, given a full-span bend.

    For a circular arc the chord of the first half points along the tangent at its
    midpoint, likewise the second half, so the angle between the two chords is
    exactly half the full-span tangent turning. Any chord-of-halves metric
    therefore has an expected slope of 0.5 against `total_bend_deg`, and that is
    geometry, not a defect.
    """
    return 0.5 * total_bend_deg


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _angle(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu < _EPS or nv < _EPS:
        return float("nan")
    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def _principal_axis(points):
    """Unit first principal axis, oriented first->last. Also returns s1/s0."""
    P = np.asarray(points, float)
    C = P - P.mean(axis=0)
    _, s, Vt = np.linalg.svd(C, full_matrices=False)
    if s[0] < _EPS:
        return None, 1.0
    ax = Vt[0]
    if np.dot(ax, P[-1] - P[0]) < 0:
        ax = -ax
    cond = float(s[1] / s[0]) if s[0] > _EPS else 1.0
    return ax / np.linalg.norm(ax), cond


def smooth_trace(ca, w):
    """Sliding-window mean of consecutive Ca.

    For a regular helix of twist theta, averaging w consecutive points shrinks the
    radial wobble by the Dirichlet factor |sin(w*theta/2) / (w*sin(theta/2))|, so
    the smoothed points sit on (or very near) the helical axis. w=4 gives 0.11 for
    an alpha-helix (2.30 A -> 0.26 A residual) and exactly 0 for a beta-strand.
    """
    P = np.asarray(ca, float)
    if w <= 1:
        return P
    k = P.shape[0] - w + 1
    if k < 1:
        raise ValueError(f"need >= {w} points to smooth with w={w}, got {P.shape[0]}")
    return np.array([P[i:i + w].mean(axis=0) for i in range(k)])


def dirichlet_shrink(w, turn_deg):
    """Analytic radial shrink factor of `smooth_trace` for a regular helix."""
    th = np.radians(turn_deg)
    num = np.sin(w * th / 2.0)
    den = w * np.sin(th / 2.0)
    return float("nan") if abs(den) < _EPS else abs(num / den)


# --------------------------------------------------------------------------- #
# candidates
# --------------------------------------------------------------------------- #
class LegacyPCA5:
    """v1 metric, unchanged: 5 Ca, PCA axis of each overlapping half, angle between.

    Kept as the baseline row so the gate compares against it rather than assuming
    it is wrong. It does NOT claim to be an axis metric -- that is the finding.
    """
    name = "v1_pca5"
    span = 5
    is_axis_metric = False
    smooth = 1

    def __call__(self, ca):
        P = np.asarray(ca, float)
        if P.shape != (5, 3):
            raise ValueError(f"expected (5,3), got {P.shape}")
        d1, _ = _principal_axis(P[0:3])
        d2, _ = _principal_axis(P[2:5])
        if d1 is None or d2 is None:
            return float("nan")
        return _angle(d1, d2)

    def conditioning(self, ca):
        P = np.asarray(ca, float)
        return max(_principal_axis(P[0:3])[1], _principal_axis(P[2:5])[1])


class SmoothedChord:
    """Smooth the Ca trace to its axis, then the angle between half-chords.

    No PCA anywhere, so none of the s1/s0 ill-conditioning of AUDIT.md F11 can
    arise. Expected slope 0.5 vs full-span bend (see `true_half_chord_bend`).
    """
    is_axis_metric = True

    def __init__(self, span=9, smooth=4):
        self.span, self.smooth = span, smooth
        self.name = f"smooth_chord_s{span}_w{smooth}"
        if span - smooth + 1 < 3:
            raise ValueError("need >= 3 smoothed points")

    def __call__(self, ca):
        P = np.asarray(ca, float)
        if P.shape != (self.span, 3):
            raise ValueError(f"expected ({self.span},3), got {P.shape}")
        s = smooth_trace(P, self.smooth)
        mid = (len(s) - 1) // 2
        return _angle(s[mid] - s[0], s[-1] - s[mid])

    def conditioning(self, ca):
        return 0.0          # no eigen-decomposition involved


class SmoothedPCA:
    """Smooth first, THEN the v1 PCA-of-halves construction.

    Included to separate the two possible causes of v1's behaviour: is the problem
    the PCA fit, or the absence of smoothing? If this candidate behaves like
    SmoothedChord, smoothing was the whole story and PCA was incidental.
    """
    is_axis_metric = True

    def __init__(self, span=9, smooth=4):
        self.span, self.smooth = span, smooth
        self.name = f"smooth_pca_s{span}_w{smooth}"

    def __call__(self, ca):
        P = np.asarray(ca, float)
        s = smooth_trace(P, self.smooth)
        mid = (len(s) - 1) // 2
        d1, _ = _principal_axis(s[:mid + 1])
        d2, _ = _principal_axis(s[mid:])
        if d1 is None or d2 is None:
            return float("nan")
        return _angle(d1, d2)

    def conditioning(self, ca):
        s = smooth_trace(np.asarray(ca, float), self.smooth)
        mid = (len(s) - 1) // 2
        return max(_principal_axis(s[:mid + 1])[1], _principal_axis(s[mid:])[1])


class CircleFitCurvature:
    """Fit a circle to the smoothed axis points; report turning over the span.

    Uses all smoothed points rather than two chords, so it averages more of the
    coordinate noise, at the cost of assuming the axis is an arc (fine locally).
    Reports arclen * kappa in degrees, comparable to the full-span bend, so its
    expected slope is 1.0, not 0.5.
    """
    is_axis_metric = True

    def __init__(self, span=9, smooth=4):
        self.span, self.smooth = span, smooth
        self.name = f"circlefit_s{span}_w{smooth}"

    def __call__(self, ca):
        s = smooth_trace(np.asarray(ca, float), self.smooth)
        if len(s) < 3:
            return float("nan")
        # project onto the best-fit plane of the smoothed points, fit a circle there
        c = s.mean(axis=0)
        _, _, Vt = np.linalg.svd(s - c, full_matrices=False)
        e1, e2 = Vt[0], Vt[1]
        x = (s - c) @ e1
        y = (s - c) @ e2
        # algebraic circle fit: x^2+y^2 + D x + E y + F = 0
        A = np.column_stack([x, y, np.ones_like(x)])
        b = -(x ** 2 + y ** 2)
        try:
            D, E, F = np.linalg.lstsq(A, b, rcond=None)[0]
        except np.linalg.LinAlgError:
            return float("nan")
        cx, cy = -D / 2.0, -E / 2.0
        r2 = cx ** 2 + cy ** 2 - F
        if r2 <= _EPS:
            return 0.0
        R = float(np.sqrt(r2))
        arclen = float(np.linalg.norm(np.diff(s, axis=0), axis=1).sum())
        return float(np.degrees(arclen / R))

    def conditioning(self, ca):
        return 0.0


class BisectorAxis:
    """Local axis from Ca bisectors (Kahn / HELANAL construction).

    For a regular helix the bisector at residue i,
        b_i = unit(Ca[i-1]-Ca[i]) + unit(Ca[i+1]-Ca[i]),
    points radially inward toward the axis. Two consecutive bisectors are two
    radial vectors at different heights, so cross(b_i, b_i+1) lies ALONG the axis
    -- exactly, for any twist. That is strictly better than sliding-window
    averaging, which only cancels the coil for twists where the Dirichlet factor
    happens to vanish.

    The known failure mode is deliberate and reported, not hidden: at twist 180
    deg (beta-strand) consecutive bisectors are antiparallel and the cross product
    vanishes. `degeneracy` returns |sin| between consecutive bisectors, so the
    gate can show where this estimator is and is not admissible.
    """
    is_axis_metric = True
    smooth = 1

    def __init__(self, span=9, min_sin=0.15):
        self.span, self.min_sin = span, min_sin
        self.name = f"bisector_s{span}"

    def _dirs(self, P):
        b = []
        for i in range(1, len(P) - 1):
            u, v = P[i - 1] - P[i], P[i + 1] - P[i]
            nu, nv = np.linalg.norm(u), np.linalg.norm(v)
            if nu < _EPS or nv < _EPS:
                return None
            bi = u / nu + v / nv
            nb = np.linalg.norm(bi)
            if nb < _EPS:
                return None
            b.append(bi / nb)
        d = []
        for k in range(len(b) - 1):
            c = np.cross(b[k], b[k + 1])
            if np.linalg.norm(c) < self.min_sin:
                return None                      # bisectors (anti)parallel
            c = c / np.linalg.norm(c)
            if np.dot(c, P[k + 2] - P[k + 1]) < 0:
                c = -c
            d.append(c)
        return d

    def __call__(self, ca):
        d = self._dirs(np.asarray(ca, float))
        return float("nan") if d is None or len(d) < 2 else _angle(d[0], d[-1])

    def degeneracy(self, ca):
        """Smallest |sin| between consecutive bisectors; ~0 means inadmissible."""
        P = np.asarray(ca, float)
        b = []
        for i in range(1, len(P) - 1):
            u, v = P[i - 1] - P[i], P[i + 1] - P[i]
            bi = u / np.linalg.norm(u) + v / np.linalg.norm(v)
            nb = np.linalg.norm(bi)
            if nb < _EPS:
                return 0.0
            b.append(bi / nb)
        return float(min(np.linalg.norm(np.cross(b[k], b[k + 1]))
                         for k in range(len(b) - 1)))

    def conditioning(self, ca):
        return 0.0


# Calibration constants: raw response per degree of true full-span axis bend.
# Measured against `helix_on_arc` ground truth; see measurement_gate.py stage 2.
# The bisector constant is essentially universal over helical twists
# (0.6248 / 0.6256 / 0.6307 / 0.6310 for alpha / pi / 3-10 / PPII, all R^2 = 1.0000),
# which is why one number suffices for that whole branch.
CAL_BISECTOR = 0.6280
CAL_SMOOTH_W4_EXTENDED = 0.3125          # exact for twist 180 deg, R^2 = 1.0000


class HybridAxisBend:
    """Calibrated axis-bend estimate, in degrees of TRUE full-span axis turning.

    Dispatches on the actual degeneracy condition rather than on an external
    secondary-structure annotation -- so it is not circular:

      * bisectors well separated (helical twist)  -> BisectorAxis / CAL_BISECTOR
      * bisectors (anti)parallel (extended twist) -> SmoothedChord(w=4) / CAL_...

    The smoothed branch has a nonzero straight-axis offset away from an exact
    180-degree strand. Switching directly at the bisector degeneracy threshold
    creates a discontinuity, so the uncertain handover band abstains.

    `branch(ca)` reports which estimator fired, which belongs in the v2 output
    schema as a QC field.
    """
    is_axis_metric = True
    smooth = 4

    def __init__(self, span=9, min_sin=0.15, smooth_max_sin=0.035):
        if not 0 <= smooth_max_sin < min_sin:
            raise ValueError("need 0 <= smooth_max_sin < min_sin")
        self.span, self.min_sin = span, min_sin
        self.smooth_max_sin = smooth_max_sin
        self.name = f"hybrid_s{span}"
        self._bis = BisectorAxis(span=span, min_sin=min_sin)
        self._sm = SmoothedChord(span=span, smooth=4)

    def branch(self, ca):
        degeneracy = self._bis.degeneracy(ca)
        if degeneracy >= self.min_sin:
            return "bisector"
        if degeneracy <= self.smooth_max_sin:
            return "smoothed"
        return "abstain"

    def __call__(self, ca):
        P = np.asarray(ca, float)
        branch = self.branch(P)
        if branch == "bisector":
            v = self._bis(P)
            return float("nan") if v != v else v / CAL_BISECTOR
        if branch == "smoothed":
            return self._sm(P) / CAL_SMOOTH_W4_EXTENDED
        return float("nan")

    def conditioning(self, ca):
        return 0.0


def default_candidates():
    """The shortlist the gate scores. Deliberately small and interpretable."""
    return [
        LegacyPCA5(),
        SmoothedChord(span=9, smooth=4),
        SmoothedPCA(span=9, smooth=4),
        CircleFitCurvature(span=9, smooth=4),
        BisectorAxis(span=9),
        HybridAxisBend(span=9),
    ]


# --------------------------------------------------------------------------- #
# per-structure coordinate uncertainty
# --------------------------------------------------------------------------- #
def sigma_from_bfactor(B):
    """Per-axis positional spread implied by an isotropic B-factor, in Angstrom.

    B = 8 pi^2 <u^2> with <u^2> the mean-square displacement along one axis, so
    sigma_1d = sqrt(B / (8 pi^2)).

    NOTE this is the atom's positional *spread* (thermal + static disorder), which
    is the right quantity for "how far apart will this Ca sit in two independent
    crystals". It is an upper bound on refinement *coordinate error*; the
    Cruickshank DPI is the lower-bound counterpart and needs reflection counts,
    which are not in a coordinate file. The gate sweeps sigma explicitly rather
    than committing to either.
    """
    return float(np.sqrt(np.asarray(B, float) / (8.0 * np.pi ** 2)))


def noise_gain(metric, ca0, sigma_xyz, n=4000, rng=None):
    """SD of `metric` (deg) under iid Gaussian Ca perturbation of `sigma_xyz` (A)."""
    rng = rng or np.random.default_rng(0)
    P = np.asarray(ca0, float)
    v = np.array([metric(P + rng.normal(0, sigma_xyz, P.shape)) for _ in range(n)])
    v = v[np.isfinite(v)]
    return float(v.std()), float(v.mean())
