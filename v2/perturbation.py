"""
Controlled perturbation of REAL coordinates, and the local Jacobian of a metric.

Why this exists
---------------
`measurement_gate.py` measures every metric's response to axis bend on IDEAL
geometry (`helix_on_arc`). That established v1's metric has slope 0.013 on an
ideal alpha-helix. It does NOT establish the same for real helices, which are
irregular: bent, frayed, variable in rise and twist, and noisy.

That transfer has to be measured, not assumed. This module does it by imposing a
KNOWN deformation on an ARBITRARY window -- real or ideal -- and reading off
d(metric)/d(theta). Run it on real T4L windows and the blindness claim either
transfers or it does not.

The deformation
---------------
`rigid_half_bend(P, theta, azimuth)` rotates the C-terminal half of the window
rigidly about an axis through the midpoint. Properties, all of which matter:

  * every Ca-Ca distance is preserved exactly (each half is rigid, and the
    junction atom lies on the rotation axis), so it is a pure bending
    deformation -- it does not stretch bonds or smuggle in twist within a half;
  * the two halves' axes turn by exactly `theta` relative to one another, so
    `theta` is ground truth by construction, with no reference to any idealised
    model of the starting conformation;
  * the bend direction is a free parameter (`azimuth`), so the Jacobian is
    averaged over it rather than being an artifact of one arbitrary plane.

What to read off
----------------
The absolute Jacobian feeds a real-data ceiling once within-variant noise is
known. But the robust quantity is the RATIO of two metrics' Jacobians on the
same window under the same perturbation: it is independent of how one chooses to
define "true bend", so it transfers the v1-vs-v2 comparison to real coordinates
without depending on a convention.
"""
from __future__ import annotations

import numpy as np

from metrics import IDEAL, helix_on_arc, smooth_trace

_EPS = 1e-12
DEFAULT_AZIMUTHS = (0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0)


def _rot(axis, ang_rad):
    """Rodrigues rotation matrix about a unit axis."""
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    K = np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])
    return np.eye(3) + np.sin(ang_rad) * K + (1 - np.cos(ang_rad)) * (K @ K)


def local_frame(P, mid=None, smooth=4):
    """(tangent, u, v) at the window midpoint.

    The tangent comes from the SMOOTHED trace, so it is the local axis direction
    rather than a coil chord -- otherwise, on a helix, "perpendicular to the
    chain" would mean perpendicular to a coil segment and the perturbation would
    mix bend with twist. u, v span the plane perpendicular to it.
    """
    P = np.asarray(P, float)
    n = len(P)
    mid = n // 2 if mid is None else mid
    s = smooth_trace(P, smooth) if n >= smooth + 2 else P
    t = s[-1] - s[0]
    nt = np.linalg.norm(t)
    if nt < _EPS:
        t = P[-1] - P[0]
        nt = np.linalg.norm(t)
    t = t / nt
    seed = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(seed, t)) > 0.9:
        seed = np.array([1.0, 0.0, 0.0])
    u = np.cross(t, seed)
    u /= np.linalg.norm(u)
    v = np.cross(t, u)
    return t, u, v


def rigid_half_bend(P, theta_deg, azimuth_deg=0.0, mid=None, smooth=4):
    """Rotate the C-terminal half rigidly about the midpoint by `theta_deg`.

    Preserves every Ca-Ca distance. `azimuth_deg` selects the bend direction in
    the plane perpendicular to the local axis.
    """
    P = np.asarray(P, float).copy()
    n = len(P)
    mid = n // 2 if mid is None else mid
    t, u, v = local_frame(P, mid, smooth)
    az = np.radians(azimuth_deg)
    axis = np.cos(az) * u + np.sin(az) * v        # perpendicular to the axis
    R = _rot(axis, np.radians(theta_deg))
    pivot = P[mid]
    P[mid + 1:] = (P[mid + 1:] - pivot) @ R.T + pivot
    return P


def check_bond_lengths(P0, P1, tol=1e-8):
    """True if the deformation preserved all consecutive Ca-Ca distances."""
    d0 = np.linalg.norm(np.diff(np.asarray(P0, float), axis=0), axis=1)
    d1 = np.linalg.norm(np.diff(np.asarray(P1, float), axis=0), axis=1)
    return bool(np.abs(d0 - d1).max() < tol)


def jacobian(metric, P, thetas=(0.0, 1.0, 2.0, 3.0, 5.0),
             azimuths=DEFAULT_AZIMUTHS, signed=False):
    """d(metric)/d(theta) at `P`, averaged over bend azimuth.

    Returns dict(slope, r2, intercept, n_ok). `slope` is in degrees of metric
    response per degree of imposed relative half-rotation.

    With `signed=False` the response is measured as |m(theta) - m(0)|, which is
    what matters for a magnitude-based mover score and avoids cancellation when
    the metric's sign convention differs from the imposed bend direction.
    """
    P = np.asarray(P, float)
    base = metric(P)
    if not np.isfinite(base):
        return dict(slope=float("nan"), r2=float("nan"),
                    intercept=float("nan"), n_ok=0)
    resp, xs = [], []
    for th in thetas:
        vals = []
        for az in azimuths:
            m = metric(rigid_half_bend(P, th, az))
            if np.isfinite(m):
                vals.append(m - base if signed else abs(m - base))
        if vals:
            resp.append(float(np.mean(vals)))
            xs.append(th)
    if len(xs) < 3:
        return dict(slope=float("nan"), r2=float("nan"),
                    intercept=float("nan"), n_ok=len(xs))
    A = np.column_stack([xs, np.ones(len(xs))])
    slope, icpt = np.linalg.lstsq(A, np.array(resp), rcond=None)[0]
    pred = A @ [slope, icpt]
    y = np.array(resp)
    sst = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ((y - pred) ** 2).sum() / sst if sst > 1e-12 else float("nan")
    return dict(slope=float(slope), r2=float(r2), intercept=float(icpt),
                n_ok=len(xs))


def repeatability(metric, P, sigma_xyz, n=2000, rng=None):
    """SD of `metric` at `P` under iid Ca perturbation -- the iid noise term only.

    NOTE this is the IID component. Real crystal-to-crystal variation also has
    systematic components (packing, space group, refinement protocol and era,
    construct/background, cryo vs RT, ligand state) that are correlated across
    crystals of one variant and therefore do NOT shrink as 1/sqrt(n). Estimating
    those needs real redundant crystals; see `seed_test.py`.
    """
    rng = rng or np.random.default_rng(0)
    P = np.asarray(P, float)
    v = np.array([metric(P + rng.normal(0, sigma_xyz, P.shape)) for _ in range(n)])
    v = v[np.isfinite(v)]
    return float(v.std()) if len(v) > 50 else float("nan")


# --------------------------------------------------------------------------- #
# self-validation: the machinery must reproduce known behaviour on ideal
# geometry BEFORE it is trusted on real coordinates.
# --------------------------------------------------------------------------- #
def _self_test():
    from metrics import HybridAxisBend, LegacyPCA5, SmoothedChord

    print("=" * 78)
    print("  perturbation.py self-test -- validate on IDEAL geometry first")
    print("=" * 78)

    print("\n  1. does rigid_half_bend preserve every Ca-Ca distance?")
    ok = True
    for ss, (r, rise, turn) in IDEAL.items():
        P0 = helix_on_arc(9, r, rise, turn, 0.0)
        for th in (1.0, 5.0, 20.0):
            for az in (0.0, 90.0, 210.0):
                good = check_bond_lengths(P0, rigid_half_bend(P0, th, az))
                ok &= good
        print(f"      {ss:12s} {'preserved' if ok else 'BROKEN'}")
    print(f"      -> {'PASS' if ok else 'FAIL'} (pure bending, no bond stretching)")

    print("\n  2. is the imposed theta recovered as relative half-axis turning?")
    P0 = helix_on_arc(9, *IDEAL['alpha_helix'], 0.0)
    for th in (0.0, 2.0, 5.0, 10.0, 20.0):
        Pb = rigid_half_bend(P0, th, 0.0)
        s0, s1 = smooth_trace(P0, 4), smooth_trace(Pb, 4)
        a0 = s0[-1] - s0[0]
        a1 = s1[-1] - s1[0]

        def ang(x, y):
            return np.degrees(np.arccos(np.clip(
                np.dot(x, y) / np.linalg.norm(x) / np.linalg.norm(y), -1, 1)))
        print(f"      imposed {th:5.1f} deg -> smoothed end-to-end axis moved "
              f"{ang(a0, a1):6.2f} deg")

    print("\n  3. Jacobian on ideal geometry: v1 vs v2 (the claim under test)")
    v1, hyb, sm = LegacyPCA5(), HybridAxisBend(9), SmoothedChord(9, 4)
    print(f"      {'conformation':>14} {'v1 slope':>10} {'hybrid':>10} "
          f"{'smooth':>10} {'hybrid/v1':>11}")
    for ss in ("alpha_helix", "beta_strand", "three_ten", "pp_ii"):
        r, rise, turn = IDEAL[ss]
        P9 = helix_on_arc(9, r, rise, turn, 0.0)
        P5 = helix_on_arc(5, r, rise, turn, 0.0)
        j1 = jacobian(v1, P5)["slope"]           # v1 needs a 5-Ca window
        jh = jacobian(hyb, P9)["slope"]
        js = jacobian(sm, P9)["slope"]
        ratio = jh / j1 if j1 and np.isfinite(j1) and abs(j1) > 1e-9 else float("nan")
        print(f"      {ss:>14} {j1:>10.4f} {jh:>10.4f} {js:>10.4f} {ratio:>11.1f}x")
    print("\n      These are the IDEAL-geometry reference values. Running the same")
    print("      function on real T4L windows is what tests whether v1's helical")
    print("      blindness is a property of ideal helices or of real ones.")
    print("\n  NOTE the v1 Jacobian here uses a 5-Ca window and the others a 9-Ca")
    print("  window, because that is each metric's own span. The perturbation is")
    print("  therefore not literally identical between them, so read the ratio as")
    print("  indicative; the like-for-like comparison is v1 vs v2 on the SAME real")
    print("  window, both spans centred on the same mutated residue.")


if __name__ == "__main__":
    _self_test()
