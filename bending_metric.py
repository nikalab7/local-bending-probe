"""
proteinX -- canonical local backbone bending metric.

THIS IS THE SINGLE SOURCE OF TRUTH for "bending". The same function MUST be
used for all three consumers:
    (1) training labels                (experimental C-alpha)
    (2) ground-truth WT/mutant deltas  (experimental C-alpha)
    (3) AlphaFold baseline scoring     (predicted  C-alpha)
If these three ever compute bending differently, the "we beat AF2" claim is
meaningless. Import this; do not reimplement.

Definition (intrinsic, superposition-free)
-------------------------------------------
Outcome window = 5 consecutive C-alpha atoms for residues i .. i+4 (matches
the v2 outcome region; entry residues i-2, i-1 are NEVER touched here, so the
label cannot leak entry geometry).

Split at the midpoint i+2 into two overlapping halves:
    H1 = {i,   i+1, i+2}
    H2 = {i+2, i+3, i+4}
Fit the best-fit line to each half (first principal axis via SVD; all three
points contribute, so it is more robust than an endpoint vector), orient each
axis along the chain, and take the angle between them, in degrees, in [0,180].

Because the value is built only from C-alpha *internal* geometry it is
invariant to global rotation/translation. Therefore delta-bending between two
structures needs NO superposition -- removing a whole class of
alignment-induced artifact. (The "no indels" rule in pair-mining guarantees
window i..i+4 maps to the same residues in WT and mutant.)

Cross-check metric (optional): unsigned sum of the three C-alpha pseudo-bond
exterior angles. Should track this within a few degrees; large disagreement
flags an S-shaped (sign-cancelling) window worth inspecting.
"""

from __future__ import annotations
import numpy as np

_EPS = 1e-9


def _principal_axis(points: np.ndarray) -> np.ndarray:
    """Unit first-principal-axis of a small point set, oriented first->last."""
    P = np.asarray(points, dtype=float)
    centered = P - P.mean(axis=0)
    # SVD: rows of Vt are the principal axes, largest singular value first.
    _, s, Vt = np.linalg.svd(centered, full_matrices=False)
    if s[0] < _EPS:                       # all points coincide
        raise ValueError("degenerate half: coincident C-alpha coordinates")
    axis = Vt[0]
    if np.dot(axis, P[-1] - P[0]) < 0:    # orient along the chain (first->last)
        axis = -axis
    return axis / np.linalg.norm(axis)


def bending_angle(ca_window: np.ndarray) -> float:
    """
    Bending angle in degrees for one outcome window.

    ca_window : (5, 3) C-alpha coordinates for residues i .. i+4, in chain
                order, no NaNs.
    """
    P = np.asarray(ca_window, dtype=float)
    if P.shape != (5, 3):
        raise ValueError(f"expected (5,3) C-alpha window, got {P.shape}")
    if not np.isfinite(P).all():
        raise ValueError("window contains missing (non-finite) coordinates")
    d1 = _principal_axis(P[0:3])          # H1 = i, i+1, i+2
    d2 = _principal_axis(P[2:5])          # H2 = i+2, i+3, i+4
    cos = float(np.clip(np.dot(d1, d2), -1.0, 1.0))
    return float(np.degrees(np.arccos(cos)))


def delta_bending(ca_window_wt: np.ndarray, ca_window_mut: np.ndarray) -> float:
    """Signed bending change for a matched WT/mutant window (mut - wt).

    Positive => the mutation adds curvature. No superposition required.
    """
    return bending_angle(ca_window_mut) - bending_angle(ca_window_wt)


# --- chain-continuity guard: run before trusting any window ----------------
_CA_CA_IDEAL = 3.8     # Angstrom, trans peptide
_CA_CA_BREAK = 4.5     # consecutive C-alpha gap above this => chain break


def is_continuous(ca_window: np.ndarray) -> bool:
    """False if a gap/chain-break makes the window's geometry meaningless."""
    P = np.asarray(ca_window, dtype=float)
    if P.shape != (5, 3) or not np.isfinite(P).all():
        return False
    d = np.linalg.norm(np.diff(P, axis=0), axis=1)
    return bool((d <= _CA_CA_BREAK).all())


if __name__ == "__main__":
    # Self-test. A window is built so its TRUE bend equals a chosen theta:
    # entry half runs into the apex along +x, exit half leaves along theta.
    def _bent(theta_deg, scale=_CA_CA_IDEAL):
        t = np.radians(theta_deg)
        d1 = np.array([1.0, 0.0, 0.0])
        d2 = np.array([np.cos(t), np.sin(t), 0.0])
        return np.vstack([-2 * d1, -d1, np.zeros(3), d2, 2 * d2]) * scale

    # Span acute -> right -> OBTUSE -> near-reversal so a PCA sign/orientation
    # bug (which would read a ~170 deg hairpin as ~10 deg) cannot slip through.
    ok = True
    for th in (0.0, 30.0, 90.0, 135.0, 160.0, 170.0):
        got = bending_angle(_bent(th))
        flag = "OK" if abs(got - th) < 0.5 else "FAIL"
        ok &= flag == "OK"
        print(f"  true {th:6.1f} deg  ->  metric {got:6.2f} deg   [{flag}]")

    # Hard orientation guard: a chain reversal MUST read obtuse, never acute.
    assert bending_angle(_bent(170.0)) > 150.0, "SIGN BUG: reversal read as acute"

    broken = np.array([[0, 0, 0], [3.8, 0, 0], [7.6, 0, 0],
                       [20, 0, 0], [24, 0, 0]], float)
    print(f"continuity: straight={is_continuous(_bent(0.0))}  broken={is_continuous(broken)}")
    print("ALL PASS" if ok else "SELF-TEST FAILED")
