"""
proteinX -- cavity follow-up to THEORY.md T2 (CAVITY.md).

Exploratory; steps, measurements and decision rules were committed in
CAVITY.md before this was run.

  python cavity.py   -> results/cavity.json
"""
from __future__ import annotations
import json
import os
import numpy as np

import delta_model as dm
import final_eval as fe
from mechanisms import path_of, chains_of, atoms, at, kabsch, window_geometry
from theory import radial_one, FAR, MIN_FAR, MIN_ATOMS

SC = dict(A="CB", R="CB CG CD NE CZ NH1 NH2", N="CB CG OD1 ND2", D="CB CG OD1 OD2", C="CB SG",
          Q="CB CG CD OE1 NE2", E="CB CG CD OE1 OE2", G="", H="CB CG ND1 CD2 CE1 NE2",
          I="CB CG1 CG2 CD1", L="CB CG CD1 CD2", K="CB CG CD CE NZ", M="CB CG SD CE",
          F="CB CG CD1 CD2 CE1 CE2 CZ", P="CB CG CD", S="CB OG", T="CB OG1 CG2",
          W="CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2", Y="CB CG CD1 CD2 CE1 CE2 CZ OH", V="CB CG1 CG2")
SC = {k: set(v.split()) for k, v in SC.items()}
CONTACT_B, N_NULL, SEED, N_BOOT = 6.0, 1500, 0, 2000


def inward_one(wt, mt, r, removed):
    """Measurement B: mean displacement of cavity-lining atoms toward the removed-atom centroid."""
    rem = [wt[r][1][n] for n in removed if r in wt and n in wt[r][1]]
    ca_r = at(wt, r, "CA")
    if not rem or ca_r is None:
        return None
    rem = np.array(rem); c = rem.mean(0)
    P, Q = [], []
    for k, (aa, a) in wt.items():
        if "CA" in a and k in mt and "CA" in mt[k][1] and np.linalg.norm(a["CA"] - ca_r) > FAR:
            P.append(mt[k][1]["CA"]); Q.append(a["CA"])
    if len(P) < MIN_FAR:
        return None
    R, t = kabsch(np.array(P), np.array(Q))
    vals = []
    for k, (aa, a) in wt.items():
        if abs(k - r) < 2 or k not in mt:
            continue
        for name, xw in a.items():
            if np.min(np.linalg.norm(rem - xw, axis=1)) > CONTACT_B:
                continue
            xm = mt[k][1].get(name)
            if xm is None or np.linalg.norm(c - xw) < 1e-6:
                continue
            v = (c - xw) / np.linalg.norm(c - xw)
            vals.append(float(np.dot(xm @ R + t - xw, v)))
    return float(np.mean(vals)) if len(vals) >= MIN_ATOMS else None


def measure(x, expect_aa, want_b):
    """(A, B) medians over the row's mutant (or held-out WT) crystals."""
    wt = atoms(x["scaffold_path"], x["scaffold_chain"])
    removed = SC[x["wt"]] - SC[x["mut"]]
    A, B = [], []
    for pdb in x["pdbs"]:
        path = path_of(pdb)
        if path is None:
            continue
        best, rm = None, np.inf
        for ch in chains_of(path):
            mt = atoms(path, ch)
            if mt.get(x["r"], ("?",))[0] != expect_aa:
                continue
            g = window_geometry(wt, mt, x["s"])
            if g and g["ca_rmsd"] < rm:
                best, rm = mt, g["ca_rmsd"]
        if best is None:
            continue
        a = radial_one(wt, best, x["r"])
        if a is not None:
            A.append(a)
        if want_b:
            b = inward_one(wt, best, x["r"], removed)
            if b is not None:
                B.append(b)
    med = lambda v: float(np.median(v)) if v else np.nan
    return med(A), med(B)


def _boot(cl, seed=SEED):
    u = np.unique(cl); by = {c: np.flatnonzero(cl == c) for c in u}
    rng = np.random.default_rng(seed)
    for _ in range(N_BOOT):
        yield np.concatenate([by[c] for c in rng.choice(u, len(u))])


def q95(v):
    v = [x for x in v if np.isfinite(x)]
    return [float(x) for x in np.percentile(v, [2.5, 97.5])] if v else None


def summary(v, cl, real, null):
    """Mean of real rows, mean of null rows, and real - null, cluster bootstrap."""
    real, null = real & np.isfinite(v), null & np.isfinite(v)
    bs_m, bs_d = [], []
    for i in _boot(cl):
        r_, n_ = real[i], null[i]
        if r_.any():
            bs_m.append(v[i][r_].mean())
            if n_.any():
                bs_d.append(v[i][r_].mean() - v[i][n_].mean())
    return dict(n=int(real.sum()), n_null=int(null.sum()), clusters=int(len(set(cl[real]))),
                mean=float(v[real].mean()) if real.any() else np.nan, ci=q95(bs_m),
                null_mean=float(v[null].mean()) if null.any() else np.nan,
                diff=float(v[real].mean() - v[null].mean()) if real.any() and null.any() else np.nan,
                diff_ci=q95(bs_d))


def neg(d):
    return d["ci"] is not None and d["ci"][1] < 0 and d["diff_ci"] is not None and d["diff_ci"][1] < 0


def pos(d):
    return d["ci"] is not None and d["ci"][0] > 0 and d["diff_ci"] is not None and d["diff_ci"][0] > 0


def asymmetry(vA, cl, ov, cav, null):
    def idx(i):
        o, c, n = ov[i] & np.isfinite(vA[i]), cav[i] & np.isfinite(vA[i]), null[i] & np.isfinite(vA[i])
        if not (o.any() and c.any() and n.any()):
            return np.nan
        m = vA[i][n].mean()
        return (vA[i][o].mean() - m) + (vA[i][c].mean() - m)
    allidx = np.arange(len(vA))
    return dict(index=float(idx(allidx)), ci=q95([idx(i) for i in _boot(cl)]),
                overpacking_minus_null=float(vA[ov & np.isfinite(vA)].mean() - vA[null & np.isfinite(vA)].mean()),
                cavity_minus_null=float(vA[cav & np.isfinite(vA)].mean() - vA[null & np.isfinite(vA)].mean()),
                n_overpacking=int((ov & np.isfinite(vA)).sum()), n_cavity=int((cav & np.isfinite(vA)).sum()),
                n_null=int((null & np.isfinite(vA)).sum()))


def fmt(d):
    return (f"mean {d['mean']:+.3f} {d['ci']} | null {d['null_mean']:+.3f} | diff {d['diff']:+.3f} {d['diff_ci']} "
            f"(n {d['n']}, null {d['n_null']}, clusters {d['clusters']})")


def main():
    rd, nd, Fd = fe.load("dev"); rl, nl, Fl = fe.load("lockbox")
    rows, F, null = rd + rl, Fd + Fl, nd + nl
    hse = dm.matrix(F, ["hse_up"])[:, 0]; med = np.nanmedian(hse)
    dv = np.array([dm.VOL[x["mut"]] - dm.VOL[x["wt"]] for x in rows], float)
    gp = np.array([x["wt"] in "GP" or x["mut"] in "GP" for x in rows])
    t4 = np.array([x["protein"] == "T4L" for x in rows])
    buried = hse >= med
    cav = buried & (dv <= -25) & ~gp
    ov = buried & (dv >= 25) & ~gp
    R = dict(note="exploratory; pre-specified in CAVITY.md")
    # ---- T4L: cavity rows + all T4L null rows; clusters = residue site
    i_r = np.flatnonzero(cav & t4)
    nl_t4 = [x for x in null if x["protein"] == "T4L"]
    mr = [measure(rows[i], rows[i]["mut"], True) for i in i_r]
    mn = [measure(x, x["wt"], True) for x in nl_t4]
    vA = np.array([m[0] for m in mr] + [m[0] for m in mn]); vB = np.array([m[1] for m in mr] + [m[1] for m in mn])
    cl = np.array([str(rows[i]["r"]) for i in i_r] + [str(x["r"]) for x in nl_t4])
    real = np.r_[np.ones(len(i_r), bool), np.zeros(len(nl_t4), bool)]
    s1 = summary(vA, cl, real, ~real)
    R["step1_T4L_measurementA"] = s1
    R["step1_collapse"] = bool(neg(s1))
    print(f"Step 1 T4L, measurement A (negative = collapse): {fmt(s1)} -> collapse: {R['step1_collapse']}")
    # ---- non-T4L measurement A (needed for the asymmetry index in either branch)
    nt_null = [x for x in null if x["protein"] != "T4L"]
    pick = np.random.default_rng(SEED).choice(len(nt_null), min(N_NULL, len(nt_null)), replace=False)
    nt_null = [nt_null[i] for i in pick]

    def non_t4l(want_b):
        idx = np.flatnonzero((cav | ov) & ~t4)
        m_r = [measure(rows[i], rows[i]["mut"], want_b) for i in idx]
        m_n = [measure(x, x["wt"], want_b) for x in nt_null]
        a = np.array([m[0] for m in m_r] + [m[0] for m in m_n])
        b = np.array([m[1] for m in m_r] + [m[1] for m in m_n])
        fam = np.array([str(rows[i]["family"]) for i in idx] + [str(x["family"]) for x in nt_null])
        isr = np.r_[np.ones(len(idx), bool), np.zeros(len(m_n), bool)]
        c_m = np.r_[cav[idx], np.zeros(len(m_n), bool)]; o_m = np.r_[ov[idx], np.zeros(len(m_n), bool)]
        return a, b, fam, isr, c_m, o_m
    if R["step1_collapse"]:
        a, _, fam, isr, c_m, o_m = non_t4l(False)
        s3 = summary(a, fam, c_m, ~isr)
        R["step3_nonT4L_measurementA"] = s3
        R["asymmetry_nonT4L"] = asymmetry(a, fam, o_m, c_m, ~isr)
        print(f"Step 3 non-T4L, measurement A: {fmt(s3)}")
    else:
        s2 = summary(vB, cl, real, ~real)
        R["step2_T4L_measurementB"] = s2
        R["B_validated"] = bool(pos(s2))
        print(f"Step 2 T4L, measurement B (positive = into cavity): {fmt(s2)} -> validated: {R['B_validated']}")
        if R["B_validated"]:
            a, b, fam, isr, c_m, o_m = non_t4l(True)
            s2b = summary(b, fam, c_m, ~isr)
            R["step2_nonT4L_measurementB"] = s2b
            R["nonT4L_collapse_B"] = bool(pos(s2b))
            R["asymmetry_nonT4L"] = asymmetry(a, fam, o_m, c_m, ~isr)
            print(f"Step 2 non-T4L, measurement B: {fmt(s2b)} -> collapse: {R['nonT4L_collapse_B']}")
        else:
            R["stopped"] = "B not validated on T4L: cavity response not measurable with a validated method"
            print(R["stopped"])
    if "asymmetry_nonT4L" in R:
        s = R["asymmetry_nonT4L"]
        print(f"Asymmetry index (non-T4L, A scale): {s['index']:+.3f} {s['ci']} = overpack-null "
              f"{s['overpacking_minus_null']:+.3f} (n {s['n_overpacking']}) + cavity-null {s['cavity_minus_null']:+.3f} "
              f"(n {s['n_cavity']})")
    with open(os.path.join("results", "cavity.json"), "w") as fh:
        json.dump(R, fh, indent=1, default=float)
    print("wrote results/cavity.json")


if __name__ == "__main__":
    main()
