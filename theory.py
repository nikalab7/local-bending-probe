"""
proteinX -- physics-derived predictions of local backbone response (THEORY.md).

Exploratory; hypotheses, measurements and decision rules were committed in
THEORY.md before this was run. Data: dev + first-lockbox rows (final_eval
caches), frozen label unchanged.

  python theory.py   -> results/theory.json
"""
from __future__ import annotations
import json
import os
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import spearmanr

import delta_model as dm
import final_eval as fe
import structure_features as sf
from mechanisms import path_of, chains_of, atoms, at, kabsch, window_geometry

FAR, CONTACT, MIN_FAR, MIN_ATOMS = 10.0, 7.0, 30, 3
N_NULL, SEED, N_BOOT = 1500, 0, 2000
BIN = 10


# ------------------------------ bootstrap -----------------------------------
def _boot(fam, seed=SEED):
    fams = np.unique(fam); by = {f: np.flatnonzero(fam == f) for f in fams}
    rng = np.random.default_rng(seed)
    for _ in range(N_BOOT):
        yield np.concatenate([by[f] for f in rng.choice(fams, len(fams))])


def q95(v):
    v = [x for x in v if np.isfinite(x)]
    return [float(x) for x in np.percentile(v, [2.5, 97.5])] if v else None


def mean_ci(v, fam, m):
    m = m & np.isfinite(v)
    if m.sum() == 0:
        return dict(n=0)
    return dict(n=int(m.sum()), families=int(len(set(fam[m]))), mean=float(v[m].mean()),
                ci=q95([v[i][m[i]].mean() for i in _boot(fam) if m[i].any()]))


def diff_ci(v, fam, ma, mb):
    ma, mb = ma & np.isfinite(v), mb & np.isfinite(v)
    if not ma.any() or not mb.any():
        return dict(n_a=int(ma.sum()), n_b=int(mb.sum()))
    bs = [v[i][ma[i]].mean() - v[i][mb[i]].mean() for i in _boot(fam) if ma[i].any() and mb[i].any()]
    return dict(n_a=int(ma.sum()), n_b=int(mb.sum()), a=float(v[ma].mean()), b=float(v[mb].mean()),
                diff=float(v[ma].mean() - v[mb].mean()), ci=q95(bs))


def spearman_ci(x, y, fam, m):
    m = m & np.isfinite(x) & np.isfinite(y)
    rho = float(spearmanr(x[m], y[m])[0])
    bs = [spearmanr(x[i][m[i]], y[i][m[i]])[0] for i in _boot(fam) if m[i].sum() > 5]
    return dict(n=int(m.sum()), families=int(len(set(fam[m]))), rho=rho, ci=q95(bs))


# ------------------------------ radial displacement -------------------------
def radial_one(wt, mt, r):
    cb = at(wt, r, "CB"); ca_r = at(wt, r, "CA")
    if cb is None or ca_r is None:
        return None
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
            if np.linalg.norm(xw - cb) > CONTACT:
                continue
            xm = mt[k][1].get(name)
            if xm is None:
                continue
            u = (xw - cb) / np.linalg.norm(xw - cb)
            vals.append(float(np.dot(xm @ R + t - xw, u)))
    return float(np.mean(vals)) if len(vals) >= MIN_ATOMS else None


def radial_row(x, expect_aa):
    wt = atoms(x["scaffold_path"], x["scaffold_chain"])
    out = []
    for pdb in x["pdbs"]:
        path = path_of(pdb)
        if path is None:
            continue
        best, best_rmsd = None, np.inf
        for ch in chains_of(path):
            mt = atoms(path, ch)
            if mt.get(x["r"], ("?",))[0] != expect_aa:
                continue
            g = window_geometry(wt, mt, x["s"])
            if g and g["ca_rmsd"] < best_rmsd:
                best, best_rmsd = mt, g["ca_rmsd"]
        if best is not None:
            v = radial_one(wt, best, x["r"])
            if v is not None:
                out.append(v)
    return float(np.median(out)) if out else np.nan


# ------------------------------ Ramachandran --------------------------------
def rclass(aa, nxt):
    if aa == "G":
        return "G"
    if aa == "P":
        return "P"
    if nxt == "P":
        return "prePro"
    if aa in "IV":
        return "IV"
    return "general"


def rama_reference(rows):
    nb = 360 // BIN
    H = {c: np.zeros((nb, nb)) for c in ("G", "P", "prePro", "IV", "general")}
    seen = set()
    for x in rows:
        k = (x["scaffold_path"], x["scaffold_chain"])
        if k in seen:
            continue
        seen.add(k)
        res = atoms(*k)
        for i, (aa, _) in res.items():
            phi, psi = sf.phi_psi(res, i)
            if np.isfinite(phi) and np.isfinite(psi):
                nxt = res.get(i + 1, ("?",))[0]
                H[rclass(aa, nxt)][int((phi + 180) // BIN) % nb, int((psi + 180) // BIN) % nb] += 1
    L = {}
    for c, h in H.items():
        h = gaussian_filter(h, 1.0, mode="wrap")
        h = h + 1e-4 * h.sum()
        L[c] = np.log(h / h.sum())
    return L, {c: int(h.sum()) for c, h in H.items()}


def logp(L, c, phi, psi):
    nb = 360 // BIN
    return L[c][int((phi + 180) // BIN) % nb, int((psi + 180) // BIN) % nb]


def strain(L, x):
    res = atoms(x["scaffold_path"], x["scaffold_chain"])
    r = x["r"]
    phi, psi = sf.phi_psi(res, r)
    if not (np.isfinite(phi) and np.isfinite(psi)):
        return np.nan
    nxt = res.get(r + 1, ("?",))[0]
    S = logp(L, rclass(x["mut"], nxt), phi, psi) - logp(L, rclass(x["wt"], nxt), phi, psi)
    if "P" in (x["wt"], x["mut"]) and (r - 1) in res:
        pf, ps = sf.phi_psi(res, r - 1)
        prev = res[r - 1][0]
        if np.isfinite(pf) and np.isfinite(ps):
            S += logp(L, rclass(prev, x["mut"]), pf, ps) - logp(L, rclass(prev, x["wt"]), pf, ps)
    return float(S)


# ------------------------------ main ----------------------------------------
def structural_class(x):
    ch = x["scaffold_chain"]; n = len(x["scaffold_res"])
    h = sum(1 for c, _ in x["scaffold_helix"] if c == ch) / max(n, 1)
    e = sum(1 for c, _ in x["scaffold_sheet"] if c == ch) / max(n, 1)
    return "all-alpha" if h >= 0.4 and e < 0.1 else ("all-beta" if e >= 0.3 and h < 0.1 else "alpha/beta")


def main():
    rd, nd, Fd = fe.load("dev"); rl, nl, Fl = fe.load("lockbox")
    rows, F = rd + rl, Fd + Fl
    n = len(rows)
    fam = np.array([str(x["family"]) for x in rows])
    z = np.abs([x["z"] for x in rows]); y = z > 2
    hse = dm.matrix(F, ["hse_up"])[:, 0]
    dv = np.array([dm.VOL[x["mut"]] - dm.VOL[x["wt"]] for x in rows], float)
    gp = np.array([x["wt"] in "GP" or x["mut"] in "GP" for x in rows])
    R = dict(note="exploratory; pre-specified in THEORY.md", n=n)
    # ---- radial displacement, real rows and null rows
    rad = np.array([radial_row(x, x["mut"]) for x in rows])
    null = nd + nl
    pick = np.random.default_rng(SEED).choice(len(null), min(N_NULL, len(null)), replace=False)
    nrows = [null[i] for i in pick]
    nrad = np.array([radial_row(x, x["wt"]) for x in nrows])
    nfam = np.array([str(x["family"]) for x in nrows])
    print(f"radial: rows {np.isfinite(rad).sum()}/{n}, null {np.isfinite(nrad).sum()}/{len(nrows)}")
    v = np.concatenate([rad, nrad]); allf = np.concatenate([fam, nfam])
    isreal = np.r_[np.ones(n, bool), np.zeros(len(nrows), bool)]
    pad = lambda m: np.r_[m, np.zeros(len(nrows), bool)]
    nullm = ~isreal
    for tag, q in (("primary", 50), ("sensitivity", 67)):
        b = hse >= np.nanpercentile(hse, q)
        G = {"overpacking": b & (dv >= 25) & ~gp, "cavity": b & (dv <= -25) & ~gp,
             "neutral": b & (np.abs(dv) < 10) & ~gp}
        out = {"null": mean_ci(v, allf, nullm)}
        for g, m in G.items():
            out[g] = mean_ci(v, allf, pad(m))
            out[f"{g}_minus_null"] = diff_ci(v, allf, pad(m), nullm)
        out["T3_spearman_dvol_radial"] = spearman_ci(dv, rad, fam, b & ~gp)
        R[f"radial_{tag}"] = out
    P = R["radial_primary"]
    ok = lambda d, sgn: d.get("ci") is not None and (d["ci"][0] > 0 if sgn > 0 else d["ci"][1] < 0)
    R["T1_supported"] = bool(ok(P["overpacking"], +1) and ok(P["overpacking_minus_null"], +1))
    R["T2_supported"] = bool(ok(P["cavity"], -1) and ok(P["cavity_minus_null"], -1))
    R["T3_supported"] = bool(ok(P["T3_spearman_dvol_radial"], +1))
    for g in ("null", "overpacking", "cavity", "neutral"):
        d = P[g]
        print(f"  [{g:11s}] mean radial {d['mean']:+.3f} A [{d['ci'][0]:+.3f},{d['ci'][1]:+.3f}] n={d['n']} fam={d['families']}")
    for g in ("overpacking", "cavity", "neutral"):
        d = P[f"{g}_minus_null"]
        print(f"  {g:11s} - null {d['diff']:+.3f} [{d['ci'][0]:+.3f},{d['ci'][1]:+.3f}]")
    s = P["T3_spearman_dvol_radial"]
    print(f"  T3 Spearman(dvol, radial) {s['rho']:+.3f} [{s['ci'][0]:+.3f},{s['ci'][1]:+.3f}] n={s['n']}")
    print(f"  T1 {R['T1_supported']}  T2 {R['T2_supported']}  T3 {R['T3_supported']}")
    # ---- Ramachandran strain
    L, counts = rama_reference(rows)
    S = np.array([strain(L, x) for x in rows])
    dih = np.array([max(abs(x["z_phi"]), abs(x["z_psi"])) for x in rows])
    allm = np.ones(n, bool)
    T4 = dict(reference_counts=counts, S_quantiles=[float(np.nanpercentile(S, q)) for q in (5, 25, 50, 75, 95)],
              spearman_negS_dihedral=spearman_ci(-S, dih, fam, allm),
              spearman_negS_combined_z=spearman_ci(-S, z, fam, allm),
              mover_strained_vs_unstrained=diff_ci(y.astype(float), fam, S <= -2, S > -0.5))
    folds = dm.site_folds(fam, dm.FAMILY_FOLDS, 0)
    w = dm.label_weights(np.array([x["z"] for x in rows]), 2.0)
    Xb = dm.matrix(F, fe.FEATURES)
    base = fe.oof(Xb, y, folds, w)
    withS = fe.oof(np.column_stack([Xb, S]), y, folds, w)
    T4["dAUC_site_subst_plus_S"] = fe.paired(y, base, withS, fam)
    R["T4"] = T4
    sp, st = T4["spearman_negS_dihedral"], T4["mover_strained_vs_unstrained"]
    R["T4_supported"] = bool(sp["ci"][0] > 0 and st.get("ci") and st["ci"][0] > 0)
    print(f"  T4 Spearman(-S, dihedral |z|) {sp['rho']:+.3f} [{sp['ci'][0]:+.3f},{sp['ci'][1]:+.3f}] n={sp['n']}")
    print(f"  T4 mover strained (S<=-2) {st.get('a', float('nan')):.3f} vs unstrained {st.get('b', float('nan')):.3f} "
          f"diff {st.get('diff', float('nan')):+.3f} {st.get('ci')} (n {st['n_a']} vs {st['n_b']})")
    d = T4["dAUC_site_subst_plus_S"]
    print(f"  T4 dAUC (+S) {d['d_auc']:+.3f} [{d['lo']:+.3f},{d['hi']:+.3f}]   T4 supported: {R['T4_supported']}")
    # ---- conditional T5 / T6
    any_sup = any(R[f"T{i}_supported"] for i in (1, 2, 3, 4))
    R["conditional_run"] = any_sup
    if any_sup:
        zs = lambda a: (a - np.nanmean(a)) / np.nanstd(a)
        strain_t = np.nan_to_num(np.maximum(-S, 0), nan=0.0)
        Pphys = zs(strain_t) + zs(np.abs(dv) * np.nan_to_num(hse, nan=np.nanmedian(hse)) / 20)
        C = zs(np.nan_to_num(dm.matrix(F, ["b_window"])[:, 0], nan=0.0))
        score = Pphys + C
        ss = fe.oof(dm.matrix(F, dm.SS), y, folds, w)
        R["T5"] = dict(auc=fe.auc_ci(y, score, fam), spearman_abs_z=spearman_ci(score, z, fam, allm),
                       vs_ss_only=fe.paired(y, ss, score, fam), frozen_model_oof=fe.auc_ci(y, base, fam))
        a = R["T5"]["auc"]
        print(f"  T5 composite AUC {a['auc']:.3f} [{a['lo']:.3f},{a['hi']:.3f}] vs frozen model OOF "
              f"{R['T5']['frozen_model_oof']['auc']:.3f}; vs SS-only {R['T5']['vs_ss_only']['d_auc']:+.3f}")
        cls = np.array([structural_class(x) for x in rows])
        ncls = np.array([structural_class(x) for x in nrows])
        allcls = np.r_[cls, ncls]
        R["T6"] = {}
        for c in ("all-alpha", "all-beta", "alpha/beta"):
            mc_, mn_ = pad(cls == c), (allcls == c) & nullm
            b = hse >= np.nanpercentile(hse, 50)
            R["T6"][c] = dict(
                n_rows=int((cls == c).sum()),
                overpacking_minus_null=diff_ci(v, allf, mc_ & pad(b & (dv >= 25) & ~gp), mn_),
                cavity_minus_null=diff_ci(v, allf, mc_ & pad(b & (dv <= -25) & ~gp), mn_),
                T3_spearman=spearman_ci(dv, rad, fam, (cls == c) & b & ~gp),
                T4_spearman=spearman_ci(-S, dih, fam, cls == c))
            t = R["T6"][c]
            print(f"  T6 {c:10s} rows={t['n_rows']:4d} overpack-null {t['overpacking_minus_null'].get('diff', float('nan')):+.3f} "
                  f"(n {t['overpacking_minus_null']['n_a']}) cavity-null {t['cavity_minus_null'].get('diff', float('nan')):+.3f} "
                  f"(n {t['cavity_minus_null']['n_a']}) T3 {t['T3_spearman']['rho']:+.2f} T4 {t['T4_spearman']['rho']:+.2f}")
    with open(os.path.join("results", "theory.json"), "w") as fh:
        json.dump(R, fh, indent=1, default=float)
    print("wrote results/theory.json")


if __name__ == "__main__":
    main()
