"""
proteinX -- DELTA-TARGETED mover model on clean labels (pairs.py).

Why a new engine
----------------
Gate 2 trained a regressor on ABSOLUTE bending (std ~38 deg, RMSE ~30 deg) and
then differenced two predictions that differ by one one-hot position. Movers
change by a median of ~2.7 deg, an order of magnitude below the model's
resolution, and a tree ensemble returns exactly 0 whenever no split touches
the central residue. That engine could not see the effect even if it existed.

This model instead predicts the label directly -- "does THIS mutation at THIS
site move the backbone" -- from two feature groups computed on the WT scaffold:

  SITE   where the mutation is: secondary structure, C-alpha B-factor (z within
         chain), burial (C-alpha neighbour count, half-sphere exposure along a
         C-alpha-only pseudo side-chain direction), non-local contacts, WT
         window bending, distance to a terminus.
  SUBST  what the mutation is: volume / hydrophobicity / charge change,
         to/from Pro and Gly, BLOSUM62, and two burial interactions
         (cavity creation = volume loss x burial; hydrophobicity change x burial).

The scientific question becomes a clean ablation:
  * SITE vs chance      -> is "where" predictive?
  * SITE - SS           -> is "where" more than secondary structure?
  * SITE+SUBST - SITE   -> does the substitution identity (the local-sequence
                           part) add anything once the site is known?
All are scored out-of-fold with LEAVE-SITE-OUT CV (no residue appears in both
train and test) and LEAVE-FAMILY-OUT CV (whole >= 30%-identity families held
out), with residue-cluster bootstrap CIs and paired delta-AUC tests.

NULL CONTROL. The same site model is trained on WT-vs-WT pseudo-mutants
(pairs.null_rows), where every "mover" is noise. If site features predict
those, the real-data site AUC can be noise structure, not movement.

Usage: python delta_model.py [--z 2.0] [--prior bfactor|pooled]
"""
from __future__ import annotations
import argparse
import json
import os
import numpy as np
from stats_utils import cluster_auc_ci, paired_delta_auc, within_auc_ci
from structure_features import CONTEXT, INTERACT

AA = "ARNDCQEGHILKMFPSTWYV"
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)                                   # Kyte-Doolittle
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4,
           G=60.1, H=153.2, I=166.7, L=166.7, K=168.6, M=162.9, F=189.9,
           P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)  # Zamyatnin, A^3
CHARGE = dict(D=-1, E=-1, K=1, R=1)
_BLOSUM62 = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
_B = np.array(_BLOSUM62.split(), dtype=int).reshape(20, 20)
BLOSUM62 = {(a, b): int(_B[i, j]) for i, a in enumerate(AA) for j, b in enumerate(AA)}

SITE = ["ss_H", "ss_E", "ss_L", "b_site", "b_window", "n_ca10", "hse_up",
        "nonlocal_contacts", "wt_bend", "log_term_dist"]
SUBST = ["d_vol", "abs_d_vol", "d_hyd", "d_charge", "to_P", "from_P", "to_G",
         "from_G", "blosum62", "cavity", "d_hyd_x_burial"]
SS = ["ss_H", "ss_E", "ss_L"]
FEATURE_SETS = {
    "ss": SS,
    "site": SITE,                                    # C-alpha "where" (v2/v3)
    "subst": SUBST,                                  # context-free "what"
    "site+subst": SITE + SUBST,
    "where": SITE + CONTEXT,                         # + full-atom site context
    "where+what": SITE + CONTEXT + SUBST + INTERACT,  # + substitution x context
}
# paired contrasts (base, richer, name)
CONTRASTS = [("site", "site+subst", "subst_over_site"),   # context-free "what" over "where"
             ("ss", "site", "site_over_ss"),              # is "where" more than SS?
             ("site", "where", "context_over_site"),      # does full-atom context help?
             ("where", "where+what", "what_over_where")]  # "what given where" over "where"
HGB_SETS = ("site", "where", "where+what")             # boosting only where informative
FAMILY_FOLDS = 10


# -------------------------------- features ----------------------------------
def site_features(res, r, s, ss, wt_bend):
    """Site descriptors from WT scaffold residues {resSeq: (aa, ca, bfactor)}."""
    nums = np.array(sorted(res))
    ca = np.array([res[n][1] for n in nums])
    bf = np.array([res[n][2] for n in nums], float)
    fin = np.isfinite(bf)
    bz = (bf - bf[fin].mean()) / (bf[fin].std() + 1e-9) if fin.sum() > 2 else bf * np.nan
    idx = {int(n): i for i, n in enumerate(nums)}
    k = idx[r]
    d = np.linalg.norm(ca - ca[k], axis=1)
    sep = np.abs(nums - r)
    # C-alpha-only pseudo side-chain direction (bisector away from neighbours)
    hse = np.nan
    if r - 1 in idx and r + 1 in idx:
        u = (ca[k] - ca[idx[r - 1]]) + (ca[k] - ca[idx[r + 1]])
        if np.linalg.norm(u) > 1e-6:
            u /= np.linalg.norm(u)
            near = (d < 13.0) & (sep > 0)
            hse = float((((ca[near] - ca[k]) @ u) > 0).sum())
    win = [idx[w] for w in range(s, s + 5) if w in idx]
    return dict(
        ss_H=float(ss == "H"), ss_E=float(ss == "E"), ss_L=float(ss == "L"),
        b_site=float(bz[k]),
        b_window=float(np.nanmean(bz[win])) if win else np.nan,
        n_ca10=float(((d < 10.0) & (sep > 0)).sum()),
        hse_up=hse,
        nonlocal_contacts=float(((d < 8.0) & (sep > 4)).sum()),
        wt_bend=float(wt_bend),
        log_term_dist=float(np.log1p(min(r - nums.min(), nums.max() - r))))


def subst_features(wt, mut, hse_up):
    burial = 0.0 if not np.isfinite(hse_up) else hse_up / 20.0
    dv = (VOL[mut] - VOL[wt]) / 100.0
    dh = KD[mut] - KD[wt]
    return dict(
        d_vol=dv, abs_d_vol=abs(dv), d_hyd=dh,
        d_charge=float(CHARGE.get(mut, 0) - CHARGE.get(wt, 0)),
        to_P=float(mut == "P"), from_P=float(wt == "P"),
        to_G=float(mut == "G"), from_G=float(wt == "G"),
        blosum62=float(BLOSUM62[(wt, mut)]),
        cavity=max(0.0, -dv) * burial, d_hyd_x_burial=dh * burial)


def featurize(rows):
    """Site + substitution features; full-atom CONTEXT/INTERACT when the row
    carries its scaffold file (pairs rows do), NaN otherwise (imputed)."""
    import structure_features
    F = []
    for x in rows:
        sf = site_features(x["scaffold_res"], x["r"], x["s"], x["scaffold_ss"], x["wt_bend"])
        f = {**sf, **subst_features(x["wt"], x["mut"], sf["hse_up"])}
        if "scaffold_path" in x:
            f.update(structure_features.row_features(x))
        else:
            f.update({k: np.nan for k in CONTEXT + INTERACT})
        F.append(f)
    return F


def matrix(F, cols):
    return np.array([[f[c] for c in cols] for f in F], float)


# ------------------------------- models / CV ---------------------------------
def make_model(kind):
    from sklearn.pipeline import make_pipeline
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import HistGradientBoostingClassifier
    if kind == "logreg":
        return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                             LogisticRegression(C=0.3, class_weight="balanced",
                                                max_iter=2000))
    return HistGradientBoostingClassifier(max_depth=3, max_iter=200, learning_rate=0.05,
                                          min_samples_leaf=10, class_weight="balanced",
                                          random_state=0)


def oof_predict(X, y, fold_of, kind):
    """Out-of-fold probabilities; folds whose train split is one-class are NaN."""
    p = np.full(len(y), np.nan)
    for f in np.unique(fold_of):
        te = fold_of == f; tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        m = make_model(kind).fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[:, 1]
    return p


def site_folds(groups, k, seed):
    """Assign whole (protein, residue) groups to k folds at random."""
    labs = sorted(set(groups))
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(labs))
    fold = {labs[i]: j % k for j, i in enumerate(order)}
    return np.array([fold[g] for g in groups])


def evaluate(F, y, groups, families, k=5, repeats=5, feature_sets=None,
             kinds=("logreg", "hgb"), proteins=None):
    """Out-of-fold AUCs + paired contrasts.

    leave_site_out   : (protein, residue) groups in k folds, averaged over repeats
    leave_family_out : whole sequence families (>= 30% identity) in
                       FAMILY_FOLDS folds (leave-one-family-out if fewer)
    within-protein AUC (if proteins given): only mover/non-mover pairs from the
    same protein are compared, so protein-level base rates cannot help.
    """
    from sklearn.metrics import average_precision_score
    feature_sets = feature_sets or FEATURE_SETS
    families = np.asarray(families)
    n_fam = len(set(families))
    results, preds = {}, {}
    for name, cols in feature_sets.items():
        X = matrix(F, cols)
        for kind in kinds:
            if kind == "hgb" and name not in HGB_SETS:
                continue
            P = np.nanmean([oof_predict(X, y, site_folds(groups, k, s), kind)
                            for s in range(repeats)], axis=0)
            Q = oof_predict(X, y, site_folds(families, min(FAMILY_FOLDS, n_fam), 0), kind) \
                if n_fam > 1 else np.full(len(y), np.nan)
            preds[(name, kind)] = (P, Q)
            res = {}
            for cv, pr in (("leave_site_out", P), ("leave_family_out", Q)):
                ok = np.isfinite(pr)
                if ok.sum() < 10 or not (0 < y[ok].sum() < ok.sum()):
                    continue
                auc, lo, hi, ncl = cluster_auc_ci(y[ok], pr[ok], groups[ok])
                res[cv] = dict(auc=auc, lo=lo, hi=hi, n=int(ok.sum()),
                               movers=int(y[ok].sum()), clusters=ncl,
                               ap=float(average_precision_score(y[ok], pr[ok])),
                               base=float(y[ok].mean()))
                if proteins is not None:
                    w, wlo, whi = within_auc_ci(y[ok], pr[ok], np.asarray(proteins)[ok],
                                                groups[ok])
                    res[cv].update(within=w, within_lo=wlo, within_hi=whi)
            results[f"{name}|{kind}"] = res
    contrasts = {}
    for base, rich, label in CONTRASTS:
        for kind in kinds:
            if (base, kind) not in preds or (rich, kind) not in preds:
                continue
            for cv, j in (("leave_site_out", 0), ("leave_family_out", 1)):
                a = preds[(base, kind)][j]; b = preds[(rich, kind)][j]
                ok = np.isfinite(a) & np.isfinite(b)
                if ok.sum() >= 10 and 0 < y[ok].sum() < ok.sum():
                    d, lo, hi, p = paired_delta_auc(y[ok], a[ok], b[ok], groups[ok])
                    contrasts[f"{label}|{kind}|{cv}"] = dict(delta=d, lo=lo, hi=hi, p_le0=p)
    return results, contrasts


def null_control(null, zthr):
    """Site model on WT-vs-WT pseudo-mutants (pairs.null_rows).

    Every pseudo-mover is noise. Returns the pseudo-mover rate (calibration of
    the z threshold) and the leave-site-out AUC of the site and SS models on
    pseudo labels: an AUC above 0.5 means site features predict WHERE the
    noise is, which would also inflate the real-data site AUC.
    """
    if not null:
        return None
    y = np.array([abs(x["z"]) > zthr for x in null])
    g = np.array([f"{x['protein']}:{x['r']}" for x in null])
    fam = [x["family"] for x in null]
    out = dict(n=len(y), movers=int(y.sum()), rate=float(y.mean()),
               sites=len(set(g)))
    if y.sum() < 10:
        return out
    F = featurize(null)
    res, con = evaluate(F, y, g, fam, repeats=2, kinds=("logreg",),
                        feature_sets={"ss": SS, "site": SITE, "where": SITE + CONTEXT})
    out.update(results=res, contrasts=con)
    return out


def noise_decomposition(F, y, groups, Fn, yn, gn, cols, k=5, seed=0):
    """How much of a feature set's real-data AUC is predictability of noise?

    Real rows and WT-vs-WT pseudo-mutants share k site folds. For each fold, a
    "noise score" model is trained on the pseudo labels of the OTHER sites and
    scores the real rows; a "real" model is trained on real labels. Reported:
      noise_auc          AUC of the noise score on real labels
      real_auc           AUC of the real-label model
      excess             paired dAUC (noise score + features) - noise score:
                         the movement-specific part the features carry
      real_model_on_null AUC of the real-label model on pseudo labels
    """
    sites = sorted(set(groups) | set(gn))
    perm = np.random.default_rng(seed).permutation(len(sites))
    fold = {sites[j]: i % k for i, j in enumerate(perm)}
    fr = np.array([fold[s] for s in groups]); fn = np.array([fold[s] for s in gn])
    X, Xn = matrix(F, cols), matrix(Fn, cols)
    noise = np.full(len(y), np.nan); real = np.full(len(y), np.nan)
    on_null = np.full(len(yn), np.nan)
    for f in range(k):
        tr_n, te = fn != f, fr == f
        if len(np.unique(yn[tr_n])) == 2:
            noise[te] = make_model("logreg").fit(Xn[tr_n], yn[tr_n]).predict_proba(X[te])[:, 1]
        m = make_model("logreg").fit(X[fr != f], y[fr != f])
        real[te] = m.predict_proba(X[te])[:, 1]
        on_null[fn == f] = m.predict_proba(Xn[fn == f])[:, 1]
    ok = np.isfinite(noise)
    lg = np.log(np.clip(noise, 1e-6, 1 - 1e-6) / (1 - np.clip(noise, 1e-6, 1 - 1e-6)))[:, None]
    base = np.full(len(y), np.nan); both = np.full(len(y), np.nan)
    for f in range(k):
        tr, te = (fr != f) & ok, (fr == f) & ok
        base[te] = make_model("logreg").fit(lg[tr], y[tr]).predict_proba(lg[te])[:, 1]
        Z = np.c_[lg, X]
        both[te] = make_model("logreg").fit(Z[tr], y[tr]).predict_proba(Z[te])[:, 1]
    out = {}
    a, lo, hi, _ = cluster_auc_ci(y[ok], noise[ok], groups[ok]); out["noise_auc"] = (a, lo, hi)
    a, lo, hi, _ = cluster_auc_ci(y, real, groups); out["real_auc"] = (a, lo, hi)
    d, lo, hi, p = paired_delta_auc(y[ok], base[ok], both[ok], groups[ok])
    out["excess"] = (d, lo, hi, p)
    a, lo, hi, _ = cluster_auc_ci(yn, on_null, gn); out["real_model_on_null"] = (a, lo, hi)
    return out


# ---------------------------------- main -------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--z", type=float, default=None,
                    help="mover threshold on |delta|/SE (default: pairs.Z_MOVER)")
    ap.add_argument("--prior", choices=("bfactor", "pooled"), default=None,
                    help="noise-floor sigma prior (default: pairs.PRIOR)")
    args = ap.parse_args()

    from pairs import load_proteins, Z_MOVER, PRIOR
    zthr = args.z if args.z is not None else Z_MOVER
    prior = args.prior or PRIOR
    null = []
    data = load_proteins(prior, null_out=null)
    rows = [x for rs, _ in data.values() for x in rs]
    used = sorted(((d.get("mutations", 0), d.get("movers", 0), n) for n, (_, d) in data.items()
                   if d.get("mutations", 0)), reverse=True)
    print(f"proteins with >= 1 clean mutation: {len(used)} / {len(data)}; largest:")
    for m, mv, name in used[:12]:
        print(f"  {name:16s} mutations={m:4d} movers={mv:3d}")
    if not rows:
        print("no clean pairs -- run the gate scripts first to fill the PDB caches")
        return
    y = np.array([abs(x["z"]) > zthr for x in rows])
    groups = np.array([f"{x['protein']}:{x['r']}" for x in rows])
    proteins = [x["protein"] for x in rows]
    families = [x["family"] for x in rows]
    print(f"pooled: n={len(y)} mutations, {int(y.sum())} movers (|z|>{zthr}), "
          f"{len(set(groups))} sites, {len(set(proteins))} proteins, "
          f"{len(set(families))} families")
    if y.sum() < 10 or (~y).sum() < 10:
        print("too few movers or non-movers to evaluate"); return

    F = featurize(rows)
    results, contrasts = evaluate(F, y, groups, families, proteins=proteins)
    nc = null_control(null, zthr)
    decomp = {}
    if nc and nc.get("movers", 0) >= 10:
        Fn = featurize(null)
        yn = np.array([abs(x["z"]) > zthr for x in null])
        gn = np.array([f"{x['protein']}:{x['r']}" for x in null])
        for name in ("site", "where"):
            decomp[name] = noise_decomposition(F, y, groups, Fn, yn, gn, FEATURE_SETS[name])

    print("\n" + "=" * 74 + "\n  DELTA MODEL -- out-of-fold mover retrieval\n" + "=" * 74)
    for key, res in results.items():
        for cv, m in res.items():
            w = (f"  within-protein={m['within']:.3f} [{m['within_lo']:.2f},{m['within_hi']:.2f}]"
                 if "within" in m else "")
            print(f"  {key:20s} {cv:17s} AUC={m['auc']:.3f} 90%CI[{m['lo']:.2f},{m['hi']:.2f}] "
                  f"AP={m['ap']:.3f} (base {m['base']:.2f}){w}")
    print("\n  paired dAUC (richer - base): subst_over_site, site_over_ss, "
          "context_over_site, what_over_where")
    for key, c in contrasts.items():
        print(f"  {key:38s} dAUC={c['delta']:+.3f} 90%CI[{c['lo']:+.2f},{c['hi']:+.2f}] "
              f"P(d<=0)={c['p_le0']:.3f}")
    if nc:
        print(f"\n  NULL CONTROL (WT-vs-WT pseudo-mutants): n={nc['n']} at {nc['sites']} sites, "
              f"pseudo-movers={nc['movers']} ({nc['rate']:.1%}; real {y.mean():.1%})")
        for key, res in nc.get("results", {}).items():
            for cv, m in res.items():
                print(f"  {key:18s} {cv:17s} AUC={m['auc']:.3f} 90%CI[{m['lo']:.2f},{m['hi']:.2f}]")
    for name, d in decomp.items():
        print(f"\n  NOISE DECOMPOSITION ({name}): noise score on real labels AUC={d['noise_auc'][0]:.3f} "
              f"[{d['noise_auc'][1]:.2f},{d['noise_auc'][2]:.2f}]; real model AUC={d['real_auc'][0]:.3f}; "
              f"excess over noise dAUC={d['excess'][0]:+.3f} [{d['excess'][1]:+.2f},{d['excess'][2]:+.2f}] "
              f"P(d<=0)={d['excess'][3]:.3f}; real model on null AUC={d['real_model_on_null'][0]:.3f}")

    os.makedirs("results", exist_ok=True)
    with open("results/delta_model.json", "w") as fh:
        json.dump(dict(z_threshold=zthr, prior=prior, n=len(y), movers=int(y.sum()),
                       results=results, contrasts=contrasts, null_control=nc,
                       noise_decomposition=decomp), fh, indent=1)
    print("\nwrote results/delta_model.json")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    keys = [k for k in results if "leave_site_out" in results[k]]
    for i, k in enumerate(keys):
        m = results[k]["leave_site_out"]
        ax.errorbar(i, m["auc"], yerr=[[m["auc"] - m["lo"]], [m["hi"] - m["auc"]]],
                    fmt="o", capsize=5, color="#4C78A8" if "logreg" in k else "#F58518")
    ax.axhline(0.5, color="#E45756", ls="--", lw=1.2, label="chance")
    ax.set_xticks(range(len(keys))); ax.set_xticklabels(keys, rotation=20, fontsize=8)
    ax.set_ylabel("leave-site-out AUC (90% CI)")
    ax.set_title(f"Delta model on clean labels (n={len(y)}, movers={int(y.sum())})")
    ax.legend(); fig.tight_layout(); fig.savefig("delta_model.png", dpi=130)
    print("plot -> delta_model.png")


if __name__ == "__main__":
    main()
