"""
proteinX -- post-lockbox EXPLORATORY analyses on dev data only.

Run after the single lockbox evaluation (results/final_lockbox.json); it does
not touch lockbox data and none of it is confirmatory. It uses the frozen
label and the cached dev rows of final_eval.py.

  1. Bounds on "what": paired family-bootstrap dAUC (95%) for adding
     mutation-identity features to the site model, and for adding further
     site descriptors to the frozen site + substitution model. The upper CI
     bound is the largest gain the dev data are compatible with.
  2. Artifact rates: mover rate of crystallographic-mismatch groups vs the
     rest, with family-bootstrap 95% CIs on each rate and on the difference.

  python posthoc_dev.py [--esm-emb EMB.npy --esm-rows ROWS.pkl]
      -> results/posthoc_dev.json
"""
from __future__ import annotations
import argparse
import json
import pickle
import numpy as np

import delta_model as dm
import structure_features as sf
import final_eval as fe

ATTENUATION = (fe.CEILING - 0.5) / 0.5      # oracle margin: 0.574 of a perfect predictor's


def contrasts(F, y, w, fam, emb=None):
    folds = dm.site_folds(fam, dm.FAMILY_FOLDS, 0)
    base_site, base = dm.SITE, dm.SITE + dm.SUBST
    sets = [  # (group, name, base cols, added cols)
        ("what", "substitution (physicochemical)", base_site, dm.SUBST),
        ("what", "ESM-2 substitution LLR", base, dm.ESM_WHAT),
        ("what", "substitution x context terms", base, sf.INTERACT),
        ("what", "all mutation-identity features", base_site, dm.SUBST + dm.ESM_WHAT + sf.INTERACT),
        ("where", "ESM-2 site terms (entropy, WT log-prob)", base, dm.ESM_WHERE),
        ("where", "elastic network (ENM)", base, sf.ENM),
        ("where", "full-atom site context", base, sf.CONTEXT),
    ]
    out, cache = [], {}

    def score(cols):
        k = tuple(cols)
        if k not in cache:
            cache[k] = fe.oof(dm.matrix(F, cols), y, folds, w)
        return cache[k]
    for g, name, b, add in sets:
        r = fe.paired(y, score(b), score(b + add), fam)
        out.append(dict(group=g, name=name, base="site" if b == base_site else "site+subst", **r))
    if emb is not None:                       # ESM-2 embeddings: PCA fitted inside each fold
        from sklearn.decomposition import PCA
        emb, have = emb                       # rows without a saved embedding are dropped (both arms)
        F, y, w, fam, folds = [F[i] for i in np.flatnonzero(have)], y[have], w[have], fam[have], folds[have]
        emb = emb[have]
        Xb = dm.matrix(F, base)
        p = np.full(len(y), np.nan)
        for f in np.unique(folds):
            te, tr = folds == f, folds != f
            pca = PCA(16, random_state=0).fit(emb[tr])
            m = dm.fit_model("logreg", np.hstack([Xb[tr], pca.transform(emb[tr])]), y[tr], w[tr])
            p[te] = m.predict_proba(np.hstack([Xb[te], pca.transform(emb[te])]))[:, 1]
        p0 = fe.oof(Xb, y, folds, w)
        out.append(dict(group="where", name=f"ESM-2 per-residue embedding (16 PCs; n={len(y)})",
                        base="site+subst", **fe.paired(y, p0, p, fam)))
    for r in out:
        r["upper_bound_true_label_rough"] = r["hi"] / ATTENUATION
    return out


GP = ["to_P", "from_P", "to_G", "from_G"]


def subst_decomposition(rows, F, y, w, fam):
    """Where does the substitution gain come from? (site -> site + X, paired)."""
    folds = dm.site_folds(fam, dm.FAMILY_FOLDS, 0)
    S = dm.SITE
    gp_row = np.array([x["wt"] in "GP" or x["mut"] in "GP" for x in rows])
    out = []

    def add(name, mask, cols, yy=None):
        yy = y if yy is None else yy
        m = np.flatnonzero(mask)
        Fm, ym, wm, fm = [F[i] for i in m], yy[m], w[m], fam[m]
        f2 = folds[m]
        a = fe.oof(dm.matrix(Fm, S), ym, f2, wm)
        b = fe.oof(dm.matrix(Fm, S + cols), ym, f2, wm)
        r = fe.paired(ym, a, b, fm)
        out.append(dict(name=name, n=int(len(m)), movers=int(ym.sum()), **r))
    allr = np.ones(len(y), bool)
    add("all SUBST", allr, dm.SUBST)
    add("Pro/Gly flags only", allr, GP)
    add("SUBST without Pro/Gly flags", allr, [c for c in dm.SUBST if c not in GP])
    add("all SUBST, rows without G or P on either side", ~gp_row, dm.SUBST)
    add("all SUBST, rows with G or P", gp_row, dm.SUBST)
    zt = lambda k, t: np.abs(np.array([x.get(k, np.nan) for x in rows], float)) > t
    for k in ("z_bend", "z_phi", "z_psi", "z_ca_tor_b"):   # single-metric labels, same rows
        add(f"all SUBST, label = {k} alone (|z| > 2)", allr, dm.SUBST, zt(k, 2.0))
    return out


def rate_ci(y, fam, mask, seed=0):
    """Mover rates (group, rest, difference) with family-bootstrap 95% CIs."""
    def rates(idx):
        m = mask[idx]
        return (y[idx][m].mean() if m.any() else np.nan,
                y[idx][~m].mean() if (~m).any() else np.nan)
    a, b = rates(np.arange(len(y)))
    bs = np.array([rates(i) for i in fe._boot(fam, seed)])
    bs = bs[np.isfinite(bs).all(1)]
    q = lambda v: [float(x) for x in np.percentile(v, [2.5, 97.5])]
    return dict(n=int(mask.sum()), families=int(len(set(fam[mask]))), rate=float(a), rate_ci=q(bs[:, 0]),
                rest_n=int((~mask).sum()), rest_rate=float(b), rest_ci=q(bs[:, 1]),
                difference=float(a - b), difference_ci=q(bs[:, 0] - bs[:, 1]))


def artifacts(rows, y, fam):
    g = lambda k: np.array([x.get(k, np.nan) for x in rows], float)
    tw, tm, rw, rm = g("temp_wt"), g("temp_mut"), g("res_wt"), g("res_mut")
    with np.errstate(invalid="ignore"):
        groups = {"mutant room temperature (>250 K), WT cryo (<150 K)": (tm > 250) & (tw < 150),
                  "temperature mismatch > 50 K": np.abs(tm - tw) > 50,
                  "mutant resolution worse by > 0.5 A": (rm - rw) > 0.5,
                  "mutant resolution worse by > 0.3 A": (rm - rw) > 0.3}
    return {k: rate_ci(y, fam, m) for k, m in groups.items()}


def load_embeddings(rows, emb_path, emb_rows_path, key="bend+ncs"):
    """Map saved per-row ESM-2 embeddings (scratch run) onto the dev rows."""
    E = np.load(emb_path)
    src = pickle.load(open(emb_rows_path, "rb"))[key]["rows"]
    idx = {(x["protein"], x["scaffold"], x["r"]): i for i, x in enumerate(src)}
    keys = [(x["protein"], x["scaffold"], x["r"]) for x in rows]
    have = np.array([k in idx for k in keys])
    out = np.full((len(rows), E.shape[1]), np.nan, np.float32)
    out[have] = E[[idx[k] for k, h in zip(keys, have) if h]]
    return out, have


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--esm-emb"); ap.add_argument("--esm-rows")
    a = ap.parse_args()
    rows, _, F = fe.load("dev")
    A = fe.arrays(rows)
    emb = load_embeddings(rows, a.esm_emb, a.esm_rows) if a.esm_emb else None
    R = dict(note="exploratory, dev only, after the lockbox run", attenuation=ATTENUATION,
             contrasts=contrasts(F, A["y"], A["w"], A["fam"], emb),
             substitution_gain=subst_decomposition(rows, F, A["y"], A["w"], A["fam"]),
             artifacts=artifacts(rows, A["y"], A["fam"]))
    with open("results/posthoc_dev.json", "w") as fh:
        json.dump(R, fh, indent=1)
    for c in R["contrasts"]:
        print(f"  [{c['group']}] +{c['name']:42s} over {c['base']:10s} dAUC {c['d_auc']:+.3f} "
              f"[{c['lo']:+.3f},{c['hi']:+.3f}]  upper/attenuation {c['upper_bound_true_label_rough']:+.3f}")
    for c in R["substitution_gain"]:
        print(f"  site -> site + {c['name']:52s} n={c['n']:4d} movers={c['movers']:3d} dAUC {c['d_auc']:+.3f} "
              f"[{c['lo']:+.3f},{c['hi']:+.3f}] p={c['p_two_sided']:.3g}")
    for k, v in R["artifacts"].items():
        print(f"  {k:52s} n={v['n']:4d} ({v['families']} fam) {v['rate']:.3f} "
              f"[{v['rate_ci'][0]:.3f},{v['rate_ci'][1]:.3f}] vs rest {v['rest_rate']:.3f} "
              f"[{v['rest_ci'][0]:.3f},{v['rest_ci'][1]:.3f}]  diff {v['difference']:+.3f} "
              f"[{v['difference_ci'][0]:+.3f},{v['difference_ci'][1]:+.3f}]")
    print("wrote results/posthoc_dev.json")


if __name__ == "__main__":
    main()
