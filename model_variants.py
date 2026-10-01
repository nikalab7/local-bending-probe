"""
proteinX -- exploratory model variants and a learning curve on the clean labels.

delta_model.py reports the pre-specified models (logistic regression, plus
shallow boosting on three feature sets). This script asks the follow-up
"can a different learner or more data buy accuracy?" with the same
leave-family-out folds, every variant paired against the baseline:

  baseline logreg C = 0.3       (delta_model's model)
  nested logreg                 C chosen by inner leave-site-out CV per fold
  ridge on |z|                  continuous target instead of the binary label
  hgb classifier                strongly regularized boosting (depth 2)
  hgb regressor on |z|
  ensemble                      rank(logreg) + rank(hgb classifier)

LEARNING CURVE: family-out AUC when only 25 / 50 / 75 / 100% of the training
families of each fold are used (5 random draws per fraction).

Picking the best of these variants after the fact is optimistic by a small
amount (6 variants x 2 feature sets); treat the best number as an upper
estimate, not a validated model.

Usage: python model_variants.py   -> results/model_variants.json
"""
from __future__ import annotations
import json
import os
import numpy as np
from scipy.stats import rankdata
import delta_model as dm
from stats_utils import cluster_auc_ci, paired_delta_auc, within_auc, fast_auc

FEATURES = ("site", "where+what")
FRACTIONS = (0.25, 0.5, 0.75, 1.0)


def _pipe(est):
    from sklearn.pipeline import make_pipeline
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), est)


def logreg(C=0.3):
    from sklearn.linear_model import LogisticRegression
    return _pipe(LogisticRegression(C=C, class_weight="balanced", max_iter=3000))


def ridge(alpha=10.0):
    from sklearn.linear_model import Ridge
    return _pipe(Ridge(alpha=alpha))


def hgb(regressor=False):
    from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
    kw = dict(max_depth=2, max_iter=150, learning_rate=0.05, min_samples_leaf=40,
              l2_regularization=1.0, random_state=0)
    return HistGradientBoostingRegressor(**kw) if regressor else \
        HistGradientBoostingClassifier(class_weight="balanced", **kw)


def oof(make, X, target, folds, clf=True):
    p = np.full(len(target), np.nan)
    for f in np.unique(folds):
        te, tr = folds == f, folds != f
        m = make().fit(X[tr], target[tr])
        p[te] = m.predict_proba(X[te])[:, 1] if clf else m.predict(X[te])
    return p


def nested_logreg(X, y, groups, folds, grid=(0.01, 0.03, 0.1, 0.3, 1.0)):
    p, chosen = np.full(len(y), np.nan), []
    for f in np.unique(folds):
        te, tr = folds == f, folds != f
        inner = dm.site_folds(groups[tr], 5, 1)
        scores = []
        for C in grid:
            q = np.full(tr.sum(), np.nan)
            for k in range(5):
                a, b = inner != k, inner == k
                q[b] = logreg(C).fit(X[tr][a], y[tr][a]).predict_proba(X[tr][b])[:, 1]
            scores.append(fast_auc(y[tr], q))
        C = grid[int(np.argmax(scores))]
        chosen.append(C)
        p[te] = logreg(C).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return p, chosen


def learning_curve(X, y, fam, prot, folds):
    out = []
    for frac in FRACTIONS:
        aucs, wins = [], []
        for rep in range(5 if frac < 1 else 1):
            rng = np.random.default_rng(rep)
            p = np.full(len(y), np.nan)
            for f in np.unique(folds):
                te = folds == f
                fams = sorted(set(fam[~te]))
                keep = rng.choice(fams, max(2, int(round(frac * len(fams)))), replace=False)
                tr = ~te & np.isin(fam, keep)
                if len(np.unique(y[tr])) == 2:
                    p[te] = logreg().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
            ok = np.isfinite(p)
            aucs.append(fast_auc(y[ok], p[ok])); wins.append(within_auc(y[ok], p[ok], prot[ok]))
        out.append(dict(fraction=frac, auc=float(np.mean(aucs)), auc_sd=float(np.std(aucs)),
                        within=float(np.mean(wins))))
    return out


def main():
    from pairs import load_proteins, Z_MOVER
    rows = [x for rs, _ in load_proteins().values() for x in rs]
    y = np.array([abs(x["z"]) > Z_MOVER for x in rows])
    az = np.clip(np.abs([x["z"] for x in rows]), 0, 6)
    g = np.array([f"{x['protein']}:{x['r']}" for x in rows])
    prot = np.array([x["protein"] for x in rows]); fam = np.array([x["family"] for x in rows])
    folds = dm.site_folds(fam, dm.FAMILY_FOLDS, 0)
    F = dm.featurize(rows)
    report = dict(n=len(y), movers=int(y.sum()), variants={}, learning_curve={})
    print(f"n={len(y)} movers={int(y.sum())}; leave-family-out, paired against baseline logreg")
    for name in FEATURES:
        X = dm.matrix(F, dm.FEATURE_SETS[name])
        V = {"baseline logreg": oof(logreg, X, y, folds)}
        V["nested logreg"], chosen = nested_logreg(X, y, g, folds)
        V["ridge on |z|"] = oof(ridge, X, az, folds, clf=False)
        V["hgb classifier"] = oof(hgb, X, y, folds)
        V["hgb regressor on |z|"] = oof(lambda: hgb(True), X, az, folds, clf=False)
        V["ensemble logreg+hgb"] = rankdata(V["baseline logreg"]) + rankdata(V["hgb classifier"])
        res = {}
        for k, p in V.items():
            a, lo, hi, _ = cluster_auc_ci(y, p, g, n_boot=1000)
            d = paired_delta_auc(y, V["baseline logreg"], p, g, n_boot=1000)
            res[k] = dict(auc=a, lo=lo, hi=hi, within=within_auc(y, p, prot),
                          d_vs_baseline=d[0], d_lo=d[1], d_hi=d[2])
            print(f"  [{name}] {k:22s} AUC={a:.3f} [{lo:.2f},{hi:.2f}] within={res[k]['within']:.3f} "
                  f"vs baseline {d[0]:+.3f} [{d[1]:+.3f},{d[2]:+.3f}]")
        res["nested_C_chosen"] = chosen
        report["variants"][name] = res
    X = dm.matrix(F, dm.FEATURE_SETS["site"])
    report["learning_curve"]["site"] = lc = learning_curve(X, y, fam, prot, folds)
    for r in lc:
        print(f"  learning curve (site): {r['fraction']:.0%} of training families -> "
              f"AUC {r['auc']:.3f} +- {r['auc_sd']:.3f}, within-protein {r['within']:.3f}")
    os.makedirs("results", exist_ok=True)
    with open("results/model_variants.json", "w") as fh:
        json.dump(report, fh, indent=1)
    print("wrote results/model_variants.json")


if __name__ == "__main__":
    main()
