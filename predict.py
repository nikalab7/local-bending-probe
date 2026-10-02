"""
proteinX -- score point mutations for backbone movement (the usable end product).

  python predict.py train
      Fits the mover model on every clean label (pairs.load_proteins), with
      confidence weights (delta_model.label_weights), reports its leave-family-
      out accuracy, calibrates probabilities on those out-of-fold scores, and
      saves results/mover_model.pkl.

  python predict.py score --pdb 2LZM --chain A --mut L99A,T26E,G28A
      Scores mutations on one WT structure (PDB ID, downloaded, or a local
      .pdb/.pdb.gz path). Output per mutation: calibrated probability that it
      moves the 5-residue window around it beyond crystal noise, its
      percentile among the training mutations, and the site drivers.

Model: logistic regression on SITE + SUBST features (C-alpha site context and
the substitution), the best pre-specified model in RESULTS.md v4/v5. It needs
only C-alpha coordinates, B-factors and HELIX/SHEET records.

Read the numbers for what they are: leave-family-out AUC ~0.62 against all
labels and ~0.70 against labels that are clearly mover / clearly not (|z| > 3
vs < 1); probabilities are calibrated to the training mix (~26% movers). Good
for ranking candidate mutations by risk of local backbone change, not for
yes/no calls on single mutations. One difference from training: WT window
bending comes from this one structure, not from the median of a crystal form.
"""
from __future__ import annotations
import argparse
import json
import os
import pickle
import urllib.request
import numpy as np

MODEL_PATH = os.path.join("results", "mover_model.pkl")
FEATURES = "site+subst"


def _xy(rows, zthr):
    import delta_model as dm
    F = dm.featurize(rows)
    X = dm.matrix(F, dm.FEATURE_SETS[FEATURES])
    z = np.array([x["z"] for x in rows])
    return X, np.abs(z) > zthr, z


def train():
    import delta_model as dm
    from pairs import load_proteins, Z_MOVER
    from sklearn.linear_model import LogisticRegression
    from stats_utils import cluster_auc_ci, within_auc
    rows = [x for rs, _ in load_proteins().values() for x in rs]
    X, y, z = _xy(rows, Z_MOVER)
    w = dm.label_weights(z, Z_MOVER)
    fam = np.array([x["family"] for x in rows])
    groups = np.array([f"{x['protein']}:{x['r']}" for x in rows])
    prot = np.array([x["protein"] for x in rows])
    oof = dm.oof_predict(X, y, dm.site_folds(fam, dm.FAMILY_FOLDS, 0), "logreg", w)
    auc = cluster_auc_ci(y, oof, groups)
    conf = (np.abs(z) > dm.CONFIDENT_HI) | (np.abs(z) < dm.CONFIDENT_LO)
    auc_c = cluster_auc_ci(y[conf], oof[conf], groups[conf])
    wit = within_auc(y, oof, prot)
    # Platt calibration on out-of-fold scores (class_weight="balanced" skews raw probabilities)
    logit = np.log(np.clip(oof, 1e-6, 1 - 1e-6) / np.clip(1 - oof, 1e-6, 1))
    cal = LogisticRegression().fit(logit[:, None], y)
    model = dm.fit_model("logreg", X, y, w)
    final_scores = model.predict_proba(X)[:, 1]
    out = dict(model=model, calibrator=cal, features=dm.FEATURE_SETS[FEATURES],
               train_scores=np.sort(final_scores), z_threshold=Z_MOVER,
               metrics=dict(n=len(y), movers=int(y.sum()), proteins=len(set(prot)),
                            families=len(set(fam)),
                            auc_family_out=auc[:3], auc_confident_family_out=auc_c[:3],
                            within_protein_auc=wit, n_confident=int(conf.sum())))
    os.makedirs("results", exist_ok=True)
    with open(MODEL_PATH, "wb") as fh:
        pickle.dump(out, fh)
    with open(MODEL_PATH.replace(".pkl", ".json"), "w") as fh:
        json.dump(dict(features=out["features"], **out["metrics"],
                       coefficients=dict(zip(out["features"],
                                             map(float, model[-1].coef_[0])))), fh, indent=1)
    m = out["metrics"]
    print(f"trained on {m['n']} mutations ({m['movers']} movers, {m['proteins']} proteins, "
          f"{m['families']} families)")
    print(f"leave-family-out AUC {auc[0]:.3f} [{auc[1]:.2f},{auc[2]:.2f}]; on confident labels "
          f"(n={m['n_confident']}) {auc_c[0]:.3f} [{auc_c[1]:.2f},{auc_c[2]:.2f}]; "
          f"within-protein {wit:.3f}")
    print(f"saved {MODEL_PATH}")


def _structure(pdb):
    if os.path.exists(pdb):
        return pdb
    path = os.path.join("results", f"_{pdb.upper()}.pdb")
    if not os.path.exists(path):
        data = urllib.request.urlopen(
            f"https://files.rcsb.org/download/{pdb.upper()}.pdb", timeout=60).read()
        with open(path, "wb") as fh:
            fh.write(data)
    return path


def parse_mutation(s):
    s = s.strip().upper()
    return s[0], int(s[1:-1]), s[-1]


def score(pdb, chain, muts):
    import delta_model as dm
    import pairs
    with open(MODEL_PATH, "rb") as fh:
        M = pickle.load(fh)
    path = _structure(pdb)
    st = pairs.parse_structure(path)
    if chain not in st["chains"]:
        raise SystemExit(f"chain {chain} not in {path}: {sorted(st['chains'])}")
    res = st["chains"][chain]
    rows, skipped = [], []
    for m in muts:
        wt, r, mut = parse_mutation(m)
        if r not in res:
            skipped.append(f"{m}: residue {r} not resolved"); continue
        if res[r][0] != wt:
            skipped.append(f"{m}: structure has {res[r][0]}{r}, not {wt}"); continue
        bend = pairs.window_bend(res, r - 2)
        if bend is None:
            skipped.append(f"{m}: 5-residue window {r-2}..{r+2} incomplete"); continue
        rows.append(dict(mutation=m, protein="query", r=r, s=r - 2, wt=wt, mut=mut,
                         scaffold_res=res, scaffold_ss=pairs.ss_of(st, chain, r),
                         wt_bend=bend))
    out = []
    if rows:
        X = dm.matrix(dm.featurize([{k: v for k, v in x.items() if k != "mutation"}
                                    for x in rows]), M["features"])
        raw = M["model"].predict_proba(X)[:, 1]
        logit = np.log(np.clip(raw, 1e-6, 1 - 1e-6) / np.clip(1 - raw, 1e-6, 1))
        prob = M["calibrator"].predict_proba(logit[:, None])[:, 1]
        pct = np.searchsorted(M["train_scores"], raw) / len(M["train_scores"]) * 100
        F = dm.featurize(rows)
        for x, p, q, f in zip(rows, prob, pct, F):
            out.append(dict(mutation=x["mutation"], p_mover=float(p), percentile=float(q),
                            ss=x["scaffold_ss"], wt_bend=round(f["wt_bend"], 1),
                            nonlocal_contacts=f["nonlocal_contacts"], burial=f["hse_up"]))
    return out, skipped


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("train")
    s = sub.add_parser("score")
    s.add_argument("--pdb", required=True, help="PDB ID or path to a WT structure")
    s.add_argument("--chain", default="A")
    s.add_argument("--mut", required=True, help="comma-separated, e.g. L99A,T26E")
    a = ap.parse_args()
    if a.cmd == "train":
        train(); return
    out, skipped = score(a.pdb, a.chain, a.mut.split(","))
    print(f"{'mutation':10s} {'P(mover)':>9s} {'pctile':>7s}  SS  WT bend  non-local  burial")
    for o in sorted(out, key=lambda o: -o["p_mover"]):
        print(f"{o['mutation']:10s} {o['p_mover']:9.2f} {o['percentile']:6.0f}%  {o['ss']:2s} "
              f"{o['wt_bend']:7.1f}  {o['nonlocal_contacts']:9.0f}  {o['burial']:6.0f}")
    for s_ in skipped:
        print("skipped:", s_)


if __name__ == "__main__":
    main()
