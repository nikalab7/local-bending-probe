"""
proteinX GATE 3 -- the 3D ingredient.

Does adding the mutated residue's 3D CONTACT ENVIRONMENT (the AA identities
packed around it in space, from the known scaffold) rescue mover prediction
that local sequence couldn't (core AUC 0.52)?

Clean A/B: same culled training set, same T4L validation, two feature sets:
   local      = entry geometry(4) + outcome one-hot(100)
   local+3D   = local + contact-AA composition(20) + n_contacts(1)

Note on differencing: the environment is identical for WT and mutant, so it
cancels except through its INTERACTION with the central residue change -- which
a gradient-boosted model can represent. So this tests whether environment-
conditioned substitution effects are learnable, i.e. structure-based ddG-style
signal.
"""
import os
import numpy as np
from collections import Counter, defaultdict
from bending_metric import bending_angle, is_continuous
from feasibility_t4l import parse_ca, window_bend, PDB_DIR
from loop_gate import parse_ss, ss_of, pick_t4l, longest, is_conservative
from gate2_model_feasibility import entry_features, outcome_onehot, AAIDX, TRAIN_DIR

CUTOFF = 8.0


def chain_geom(res):
    nums = np.array(sorted(res))
    coords = np.array([res[n][1] for n in nums])
    aai = np.array([AAIDX.get(res[n][0], -1) for n in nums])
    pos = {int(n): i for i, n in enumerate(nums)}
    D = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1))
    return dict(nums=nums, aai=aai, pos=pos, D=D)


def env_at(g, r):
    """contact-AA composition(20) + n_contacts, for non-local neighbors (<8A, |dj|>4)."""
    if r not in g["pos"]:
        return None
    k = g["pos"][r]
    mask = (np.abs(g["nums"] - r) > 4) & (g["D"][k] < CUTOFF)
    valid = mask & (g["aai"] >= 0)
    cnt = np.bincount(g["aai"][valid], minlength=20)[:20].astype(float)
    return list(cnt) + [float(mask.sum())]


def build_training(exclude):
    files = [f[:-4] for f in os.listdir(TRAIN_DIR) if f.endswith(".pdb")]
    Xl, Xe, y, grp = [], [], [], []
    for gi, pid in enumerate(files):
        if pid.upper() in exclude:
            continue
        p = os.path.join(TRAIN_DIR, f"{pid}.pdb")
        try:
            ch, res = longest(parse_ca(p))
        except Exception:
            ch = None
        if ch is None:
            continue
        g = chain_geom(res)
        lo, hi = int(g["nums"].min()), int(g["nums"].max())
        for s in range(lo + 3, hi - 4):
            r = s + 2
            ef = entry_features(res, s)
            win = [s + k for k in range(5)]
            if ef is None or any(w not in res for w in win):
                continue
            pts = np.array([res[w][1] for w in win])
            if not is_continuous(pts):
                continue
            e = env_at(g, r)
            if e is None:
                continue
            seq5 = "".join(res[w][0] for w in win)
            Xl.append(ef + list(outcome_onehot(seq5))); Xe.append(e)
            y.append(bending_angle(pts)); grp.append(gi)
    return (np.array(Xl), np.array(Xe), np.array(y), np.array(grp))


def build_t4l_val():
    files = [f[:-4] for f in os.listdir(PDB_DIR) if f.endswith(".pdb")]
    S = {}
    for pid in files:
        p = os.path.join(PDB_DIR, f"{pid}.pdb")
        try:
            ch, res = pick_t4l(parse_ca(p))
        except Exception:
            ch = None
        if ch is None:
            continue
        helix, sheet, has = parse_ss(p)
        if has:
            S[pid] = dict(chid=ch, res=res, helix=helix, sheet=sheet)
    counts = defaultdict(Counter)
    for d in S.values():
        for rs, (aa, _) in d["res"].items():
            counts[rs][aa] += 1
    cons = {rs: c.most_common(1)[0][0] for rs, c in counts.items()
            if sum(c.values()) >= 10}
    wt_ids, mut1 = [], []
    for pid, d in S.items():
        m = tuple(sorted((rs, cons[rs], aa) for rs, (aa, _) in d["res"].items()
                         if rs in cons and aa != cons[rs]))
        if not m:
            wt_ids.append(pid)
        elif len(m) == 1:
            mut1.append((pid, m[0]))
    wt_med, wt_sig = {}, {}
    for s in range(min(cons), max(cons) - 4):
        vals = np.array([v for v in (window_bend(S[p]["res"], s) for p in wt_ids)
                         if v is not None])
        if len(vals) >= 3:
            wt_med[s] = float(np.median(vals))
            wt_sig[s] = float(1.4826 * np.median(np.abs(vals - np.median(vals))))
    pooled = float(np.median([v for v in wt_sig.values() if v > 0]))
    scaffold = max(wt_ids, key=lambda p: len(S[p]["res"]))
    gsc = chain_geom(S[scaffold]["res"])
    rows = []
    for pid, (r, wtaa, mutaa) in mut1:
        s = r - 2
        if s not in wt_med:
            continue
        b = window_bend(S[pid]["res"], s)
        ef = entry_features(S[scaffold]["res"], s)
        win = [s + k for k in range(5)]
        env = env_at(gsc, r)
        if b is None or ef is None or env is None or any(w not in cons for w in win):
            continue
        delta = b - wt_med[s]
        if is_conservative(wtaa, mutaa) and abs(delta) > 10:        # QC
            continue
        ss = ss_of(S[pid]["helix"], S[pid]["sheet"], S[pid]["chid"], r)
        wt_seq = "".join(cons[w] for w in win)
        rows.append(dict(r=r, delta=delta, mover=abs(delta) > 2 * (wt_sig.get(s, pooled) or pooled),
                         ef=ef, env=env, ss=ss,
                         wt_oh=list(outcome_onehot(wt_seq)),
                         mut_oh=list(outcome_onehot(wt_seq[:2] + mutaa + wt_seq[3:]))))
    return rows


def auc_ci(mov, score):
    from sklearn.metrics import roc_auc_score
    if not (0 < mov.sum() < len(mov)):
        return float("nan"), float("nan"), float("nan")
    auc = roc_auc_score(mov, score)
    rng = np.random.default_rng(0); idx = np.arange(len(mov)); b = []
    for _ in range(3000):
        bi = rng.choice(idx, len(idx), replace=True)
        if 0 < mov[bi].sum() < len(bi):
            b.append(roc_auc_score(mov[bi], score[bi]))
    lo, hi = np.percentile(b, [5, 95])
    return auc, lo, hi


def main():
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.metrics import mean_squared_error, average_precision_score

    t4l = set(f[:-4].upper() for f in os.listdir(PDB_DIR) if f.endswith(".pdb"))
    print("building training set (local + 3D env) ...")
    Xl, Xe, y, grp = build_training(t4l)
    print(f"  {len(y)} windows from {len(set(grp))} chains")
    X3 = np.hstack([Xl, Xe])

    gids = np.array(sorted(set(grp))); np.random.default_rng(0).shuffle(gids)
    te = np.array([g in set(gids[:max(1, len(gids)//5)]) for g in grp]); tr = ~te
    def fit(X):
        m = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.08,
                                          max_depth=6, random_state=0)
        m.fit(X[tr], y[tr]); return m
    Ml, M3 = fit(Xl), fit(X3)
    rmse_l = np.sqrt(mean_squared_error(y[te], Ml.predict(Xl[te])))
    rmse_3 = np.sqrt(mean_squared_error(y[te], M3.predict(X3[te])))
    print(f"  absolute-bending RMSE: local={rmse_l:.1f}  local+3D={rmse_3:.1f} deg "
          f"(3D cuts {100*(1-rmse_3/rmse_l):.0f}% of remaining error)")

    rows = build_t4l_val()
    mov = np.array([x["mover"] for x in rows]); n = len(rows); k = int(mov.sum())
    EFL = np.array([x["ef"] + x["wt_oh"] for x in rows])
    EFLm = np.array([x["ef"] + x["mut_oh"] for x in rows])
    ENV = np.array([x["env"] for x in rows])
    dl = np.abs(Ml.predict(EFLm) - Ml.predict(EFL))
    d3 = np.abs(M3.predict(np.hstack([EFLm, ENV])) - M3.predict(np.hstack([EFL, ENV])))

    print("\n" + "=" * 66 + "\n  GATE 3 -- 3D CONTACT ENVIRONMENT vs LOCAL-ONLY\n" + "=" * 66)
    print(f"T4L validation: n={n} single-mutant pairs, movers={k} ({k/n:.0%})")
    for lab, d in [("local-only ", dl), ("local + 3D ", d3)]:
        auc, lo, hi = auc_ci(mov, d)
        ap = average_precision_score(mov, d)
        print(f"  {lab}: AUC={auc:.3f}  90%CI[{lo:.2f},{hi:.2f}]  AP={ap:.3f} (base {k/n:.2f})")

    # core-only (helix/sheet) -- where packing dominates and local died
    core = np.array([x["ss"] in "HE" for x in rows])
    if core.sum() >= 20 and 0 < mov[core].sum() < core.sum():
        al, _, _ = auc_ci(mov[core], dl[core])
        a3, _, _ = auc_ci(mov[core], d3[core])
        print(f"  core subset (n={int(core.sum())}, movers={int(mov[core].sum())}): "
              f"local AUC={al:.3f} -> local+3D AUC={a3:.3f}")

    aL, loL, hiL = auc_ci(mov, dl)
    a3, lo3, hi3 = auc_ci(mov, d3)
    print("\n" + "-" * 66)
    if lo3 > hiL and a3 >= 0.62:
        v = ("3D HELPS -- environment-conditioned signal is real; this is the lever. "
             "But note: you've now entered structure-based ddG territory (FoldX/Rosetta/"
             "ThermoMPNN), not novel local-sequence mining.")
    elif a3 > aL + 0.04:
        v = (f"3D helps WEAKLY (AUC {aL:.2f}->{a3:.2f}); real but modest, CIs overlap. "
             "The environment matters but the local-bending-delta signal stays weak.")
    else:
        v = (f"3D does NOT rescue it (AUC {aL:.2f}->{a3:.2f}). Even the tertiary environment "
             "doesn't make the local-bending DELTA predictable -- the deepest negative: the "
             "mutation's effect on local bending isn't a learnable function of its neighbors either.")
    print("VERDICT:", v)

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve
    plt.figure(figsize=(5.4, 5))
    for lab, d, c in [("local-only", dl, "#888888"), ("local + 3D", d3, "#4C78A8")]:
        fpr, tpr, _ = roc_curve(mov, d); a, _, _ = auc_ci(mov, d)
        plt.plot(fpr, tpr, color=c, label=f"{lab} (AUC {a:.2f})")
    plt.plot([0, 1], [0, 1], "k--", lw=.8, label="chance")
    plt.xlabel("FPR"); plt.ylabel("TPR"); plt.title("Does the 3D environment rescue it?")
    plt.legend(); plt.tight_layout(); plt.savefig("gate3_3d.png", dpi=130)
    print("plot -> gate3_3d.png")


if __name__ == "__main__":
    main()
