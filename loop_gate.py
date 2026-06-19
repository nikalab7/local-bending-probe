"""
proteinX LOOP GATE -- does local-bending signal survive in loops/turns?

Cores failed (gate2 AUC 0.52) because tertiary packing dominates there. Loops
aren't packing-locked, so local sequence should carry more signal -- the last
place the local thesis can be alive. Four checks, all read-only on cached data:

  1. Split T4L movers/non-movers by secondary structure (PDB HELIX/SHEET records;
     else loop). Are movers enriched in loops?
  2. Loop-only corrected gate2: train absolute bending on LOOP windows from the
     culled non-T4L set, difference on held-out T4L loop pairs, score full-set
     mover retrieval (AUC/AP/p@k). Compare to the global 0.52.
  3. Dataset viability: after the non-local-contact filter, how many clean LOOP
     windows actually survive to train on?
  4. QC: flag conservative-core-substitution outliers (V111M/V111I type) so they
     don't poison loop labels.

GO  : loop AUC comfortably > 0.52 AND enough clean loop windows -> re-scope here.
NO-GO: still ~0.5 OR too few windows -> local premise dead; re-aim target.
"""
import os
import numpy as np
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
from bending_metric import bending_angle, is_continuous
from feasibility_t4l import parse_ca, window_bend, PDB_DIR
from gate2_model_feasibility import entry_features, outcome_onehot, TRAIN_DIR
from mover_composition import has_nonlocal, window_nonlocal

CONS_GROUPS = [set("AVLIM"), set("FYW"), set("ST"), set("DE"),
               set("NQ"), set("KR")]


def is_conservative(a, b):
    return any(a in g and b in g for g in CONS_GROUPS)


def parse_ss(path):
    """(helix_set, sheet_set, has_records) of (chain,resnum) from header."""
    helix, sheet, has = set(), set(), False
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM", "MODEL")):
                break
            if line.startswith("HELIX"):
                try:
                    ch, a, b = line[19], int(line[21:25]), int(line[33:37])
                    has = True
                    helix.update((ch, r) for r in range(a, b + 1))
                except ValueError:
                    pass
            elif line.startswith("SHEET"):
                try:
                    ch, a, b = line[21], int(line[22:26]), int(line[33:37])
                    has = True
                    sheet.update((ch, r) for r in range(a, b + 1))
                except ValueError:
                    pass
    return helix, sheet, has


def ss_of(helix, sheet, ch, r):
    if (ch, r) in helix:
        return "H"
    if (ch, r) in sheet:
        return "E"
    return "L"


def pick_t4l(chains):
    best, sc = (None, None), 1e9
    for ch, res in chains.items():
        n = len(res)
        if 150 <= n <= 175 and max(res) <= 200 and abs(n - 164) < sc:
            best, sc = (ch, res), abs(n - 164)
    return best


def longest(chains):
    best, n = (None, None), 0
    for ch, res in chains.items():
        if len(res) > n and max(res) < 5000:
            best, n = (ch, res), len(res)
    return best if n >= 40 else (None, None)


# --------------------- T4L: rebuild movers w/ SS + QC -----------------------
def build_t4l():
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
        if not has:
            continue
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

    rows = []
    for pid, (r, wtaa, mutaa) in mut1:
        s = r - 2
        if s not in wt_med:
            continue
        b = window_bend(S[pid]["res"], s)
        ef = entry_features(S[scaffold]["res"], s)
        win = [s + k for k in range(5)]
        if b is None or ef is None or any(w not in cons for w in win):
            continue
        delta = b - wt_med[s]
        floor = wt_sig.get(s, pooled) or pooled
        ss = ss_of(S[pid]["helix"], S[pid]["sheet"], S[pid]["chid"], r)
        nl = has_nonlocal(S[pid]["res"], r)
        # conservative swap + implausibly large bend = artifact, contact or not
        qc_bad = is_conservative(wtaa, mutaa) and abs(delta) > 10
        wt_seq = "".join(cons[w] for w in win)
        rows.append(dict(pid=pid, r=r, wt=wtaa, mut=mutaa, delta=delta,
                         floor=floor, mover=abs(delta) > 2 * floor, ss=ss,
                         qc_bad=qc_bad, ef=ef, wt_seq=wt_seq,
                         mut_seq=wt_seq[:2] + mutaa + wt_seq[3:]))
    return rows


# --------------------- cull set: loop windows + viability -------------------
def build_training(require_clean=False, max_struct=None):
    files = [f[:-4] for f in os.listdir(TRAIN_DIR) if f.endswith(".pdb")]
    if max_struct:
        files = files[:max_struct]
    X, y, grp = [], [], []
    n_loop_total = n_loop_clean = n_skip_no_ss = 0
    for gi, pid in enumerate(files):
        p = os.path.join(TRAIN_DIR, f"{pid}.pdb")
        try:
            ch, res = longest(parse_ca(p))
        except Exception:
            ch = None
        if ch is None:
            continue
        helix, sheet, has = parse_ss(p)
        if not has:
            n_skip_no_ss += 1
            continue
        lo, hi = min(res), max(res)
        for i in range(lo + 3, hi - 4):
            r = i + 2                                   # central residue
            if ss_of(helix, sheet, ch, r) != "L":
                continue
            ef = entry_features(res, i)
            win = [i + k for k in range(5)]
            if ef is None or any(w not in res for w in win):
                continue
            pts = np.array([res[w][1] for w in win])
            if not is_continuous(pts):
                continue
            n_loop_total += 1
            clean = not window_nonlocal(res, r)
            if clean:
                n_loop_clean += 1
            if require_clean and not clean:
                continue
            seq5 = "".join(res[w][0] for w in win)
            X.append(ef + list(outcome_onehot(seq5)))
            y.append(bending_angle(pts)); grp.append(gi)
    return (np.array(X), np.array(y), np.array(grp),
            n_loop_total, n_loop_clean, n_skip_no_ss)


def main():
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.metrics import (roc_auc_score, average_precision_score,
                                 mean_squared_error, roc_curve)

    rows = build_t4l()
    flagged = [x for x in rows if x["qc_bad"]]
    clean = [x for x in rows if not x["qc_bad"]]
    print("=" * 68 + "\n  LOOP GATE\n" + "=" * 68)
    print(f"[QC] conservative-core outliers flagged & removed: {len(flagged)}")
    for x in flagged:
        print(f"      {x['pid']} {x['wt']}{x['r']}{x['mut']} delta={x['delta']:+.1f}  ss={x['ss']}")

    # ---- 1. SS distribution of movers ----
    print("\n[1] T4L mover distribution by secondary structure (central residue):")
    by = {c: [x for x in clean if x["ss"] == c] for c in "HEL"}
    for c, name in [("H", "helix"), ("E", "sheet"), ("L", "loop/turn")]:
        g = by[c]
        mv = sum(x["mover"] for x in g)
        rate = mv / len(g) if g else float("nan")
        print(f"    {name:10s}: {len(g):3d} windows, {mv:3d} movers ({rate:5.0%})")
    loop = by["L"]; nonloop = by["H"] + by["E"]
    a = sum(x["mover"] for x in loop); b = len(loop) - a
    c_ = sum(x["mover"] for x in nonloop); d = len(nonloop) - c_
    orr, pf = fisher_exact([[a, b], [c_, d]], alternative="greater")
    print(f"    loop-mover enrichment vs non-loop: OR={orr:.2f}  Fisher p={pf:.3g}")

    # ---- 3. dataset viability ----
    Xtr, ytr, grp, n_lt, n_lc, n_noss = build_training(require_clean=False)
    print(f"\n[3] loop-window viability (culled non-T4L set):")
    print(f"    structures skipped (no SS records): {n_noss}")
    print(f"    loop-central windows total        : {n_lt}")
    print(f"    clean loop windows (no non-local) : {n_lc} ({100*n_lc/max(n_lt,1):.0f}%)")

    # ---- 2. loop-only corrected gate ----
    gids = np.array(sorted(set(grp))); np.random.default_rng(0).shuffle(gids)
    test_g = set(gids[: max(1, len(gids) // 5)])
    te = np.array([g in test_g for g in grp]); tr = ~te
    model = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.08,
                                          max_depth=6, random_state=0)
    model.fit(Xtr[tr], ytr[tr])
    rmse = float(np.sqrt(mean_squared_error(ytr[te], model.predict(Xtr[te]))))
    std = float(np.std(ytr[te]))
    print(f"\n[2] loop-only model: absolute bending std={std:.1f}  RMSE={rmse:.1f} deg "
          f"(skill {100*(1-rmse/std):.0f}%)")

    val = [x for x in clean if x["ss"] == "L"]
    if len(val) < 8:
        print(f"    held-out T4L loop pairs = {len(val)}  -> TOO FEW to score.")
        auc = ap = patk = float("nan"); movers_n = sum(x["mover"] for x in val)
    else:
        EF = np.array([x["ef"] for x in val])
        Xwt = np.hstack([EF, np.array([outcome_onehot(x["wt_seq"]) for x in val])])
        Xmt = np.hstack([EF, np.array([outcome_onehot(x["mut_seq"]) for x in val])])
        dpred = np.abs(model.predict(Xmt) - model.predict(Xwt))
        mov = np.array([x["mover"] for x in val]); movers_n = int(mov.sum())
        base = mov.mean()
        auc = roc_auc_score(mov, dpred) if 0 < movers_n < len(val) else float("nan")
        ap = average_precision_score(mov, dpred)
        k = movers_n
        patk = mov[np.argsort(-dpred)[:k]].mean() if k else float("nan")
        # bootstrap CI on AUC -- n is tiny, so quantify how shaky 0.70 is
        rng = np.random.default_rng(0); idx = np.arange(len(val)); boots = []
        for _ in range(3000):
            bi = rng.choice(idx, len(idx), replace=True)
            if 0 < mov[bi].sum() < len(bi):
                boots.append(roc_auc_score(mov[bi], dpred[bi]))
        lo, hi = np.percentile(boots, [5, 95])
        print(f"    held-out T4L loop pairs = {len(val)} (movers {movers_n}, base {base:.0%})")
        print(f"    loop-only mover retrieval:  AUC={auc:.3f}  90% CI [{lo:.2f}, {hi:.2f}]  (core was 0.52)")
        print(f"                                AP ={ap:.3f}  (base {base:.3f})")
        print(f"                                p@{k}={patk:.3f}")

    # ---- verdict: separate SIGNAL, POWER, and TRAIN-viability ----
    auc_ok = (not np.isnan(auc)) and auc >= 0.60
    train_ok = n_lc >= 2000
    val_powered = movers_n >= 20
    print("\n" + "-" * 68)
    print(f"signal: loop AUC {auc:.3f} (vs core 0.52), enrichment OR {orr:.2f} (p={pf:.3f})")
    print(f"power : {movers_n} loop movers in T4L  -> {'adequate' if val_powered else 'TOO FEW to certify'}")
    print(f"train : {n_lc} clean loop windows       -> {'trainable' if train_ok else 'too few'}")
    if auc_ok and train_ok and val_powered:
        v = "GO -- local thesis SURVIVES in loops. Re-scope around flexible regions."
    elif auc_ok and train_ok and not val_powered:
        v = ("UNDERPOWERED-PROMISING -- signal is consistently POSITIVE in loops (AUC + "
             "enrichment), unlike cores, but T4L has too few loop mutations to certify it. "
             "NOT dead, NOT build-ready: the decisive test needs a LOOP-RICH mutational "
             "dataset, not T4L. Do NOT build the miner on this alone.")
    elif not auc_ok:
        v = (f"DEAD -- loop AUC {auc:.3f} not above core's 0.52. Local premise dead even in "
             "loops; re-aim (flexibility/B-factor or loop mutation-ranking).")
    else:
        v = f"NO-GO -- insufficient clean loop training data (clean={n_lc})."
    print("VERDICT:", v)

    # plot
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    names = ["helix", "sheet", "loop"]
    rates = [np.mean([x["mover"] for x in by[c]]) if by[c] else 0 for c in "HEL"]
    ax[0].bar(names, rates, color=["#54A24B", "#F58518", "#4C78A8"])
    ax[0].set_ylabel("mover fraction"); ax[0].set_title("T4L mover rate by SS")
    for i, (c, rt) in enumerate(zip("HEL", rates)):
        ax[0].text(i, rt + 0.005, f"{rt:.0%}\n(n={len(by[c])})", ha="center")
    if not np.isnan(auc) and len(val) >= 8 and 0 < movers_n < len(val):
        mov = np.array([x["mover"] for x in val])
        EF = np.array([x["ef"] for x in val])
        dp = np.abs(model.predict(np.hstack([EF, np.array([outcome_onehot(x["mut_seq"]) for x in val])]))
                    - model.predict(np.hstack([EF, np.array([outcome_onehot(x["wt_seq"]) for x in val])])))
        fpr, tpr, _ = roc_curve(mov, dp)
        ax[1].plot(fpr, tpr, color="#4C78A8", label=f"loop AUC {auc:.2f}")
    ax[1].plot([0, 1], [0, 1], "k--", lw=.8, label="chance / core 0.52")
    ax[1].set_xlabel("FPR"); ax[1].set_ylabel("TPR")
    ax[1].set_title("loop-only mover retrieval"); ax[1].legend()
    fig.tight_layout(); fig.savefig("loop_gate.png", dpi=130)
    print("plot -> loop_gate.png")


if __name__ == "__main__":
    main()
