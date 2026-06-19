"""
proteinX POWERED LOOP GATE -- replicate the T4L loop signal (AUC 0.70, CI
[0.49,0.89]) on INDEPENDENT loop-rich proteins with enough movers to tighten
the interval.

Same corrected gate as loop_gate.py: train absolute bending on LOOP windows of
the culled non-T4L set (now also excluding these validation proteins -> no
leakage), difference the model on held-out loop mutation pairs, score full-set
mover retrieval (AUC/AP/p@k) with a bootstrap CI.

Validation proteins (loop-rich, single-mutant-rich, NOT T4L):
  SNase (P00644), barnase (P00648), human lysozyme (P61626), RNase A (P61823).
"""
import os, json, urllib.request, concurrent.futures as cf
import numpy as np
from collections import Counter, defaultdict
from bending_metric import bending_angle, is_continuous
from feasibility_t4l import parse_ca, window_bend
from loop_gate import parse_ss, ss_of, is_conservative
from gate2_model_feasibility import entry_features, outcome_onehot, TRAIN_DIR
from mover_composition import has_nonlocal

VAL_DIR = "val_pdb"; os.makedirs(VAL_DIR, exist_ok=True)
PROTEINS = [("P00644", "SNase", 149), ("P00648", "barnase", 110),
            ("P61626", "human_lysozyme", 130), ("P61823", "RNaseA", 124)]
MAX_PER = 400


def fetch_ids(uniprot):
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": uniprot}},
         "return_type": "entry",
         "request_options": {"return_all_hits": True}}
    req = urllib.request.Request(
        "https://search.rcsb.org/rcsbsearch/v2/query",
        data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    try:
        res = json.load(urllib.request.urlopen(req, timeout=50))
        return [r["identifier"] for r in res["result_set"]]
    except Exception:
        return []


def fetch_pdb(pid):
    path = os.path.join(VAL_DIR, f"{pid}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return pid
    try:
        data = urllib.request.urlopen(
            f"https://files.rcsb.org/download/{pid}.pdb", timeout=30).read()
        open(path, "wb").write(data)
        return pid
    except Exception:
        return None


def pick_chain(chains, hint, tol=30):
    best, sc = (None, None), 1e9
    for ch, res in chains.items():
        n = len(res)
        if abs(n - hint) <= tol and max(res) < hint + tol + 50 and n < 250 \
                and abs(n - hint) < sc:
            best, sc = (ch, res), abs(n - hint)
    return best


def build_protein(uniprot, name, hint, val_ids):
    ids = fetch_ids(uniprot)[:MAX_PER]
    with cf.ThreadPoolExecutor(max_workers=12) as ex:
        got = [p for p in ex.map(fetch_pdb, ids) if p]
    val_ids.update(p.upper() for p in got)

    S = {}
    for pid in got:
        p = os.path.join(VAL_DIR, f"{pid}.pdb")
        try:
            ch, res = pick_chain(parse_ca(p), hint)
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
            if sum(c.values()) >= 5}

    wt_ids, mut1 = [], []
    for pid, d in S.items():
        m = tuple(sorted((rs, cons[rs], aa) for rs, (aa, _) in d["res"].items()
                         if rs in cons and aa != cons[rs]))
        if not m:
            wt_ids.append(pid)
        elif len(m) == 1:
            mut1.append((pid, m[0]))

    diag = dict(name=name, structs=len(S), wt=len(wt_ids), single=len(mut1))
    if len(wt_ids) < 3:
        diag["note"] = "skip: <3 WT crystals for floor"
        return diag, []

    wt_med, wt_sig = {}, {}
    for s in range(min(cons), max(cons) - 4):
        vals = np.array([v for v in (window_bend(S[p]["res"], s) for p in wt_ids)
                         if v is not None])
        if len(vals) >= 3:
            wt_med[s] = float(np.median(vals))
            wt_sig[s] = float(1.4826 * np.median(np.abs(vals - np.median(vals))))
    if not wt_sig:
        diag["note"] = "skip: no per-window floor"
        return diag, []
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
        if ss_of(S[pid]["helix"], S[pid]["sheet"], S[pid]["chid"], r) != "L":
            continue
        delta = b - wt_med[s]
        if is_conservative(wtaa, mutaa) and abs(delta) > 10:      # QC
            continue
        floor = wt_sig.get(s, pooled) or pooled
        wt_seq = "".join(cons[w] for w in win)
        rows.append(dict(prot=name, pid=pid, r=r, delta=delta,
                         mover=abs(delta) > 2 * floor, ef=ef, wt_seq=wt_seq,
                         mut_seq=wt_seq[:2] + mutaa + wt_seq[3:]))
    diag["loop_pairs"] = len(rows)
    diag["loop_movers"] = sum(x["mover"] for x in rows)
    return diag, rows


def train_loop_model(exclude_ids):
    from sklearn.ensemble import HistGradientBoostingRegressor
    files = [f[:-4] for f in os.listdir(TRAIN_DIR) if f.endswith(".pdb")]
    X, y = [], []
    for pid in files:
        if pid.upper() in exclude_ids:
            continue
        p = os.path.join(TRAIN_DIR, f"{pid}.pdb")
        try:
            chains = parse_ca(p)
        except Exception:
            continue
        ch, res = max(((c, r) for c, r in chains.items()),
                      key=lambda cr: len(cr[1]), default=(None, None))
        if res is None or len(res) < 40 or max(res) >= 5000:
            continue
        helix, sheet, has = parse_ss(p)
        if not has:
            continue
        lo, hi = min(res), max(res)
        for i in range(lo + 3, hi - 4):
            if ss_of(helix, sheet, ch, i + 2) != "L":
                continue
            ef = entry_features(res, i)
            win = [i + k for k in range(5)]
            if ef is None or any(w not in res for w in win):
                continue
            pts = np.array([res[w][1] for w in win])
            if not is_continuous(pts):
                continue
            seq5 = "".join(res[w][0] for w in win)
            X.append(ef + list(outcome_onehot(seq5)))
            y.append(bending_angle(pts))
    m = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.08,
                                      max_depth=6, random_state=0)
    m.fit(np.array(X), np.array(y))
    return m, len(y)


def score(rows, model, label):
    from sklearn.metrics import roc_auc_score, average_precision_score
    EF = np.array([x["ef"] for x in rows])
    Xwt = np.hstack([EF, np.array([outcome_onehot(x["wt_seq"]) for x in rows])])
    Xmt = np.hstack([EF, np.array([outcome_onehot(x["mut_seq"]) for x in rows])])
    dp = np.abs(model.predict(Xmt) - model.predict(Xwt))
    mov = np.array([x["mover"] for x in rows]); n = len(rows); k = int(mov.sum())
    if not (0 < k < n):
        print(f"  {label}: n={n} movers={k} -> cannot score"); return
    auc = roc_auc_score(mov, dp); ap = average_precision_score(mov, dp)
    patk = mov[np.argsort(-dp)[:k]].mean()
    rng = np.random.default_rng(0); idx = np.arange(n); b = []
    for _ in range(3000):
        bi = rng.choice(idx, n, replace=True)
        if 0 < mov[bi].sum() < len(bi):
            b.append(roc_auc_score(mov[bi], dp[bi]))
    lo, hi = np.percentile(b, [5, 95])
    print(f"  {label}: n={n} movers={k} (base {k/n:.0%})  "
          f"AUC={auc:.3f} 90%CI[{lo:.2f},{hi:.2f}]  AP={ap:.3f}  p@{k}={patk:.3f}")
    return auc, lo, hi, n, k


def main():
    val_ids = set()
    pooled, diags = [], []
    print("building independent loop-pair validation set ...")
    for up, name, hint in PROTEINS:
        d, rows = build_protein(up, name, hint, val_ids)
        diags.append(d); pooled += rows
        print(f"  {name:16s} structs={d.get('structs',0):3d} WT={d.get('wt',0):3d} "
              f"single={d.get('single',0):3d} loop_pairs={d.get('loop_pairs',0):3d} "
              f"loop_movers={d.get('loop_movers',0):3d} {d.get('note','')}")

    print(f"\ntraining loop model on culled set (excluding {len(val_ids)} val ids + T4L) ...")
    model, ntrain = train_loop_model(val_ids)
    print(f"  trained on {ntrain} loop windows")

    print("\n" + "=" * 70 + "\n  POWERED LOOP RETRIEVAL\n" + "=" * 70)
    res = score(pooled, model, "POOLED (independent, non-T4L)")

    print("\n" + "-" * 70)
    if res is None:
        print("VERDICT: INCONCLUSIVE -- not enough pooled loop movers; add more proteins.")
    else:
        auc, lo, hi, n, k = res
        if k < 30:
            print(f"VERDICT: STILL UNDERPOWERED -- only {k} pooled loop movers; add proteins.")
        elif lo > 0.55 and auc >= 0.62:
            print(f"VERDICT: CONFIRMED -- T4L's loop signal REPLICATES on independent data "
                  f"(AUC {auc:.2f}, CI excludes chance). Local thesis lives in loops; a "
                  f"loop-restricted miner is now justified.")
        elif hi < 0.6:
            print(f"VERDICT: COLLAPSED -- loop AUC {auc:.2f}, CI [{lo:.2f},{hi:.2f}] near chance. "
                  f"T4L's 0.70 was a small-sample fluke. Local premise dead in loops too.")
        else:
            print(f"VERDICT: WEAK/AMBIGUOUS -- AUC {auc:.2f}, CI [{lo:.2f},{hi:.2f}]. Real but "
                  f"marginal; at best a weak ranker. Decide if weak-and-interpretable is worth it.")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve
    if res:
        EF = np.array([x["ef"] for x in pooled])
        dp = np.abs(model.predict(np.hstack([EF, np.array([outcome_onehot(x["mut_seq"]) for x in pooled])]))
                    - model.predict(np.hstack([EF, np.array([outcome_onehot(x["wt_seq"]) for x in pooled])])))
        mov = np.array([x["mover"] for x in pooled])
        fpr, tpr, _ = roc_curve(mov, dp)
        plt.figure(figsize=(5.2, 5))
        plt.plot(fpr, tpr, color="#4C78A8", label=f"pooled AUC {res[0]:.2f} (n={res[3]})")
        plt.plot([0, 1], [0, 1], "k--", lw=.8, label="chance")
        plt.xlabel("FPR"); plt.ylabel("TPR")
        plt.title("powered loop mover retrieval"); plt.legend()
        plt.tight_layout(); plt.savefig("powered_loop_gate.png", dpi=130)
        print("plot -> powered_loop_gate.png")


if __name__ == "__main__":
    main()
