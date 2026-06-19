"""
proteinX mover-composition diagnostic.

Gate 2 showed a local predict-absolute-then-difference model can't retrieve the
movers. This asks WHY, and whether ANY local model could: of the observed
movers, how many are explained by causes a LOCAL model cannot or trivially can
handle --
  * Pro/Gly involvement  (textbook local effect; recovering it is not novelty)
  * non-local tertiary contact at the window (>15 seq apart, <8A Ca-Ca;
    a cause a LOCAL model cannot see by construction)

The decisive number is the RESIDUAL: movers that are NEITHER -> the only set the
v2 higher-order conditioned-combination thesis could own. Rates are compared to
NON-movers (a folded two-domain protein has tertiary contacts almost
everywhere, so absolute rates mislead; enrichment is what matters).
"""
import os
import numpy as np
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
from feasibility_t4l import parse_ca, pick_t4l_chain, window_bend, PDB_DIR


def has_nonlocal(res, q, cutoff=8.0, seqsep=15):
    if q not in res:
        return False
    cq = res[q][1]
    return any(abs(j - q) > seqsep and np.linalg.norm(cq - cj) < cutoff
               for j, (_, cj) in res.items())


def window_nonlocal(res, r):
    return any(has_nonlocal(res, q) for q in range(r - 2, r + 3) if q in res)


def main():
    files = [f[:-4] for f in os.listdir(PDB_DIR) if f.endswith(".pdb")]
    structs = {}
    for pid in files:
        try:
            ch = pick_t4l_chain(parse_ca(os.path.join(PDB_DIR, f"{pid}.pdb")))
        except Exception:
            ch = None
        if ch:
            structs[pid] = ch

    counts = defaultdict(Counter)
    for res in structs.values():
        for rs, (aa, _) in res.items():
            counts[rs][aa] += 1
    consensus = {rs: c.most_common(1)[0][0] for rs, c in counts.items()
                 if sum(c.values()) >= 10}

    wt_ids, mut1 = [], []
    for pid, res in structs.items():
        m = tuple(sorted((rs, consensus[rs], aa) for rs, (aa, _) in res.items()
                         if rs in consensus and aa != consensus[rs]))
        if not m:
            wt_ids.append(pid)
        elif len(m) == 1:
            mut1.append((pid, m[0]))

    wt_med, wt_sig = {}, {}
    for s in range(min(consensus), max(consensus) - 4):
        vals = np.array([v for v in (window_bend(structs[p], s) for p in wt_ids)
                         if v is not None])
        if len(vals) >= 3:
            wt_med[s] = float(np.median(vals))
            wt_sig[s] = float(1.4826 * np.median(np.abs(vals - np.median(vals))))
    pooled = float(np.median([v for v in wt_sig.values() if v > 0]))

    rows = []
    for pid, (r, wtaa, mutaa) in mut1:
        s = r - 2
        if s not in wt_med:
            continue
        b = window_bend(structs[pid], s)
        if b is None:
            continue
        floor = wt_sig.get(s, pooled) or pooled
        delta = b - wt_med[s]
        rows.append(dict(
            pid=pid, r=r, wt=wtaa, mut=mutaa, delta=delta,
            mover=abs(delta) > 2 * floor,
            progly=(wtaa in "PG") or (mutaa in "PG"),
            nonlocal_win=window_nonlocal(structs[pid], r)))

    mv = [x for x in rows if x["mover"]]
    nm = [x for x in rows if not x["mover"]]

    def frac(group, key):
        return np.mean([x[key] for x in group]) if group else float("nan")

    def fisher(key):
        a = sum(x[key] for x in mv); b = len(mv) - a
        c = sum(x[key] for x in nm); d = len(nm) - c
        orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        return a, c, orr, p

    print("=" * 66 + "\n  MOVER COMPOSITION DIAGNOSTIC (T4 lysozyme)\n" + "=" * 66)
    print(f"single-mutation windows: {len(rows)}   movers: {len(mv)}   "
          f"non-movers: {len(nm)}\n")
    for key, label in [("progly", "Pro/Gly created or destroyed"),
                       ("nonlocal_win", "non-local contact in window (<8A,>15 sep)")]:
        a, c, orr, p = fisher(key)
        print(f"{label}:")
        print(f"    movers     {a:3d}/{len(mv)} = {frac(mv,key):5.0%}")
        print(f"    non-movers {c:3d}/{len(nm)} = {frac(nm,key):5.0%}")
        print(f"    enrichment OR = {orr:.2f}   Fisher p(greater) = {p:.3g}\n")

    residual = [x for x in mv if not x["progly"] and not x["nonlocal_win"]]
    expl_pg = [x for x in mv if x["progly"]]
    expl_nl = [x for x in mv if x["nonlocal_win"] and not x["progly"]]
    print("-" * 66)
    print(f"MOVER BREAKDOWN ({len(mv)} total):")
    print(f"    Pro/Gly explained          : {len(expl_pg):3d} ({len(expl_pg)/len(mv):.0%})")
    print(f"    non-local (not Pro/Gly)    : {len(expl_nl):3d} ({len(expl_nl)/len(mv):.0%})")
    print(f"    RESIDUAL (neither)         : {len(residual):3d} ({len(residual)/len(mv):.0%})"
          "   <-- the only set v2's local higher-order thesis could own")
    print("-" * 66)
    if residual:
        print("residual movers (local, non-Pro/Gly, no detected non-local contact):")
        for x in sorted(residual, key=lambda z: -abs(z["delta"]))[:20]:
            print(f"    {x['pid']}  {x['wt']}{x['r']}{x['mut']}  "
                  f"delta={x['delta']:+.1f} deg")

    # interpretation
    print("\nREAD:")
    nl_or = fisher("nonlocal_win")[2]
    if len(residual) <= 0.15 * len(mv) or (nl_or > 2 and len(residual) < 15):
        print("  Movers are dominated by Pro/Gly + non-local contacts. The residual")
        print("  local-explainable set is too small to support the v2 contribution.")
        print("  -> Local-only premise is the binding constraint. (B) and per-cluster")
        print("     conditioning are unlikely to help; reconsider scope.")
    else:
        print("  A non-trivial residual of local, non-Pro/Gly movers exists -> there")
        print("  IS territory for a local higher-order model. (B) predict-delta-directly")
        print("  on this residual is worth a bounded try.")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cats = ["Pro/Gly", "non-local\n(not P/G)", "residual"]
    mv_vals = [len(expl_pg), len(expl_nl), len(residual)]
    fig, ax = plt.subplots(figsize=(6, 4.3))
    ax.bar(cats, mv_vals, color=["#54A24B", "#E45756", "#4C78A8"])
    for i, v in enumerate(mv_vals):
        ax.text(i, v + 0.3, f"{v}\n{v/len(mv):.0%}", ha="center", va="bottom")
    ax.set_ylabel("mover count"); ax.set_title(f"Composition of {len(mv)} T4L movers")
    fig.tight_layout(); fig.savefig("mover_composition.png", dpi=130)
    print("\nplot -> mover_composition.png")


if __name__ == "__main__":
    main()
