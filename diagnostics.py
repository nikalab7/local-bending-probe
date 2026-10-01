"""
proteinX -- label diagnostics: is the mover label calibrated, and how reliable is it?

  1. NULL CALIBRATION  WT-vs-WT pseudo-mutants (pairs.null_rows) carry no
     mutation, so |z| > Z_MOVER should happen ~4.6% of the time. Reported by
     number of WT crystals behind the floor, where small-sample bias shows.
  2. SPLIT-HALF RELIABILITY  each protein's crystals (WT and mutant alike)
     are split at random into halves; labels are rebuilt on each half
     independently. Mutations labelled in both halves give two independent
     measurements of one effect: correlation of z, mover agreement, kappa,
     and the AUC of |z_B| for predicting mover_A -- what a predictor that
     knew each effect only as well as half the data could reach. This bounds
     the accuracy any model can show against these labels.

Only mutations with >= 2 clean mutant crystals in forms with >= 2 * MIN_WT
WT crystals survive the split, so the reliability numbers describe the
best-measured part of the data; single-crystal labels are noisier.

Usage: python diagnostics.py [--seeds 3]  -> results/diagnostics.json
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import numpy as np
import pairs
from stats_utils import cluster_auc_ci

BINS = ((2, 4), (5, 9), (10, 19), (20, 10 ** 6))


def memoize_parsing():
    """Parse each PDB file once per process (split-half rebuilds labels 2x/seed)."""
    cache, parse = {}, pairs.parse_structure

    def cached(path):
        if path not in cache:
            cache[path] = parse(path)
        return cache[path]
    pairs.parse_structure = cached


def null_calibration(zthr):
    null = []
    data = pairs.load_proteins(null_out=null)
    rows = [x for rs, _ in data.values() for x in rs]
    zn = np.abs([x["z"] for x in null]); nn = np.array([x["n_wt"] for x in null])
    zr = np.abs([x["z"] for x in rows]); nr = np.array([x["n_wt"] for x in rows])
    out = dict(n_null=len(null), fp=float((zn > zthr).mean()), fp_3=float((zn > 3).mean()),
               n_real=len(rows), real_rate=float((zr > zthr).mean()), by_n_wt=[])
    for a, b in BINS:
        k, kr = (nn >= a) & (nn <= b), (nr >= a) & (nr <= b)
        out["by_n_wt"].append(dict(n_wt=f"{a}-{b if b < 10 ** 6 else ''}",
                                   null_n=int(k.sum()), fp=float((zn[k] > zthr).mean()) if k.any() else None,
                                   real_n=int(kr.sum()),
                                   real_rate=float((zr[kr] > zthr).mean()) if kr.any() else None))
    return out


def split_half(seeds, zthr):
    orig = pairs.build_pairs
    state = {}

    def half_build(pdb_paths, *a, **kw):
        keep = {p: f for p, f in pdb_paths.items()
                if int(hashlib.md5(f"{state['seed']}:{p}".encode()).hexdigest(), 16) % 2 == state["half"]}
        return orig(keep, *a, **kw)

    pairs.build_pairs = half_build
    res = []
    try:
        for seed in range(seeds):
            lab = {}
            for half in (0, 1):
                state.update(seed=seed, half=half)
                data = pairs.load_proteins()
                lab[half] = {(x["protein"], x["r"], x["mut"]): x for rs, _ in data.values() for x in rs}
            keys = sorted(set(lab[0]) & set(lab[1]))
            if len(keys) < 10:
                continue
            zA = np.array([lab[0][k]["z"] for k in keys]); zB = np.array([lab[1][k]["z"] for k in keys])
            yA, yB = np.abs(zA) > zthr, np.abs(zB) > zthr
            g = np.array([f"{k[0]}:{k[1]}" for k in keys])
            po = float((yA == yB).mean())
            pe = yA.mean() * yB.mean() + (1 - yA.mean()) * (1 - yB.mean())
            both = yA & yB
            aucA = cluster_auc_ci(yA, np.abs(zB), g, n_boot=500)
            aucB = cluster_auc_ci(yB, np.abs(zA), g, n_boot=500)
            res.append(dict(
                seed=seed, n=len(keys), proteins=len({k[0] for k in keys}),
                pearson_z=float(np.corrcoef(zA, zB)[0, 1]),
                mover_agreement=po, kappa=float((po - pe) / (1 - pe)),
                p_B_given_A=float(yB[yA].mean()) if yA.any() else None,
                direction_agreement=float(np.mean(np.sign(zA[both]) == np.sign(zB[both]))) if both.any() else None,
                oracle_auc=float((aucA[0] + aucB[0]) / 2)))
    finally:
        pairs.build_pairs = orig
    keys = ["n", "pearson_z", "mover_agreement", "kappa", "p_B_given_A", "direction_agreement", "oracle_auc"]
    mean = {k: float(np.mean([r[k] for r in res if r[k] is not None])) for k in keys} if res else {}
    if mean:   # Spearman-Brown: reliability of z built from all crystals
        r = mean["pearson_z"]
        mean["pearson_z_full_data_spearman_brown"] = 2 * r / (1 + r)
    return dict(per_seed=res, mean=mean)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()
    memoize_parsing()
    zthr = pairs.Z_MOVER
    cal = null_calibration(zthr)
    print(f"NULL CALIBRATION (WT-vs-WT pseudo-mutants, |z| > {zthr}; nominal 4.6%): "
          f"{cal['fp']:.1%} of {cal['n_null']}  (real mutations: {cal['real_rate']:.1%} of {cal['n_real']})")
    for b in cal["by_n_wt"]:
        fp = "-" if b["fp"] is None else f"{b['fp']:.1%}"
        rr = "-" if b["real_rate"] is None else f"{b['real_rate']:.1%}"
        print(f"  n_wt {b['n_wt']:>6s}: null {fp:>6s} (n={b['null_n']:4d})   real {rr:>6s} (n={b['real_n']:4d})")
    sh = split_half(args.seeds, zthr)
    m = sh["mean"]
    if m:
        print(f"SPLIT-HALF ({args.seeds} seeds, ~{m['n']:.0f} mutations measured in both halves): "
              f"corr(z) {m['pearson_z']:.2f} (all-crystal reliability ~{m['pearson_z_full_data_spearman_brown']:.2f}), "
              f"mover agreement {m['mover_agreement']:.0%}, kappa {m['kappa']:.2f}, "
              f"P(mover in B | mover in A) {m['p_B_given_A']:.0%}, direction agreement {m['direction_agreement']:.0%}, "
              f"oracle AUC {m['oracle_auc']:.2f}")
    os.makedirs("results", exist_ok=True)
    settings = {k: getattr(pairs, k) for k in ("SIGMA_EST", "SHRINK_K", "WT_LIG_FILTER",
                                                "IGNORE_ADDITIVES", "PRIOR", "Z_MOVER")}
    with open("results/diagnostics.json", "w") as fh:
        json.dump(dict(settings=settings, null_calibration=cal, split_half=sh), fh, indent=1)
    print("wrote results/diagnostics.json")


if __name__ == "__main__":
    main()
