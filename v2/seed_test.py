"""
v2 EMPIRICAL SEED TEST -- the decision the ideal-geometry gate cannot make.

`measurement_gate.py` established that the replacement metric construction is
sound on ideal geometry. It did NOT establish that the label is usable on real
crystal data (see CORRECTIONS.md C3). This module is the test that does, on a
small curated seed rather than a full RCSB/SIFTS miner -- because if it fails,
the miner is wasted work.

What it measures, for v1 AND v2 metrics side by side
----------------------------------------------------
  1  within-variant repeatability on redundant crystals of the SAME variant,
     decomposed into an iid component and a systematic (between-crystal-form)
     component -- the decomposition CORRECTIONS.md C2 says is required, because
     only the iid part shrinks as 1/sqrt(n)
  2  the observed WT -> mutant delta distribution, one biological variant = one
     observation (AUDIT.md F3)
  3  tau, the true-effect SD, deconvolved from the observed spread
  4  the null mover FPR at a calibrated threshold
  5  the resulting oracle AUC ceiling -- does headroom for ML exist at all
  6  the real-coordinate Jacobian (perturbation.py), which tests whether the
     ideal-geometry slopes transfer to actual windows

Decision rule, fixed in advance
-------------------------------
  * v2 within-variant noise small vs the delta distribution, and oracle ceiling
    >= 0.75  ->  GO: build the full RCSB/SIFTS miner.
  * noise comparable to effect size, ceiling ~0.5-0.6  ->  NO-GO: the metric or
    the experimental target has to change. Do not build the miner.

Usage
-----
    python3 seed_test.py --selftest
        Validate the ESTIMATORS on synthetic variants with known tau and known
        iid/systematic split. No data needed. Run this first.

    python3 seed_test.py --pdb-dir ./seed_pdb --variants variants.csv
        Run on real data.

`variants.csv` columns (curated by hand for the first iteration -- exact
entity/chain mapping matters more than automation at n=20-50):

    pdb_id,chain,variant,resnum,wt_aa,mut_aa,crystal_form
    2LZM,A,WT,,,,P3221
    1L63,A,WT,,,,P3221
    1L90,A,L99A,99,L,A,P3221

  * `variant` is the biological variant label; all rows sharing it are redundant
    crystals of the same thing and get aggregated to ONE observation.
  * `WT` is the reference variant.
  * `resnum/wt_aa/mut_aa` describe the single substitution (blank for WT).
  * `crystal_form` is any grouping you trust to carry systematic effects
    (space group, crystal form, deposition batch). Replicates within at least
    one form and at least two forms at a site are needed to estimate both
    components. Otherwise the verdict is inconclusive.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

import numpy as np

from metrics import HybridAxisBend, LegacyPCA5
from perturbation import jacobian

TARGET_FPR = 0.05
THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'}


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #
def parse_ca(path, want_chain=None):
    """Parse fixed-column PDB/ENT into {resnum: (aa, xyz, bfactor)}.

    Differs from v1's parser in two ways that matter here: it keeps the B-factor
    (needed for a per-structure uncertainty estimate) and counts skipped
    insertion codes. Windows touching an insertion code must be checked by
    sequence registration before analysis.
    """
    if path.lower().endswith(".cif"):
        raise ValueError("mmCIF parsing is not implemented; provide .pdb or .ent files")
    best, ins = {}, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            aa = THREE2ONE.get(line[17:20].strip())
            if aa is None:
                continue
            if line[26] != " ":
                ins += 1
                continue
            ch = line[21]
            if want_chain and ch != want_chain:
                continue
            try:
                rs = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]),
                                float(line[46:54])])
                occ = float(line[54:60]) if line[54:60].strip() else 1.0
                bf = float(line[60:66]) if line[60:66].strip() else np.nan
            except ValueError:
                continue
            if rs not in best or occ > best[rs][0]:
                best[rs] = (occ, aa, xyz, bf)
    return {rs: (aa, xyz, bf) for rs, (_, aa, xyz, bf) in best.items()}, ins


def window(res, centre, span):
    """(span,3) Ca centred on `centre`, or None if incomplete/discontinuous."""
    half = span // 2
    idx = [centre + k for k in range(-half, span - half)]
    if any(i not in res for i in idx):
        return None
    P = np.array([res[i][1] for i in idx])
    d = np.linalg.norm(np.diff(P, axis=0), axis=1)
    if d.max() > 4.5 or d.min() < 2.5:
        return None
    return P


def window_bfactors(res, centre, span):
    half = span // 2
    idx = [centre + k for k in range(-half, span - half)]
    bs = [res[i][2] for i in idx if i in res]
    bs = [b for b in bs if np.isfinite(b)]
    return float(np.mean(bs)) if bs else float("nan")


def window_sequence(res, centre, span):
    half = span // 2
    idx = range(centre - half, centre + span - half)
    return tuple(res[i][0] for i in idx) if all(i in res for i in idx) else None


# --------------------------------------------------------------------------- #
# estimators  (validated by --selftest)
# --------------------------------------------------------------------------- #
def variance_components(groups):
    """Split within-variant spread into iid and systematic parts.

    `groups` maps crystal_form -> list of metric values for ONE variant.
    One-way random-effects decomposition:
        sigma^2_iid   = pooled within-form variance
        sigma^2_syst  = between-form variance of the form means, corrected for
                        the sampling variance of those means
    Returns (sigma_iid, sigma_syst, n_forms, n_total). sigma_syst is nan when
    there is only one form -- unestimated, not zero.
    """
    forms = {k: np.asarray(v, float) for k, v in groups.items() if len(v) >= 1}
    n_total = sum(len(v) for v in forms.values())
    multi = {k: v for k, v in forms.items() if len(v) >= 2}
    if multi:
        ss = sum(((v - v.mean()) ** 2).sum() for v in multi.values())
        df = sum(len(v) - 1 for v in multi.values())
        s2_iid = ss / df if df > 0 else float("nan")
    else:
        s2_iid = float("nan")
    if len(forms) >= 2:
        means = np.array([v.mean() for v in forms.values()])
        ns = np.array([len(v) for v in forms.values()], float)
        s2_between = means.var(ddof=1)
        corr = np.nanmean(s2_iid / ns) if np.isfinite(s2_iid) else 0.0
        s2_syst = max(s2_between - corr, 0.0)
    else:
        s2_syst = float("nan")
    return (float(np.sqrt(s2_iid)) if np.isfinite(s2_iid) else float("nan"),
            float(np.sqrt(s2_syst)) if np.isfinite(s2_syst) else float("nan"),
            len(forms), n_total)


def deconvolve_tau(deltas, sigma_delta, slope=1.0):
    """True-effect SD from the observed delta spread.

        var(delta_obs) = slope^2 * tau^2 + sigma_delta^2
    Returns 0.0 when the observed spread does not exceed the noise -- which is
    itself the NO-GO answer, not a failure to compute.
    """
    d = np.asarray([x for x in deltas if np.isfinite(x)], float)
    if len(d) < 3 or not np.isfinite(sigma_delta) or slope <= 1e-9:
        return float("nan")
    excess = d.var(ddof=1) - sigma_delta ** 2
    return float(np.sqrt(max(excess, 0.0)) / slope)


def _null_deltas(sigma_iid, n_wt, n_mut, T, rng, sigma_syst=0.0):
    """Independent crystal noise plus one persistent offset per variant."""
    wt = np.median(rng.normal(0, sigma_iid, (T, max(n_wt, 1))), axis=1)
    mut = np.median(rng.normal(0, sigma_iid, (T, max(n_mut, 1))), axis=1)
    if sigma_syst:
        wt += rng.normal(0, sigma_syst, T)
        mut += rng.normal(0, sigma_syst, T)
    return mut - wt


def calibrated_threshold(sigma_m, n_wt, n_mut, fpr=TARGET_FPR, T=200_000,
                         rng=None, sigma_syst=0.0):
    """|delta| cut giving `fpr`; systematic offsets survive aggregation."""
    rng = rng or np.random.default_rng(0)
    if not np.isfinite(sigma_m) or not np.isfinite(sigma_syst):
        return float("nan"), float("nan")
    d = _null_deltas(sigma_m, n_wt, n_mut, T, rng, sigma_syst)
    return float(np.quantile(np.abs(d), 1 - fpr)), float(np.sqrt(np.mean(d ** 2)))


def oracle_ceiling(tau, slope, sigma_m, n_wt, n_mut, fpr=TARGET_FPR,
                   T=80_000, rng=None, sigma_syst=0.0):
    """Oracle AUC: labels from noisy aggregated measurement, score = |true|."""
    from sklearn.metrics import roc_auc_score
    rng = rng or np.random.default_rng(7)
    if not all(np.isfinite([tau, slope, sigma_m, sigma_syst])) or slope <= 1e-9 or tau <= 0:
        return float("nan")
    t = rng.normal(0, tau, T)
    agg = _null_deltas(sigma_m, n_wt, n_mut, T, rng, sigma_syst)
    thr, _ = calibrated_threshold(sigma_m, n_wt, n_mut, fpr, T=T, rng=rng,
                                  sigma_syst=sigma_syst)
    mov = np.abs(slope * t + agg) > thr
    if not (0 < mov.sum() < T):
        return float("nan")
    return float(roc_auc_score(mov, np.abs(t)))


# --------------------------------------------------------------------------- #
# the real-data run
# --------------------------------------------------------------------------- #
def run(pdb_dir, variants_csv, fpr=TARGET_FPR):
    with open(variants_csv, newline="") as source:
        rows = list(csv.DictReader(source))
    if not rows:
        sys.exit("variants.csv is empty")
    has_form = any(r.get("crystal_form") for r in rows)

    metrics = [("v1_pca5", LegacyPCA5()), ("v2_hybrid", HybridAxisBend(9))]

    # variant -> list of (pdb, chain, form, resnum, wt_aa, mut_aa)
    by_variant = defaultdict(list)
    for r in rows:
        by_variant[r["variant"].strip()].append(r)
    if "WT" not in by_variant:
        sys.exit("variants.csv must contain rows with variant=WT")

    print("=" * 78)
    print("  v2 EMPIRICAL SEED TEST")
    print("=" * 78)
    print(f"  variants: {len(by_variant)}  (WT crystals: {len(by_variant['WT'])})")
    print(f"  crystal_form column present: {has_form}"
          + ("" if has_form else "  -> systematic component NOT estimable"))

    # load every structure once
    loaded, ins_total = {}, 0
    for r in rows:
        key = (r["pdb_id"], r["chain"])
        if key in loaded:
            continue
        path = None
        for ext in (".pdb", ".ent"):
            cand = os.path.join(pdb_dir, r["pdb_id"] + ext)
            if os.path.exists(cand):
                path = cand
                break
        if path is None:
            if os.path.exists(os.path.join(pdb_dir, r["pdb_id"] + ".cif")):
                print(f"    UNSUPPORTED {r['pdb_id']}.cif -- convert to PDB/ENT first")
            else:
                print(f"    MISSING {r['pdb_id']} -- skipped")
            continue
        res, ins = parse_ca(path, r["chain"])
        ins_total += ins
        loaded[key] = res
    print(f"  structures loaded: {len(loaded)}   insertion-coded Ca skipped: "
          f"{ins_total}")

    # Screen local registration at the widest metric span. A PDB author's
    # residue number alone is not proof that two windows refer to the same site;
    # full entity/SIFTS mapping is still required before a benchmark is built.
    site_meta = {}
    for var, vrows in by_variant.items():
        if var == "WT":
            continue
        meta = {(r.get("resnum", "").strip(), r.get("wt_aa", "").strip().upper(),
                 r.get("mut_aa", "").strip().upper()) for r in vrows}
        if len(meta) != 1 or not all(next(iter(meta))):
            sys.exit(f"incomplete or inconsistent mutation metadata for {var}")
        num, wt_aa, mut_aa = next(iter(meta))
        if wt_aa == mut_aa:
            sys.exit(f"variant {var} has identical WT and mutant amino acids")
        centre = int(num)
        if centre in site_meta and site_meta[centre] != wt_aa:
            sys.exit(f"conflicting WT amino acid at residue {centre}")
        site_meta[centre] = wt_aa

    valid_at = set()
    rejected = 0
    for centre, wt_aa in site_meta.items():
        wt_sequences = []
        for r in by_variant["WT"]:
            res = loaded.get((r["pdb_id"], r["chain"]))
            seq = window_sequence(res, centre, 9) if res else None
            if seq is not None and seq[4] == wt_aa and window(res, centre, 9) is not None:
                wt_sequences.append(seq)
        counts = Counter(wt_sequences)
        if not counts:
            rejected += len(by_variant["WT"])
            continue
        highest = max(counts.values())
        modes = [seq for seq, count in counts.items() if count == highest]
        if len(modes) != 1:
            print(f"    AMBIGUOUS WT sequence at {centre} -- site skipped")
            rejected += len(by_variant["WT"])
            continue
        reference = modes[0]
        for r in by_variant["WT"]:
            res = loaded.get((r["pdb_id"], r["chain"]))
            seq = window_sequence(res, centre, 9) if res else None
            if seq == reference and window(res, centre, 9) is not None:
                valid_at.add((centre, id(r)))
            else:
                rejected += 1
        for var, vrows in by_variant.items():
            if var == "WT":
                continue
            for r in vrows:
                if int(r["resnum"]) != centre:
                    continue
                res = loaded.get((r["pdb_id"], r["chain"]))
                seq = window_sequence(res, centre, 9) if res else None
                expected = reference[:4] + (r["mut_aa"].strip().upper(),) + reference[5:]
                if seq == expected and window(res, centre, 9) is not None:
                    valid_at.add((centre, id(r)))
                else:
                    rejected += 1
    print(f"  local sequence QC: {len(valid_at)} site/structure matches; "
          f"{rejected} rejected (9-Ca window, WT/mutant identity and flanks)")

    results = {}
    for mname, m in metrics:
        span = m.span
        print("\n" + "-" * 78)
        print(f"  METRIC: {mname}  (span {span})")
        print("-" * 78)

        # ---- per-variant values at each mutated site ----
        sites = sorted({(int(r["resnum"]), r["variant"])
                        for r in rows if r.get("resnum")})
        site_nums = sorted({s[0] for s in sites})

        # 1. within-variant repeatability, per site, decomposed
        iid_all, syst_all = [], []
        for var, vrows in by_variant.items():
            centres = ([int(vrows[0]["resnum"])] if vrows[0].get("resnum")
                       else site_nums)
            for centre in centres:
                groups = defaultdict(list)
                for r in vrows:
                    if (centre, id(r)) not in valid_at:
                        continue
                    res = loaded.get((r["pdb_id"], r["chain"]))
                    if not res:
                        continue
                    P = window(res, centre, span)
                    if P is None:
                        continue
                    v = m(P)
                    if np.isfinite(v):
                        groups[r.get("crystal_form") or "unknown"].append(v)
                if sum(len(x) for x in groups.values()) < 2:
                    continue
                s_iid, s_syst, nf, nt = variance_components(groups)
                if np.isfinite(s_iid):
                    iid_all.append(s_iid)
                if np.isfinite(s_syst):
                    syst_all.append(s_syst)
        sigma_iid = float(np.median(iid_all)) if iid_all else float("nan")
        sigma_syst = float(np.median(syst_all)) if syst_all else float("nan")
        sigma_m = float(np.hypot(sigma_iid, sigma_syst)) \
            if np.isfinite(sigma_iid) and np.isfinite(sigma_syst) else float("nan")
        print(f"  [1] within-variant repeatability  (n groups: iid {len(iid_all)},"
              f" syst {len(syst_all)})")
        print(f"      sigma_iid        = {sigma_iid:.3f} deg   (shrinks as 1/sqrt(n))")
        print(f"      sigma_systematic = {sigma_syst:.3f} deg   (does NOT shrink)")
        print(f"      sigma_total      = {sigma_m:.3f} deg")

        # 2. WT -> mutant delta, one variant = one observation
        wt_ref, wt_counts = {}, {}
        for centre in site_nums:
            vals = []
            for r in by_variant["WT"]:
                if (centre, id(r)) not in valid_at:
                    continue
                res = loaded.get((r["pdb_id"], r["chain"]))
                if not res:
                    continue
                P = window(res, centre, span)
                if P is not None:
                    v = m(P)
                    if np.isfinite(v):
                        vals.append(v)
            if vals:
                wt_ref[centre] = float(np.median(vals))
                wt_counts[centre] = len(vals)
        deltas, n_wt_per, n_mut_per = [], [], []
        for var, vrows in by_variant.items():
            if var == "WT" or not vrows[0].get("resnum"):
                continue
            centre = int(vrows[0]["resnum"])
            if centre not in wt_ref:
                continue
            vals = []
            for r in vrows:
                if (centre, id(r)) not in valid_at:
                    continue
                res = loaded.get((r["pdb_id"], r["chain"]))
                if not res:
                    continue
                P = window(res, centre, span)
                if P is not None:
                    v = m(P)
                    if np.isfinite(v):
                        vals.append(v)
            if vals:
                deltas.append(float(np.median(vals)) - wt_ref[centre])
                n_wt_per.append(wt_counts[centre])
                n_mut_per.append(len(vals))
        deltas = np.array(deltas, float)
        print(f"\n  [2] WT -> mutant delta   (one biological variant = one obs)")
        if len(deltas) == 0:
            print("      no scorable variants -- cannot continue for this metric")
            continue
        print(f"      variants scored : {len(deltas)}")
        print(f"      median |delta|  : {np.median(np.abs(deltas)):.3f} deg")
        delta_sd = float(deltas.std(ddof=1)) if len(deltas) >= 2 else float("nan")
        print(f"      SD of delta     : {delta_sd:.3f} deg")
        print(f"      90th pct |delta|: {np.percentile(np.abs(deltas), 90):.3f} deg")

        # 3-5
        n_wt = int(np.median(n_wt_per)) if n_wt_per else 1
        n_mut = int(np.median(n_mut_per)) if n_mut_per else 1
        thr, sd_delta = calibrated_threshold(sigma_iid, n_wt, n_mut, fpr,
                                             sigma_syst=sigma_syst)
        slope = 1.0 if mname == "v2_hybrid" else 0.316   # ideal-geometry value
        tau = deconvolve_tau(deltas, sd_delta, slope)
        print(f"\n  [3] tau (true-effect SD, deconvolved)  = {tau:.3f} deg"
              f"   [slope assumed {slope}]")
        print(f"  [4] approximate mover threshold ({fpr:.0%} FPR, median counts"
              f" {n_wt} WT / {n_mut} mut) = {thr:.3f} deg")
        movers = int((np.abs(deltas) > thr).sum())
        print(f"      movers by that cut : {movers}/{len(deltas)} "
              f"({movers/len(deltas):.0%})")
        ceil = oracle_ceiling(tau, slope, sigma_iid, n_wt, n_mut, fpr,
                              sigma_syst=sigma_syst)
        print(f"  [5] ORACLE AUC CEILING = {ceil:.3f}")
        counts_uniform = (len(set(zip(n_wt_per, n_mut_per))) == 1)
        results[mname] = dict(sigma_iid=sigma_iid, sigma_syst=sigma_syst,
                              sigma=sigma_m, n=len(deltas), tau=tau,
                              ceiling=ceil, movers=movers,
                              counts_uniform=counts_uniform)

        # 6. real-coordinate Jacobian
        js = []
        for centre in site_nums[:40]:
            for r in by_variant["WT"]:
                if (centre, id(r)) not in valid_at:
                    continue
                res = loaded.get((r["pdb_id"], r["chain"]))
                if not res:
                    continue
                P = window(res, centre, span)
                if P is None:
                    continue
                j = jacobian(m, P)["slope"]
                if np.isfinite(j):
                    js.append(j)
        if js:
            js = np.array(js)
            print(f"\n  [6] real-coordinate Jacobian (n={len(js)} windows)")
            print(f"      median {np.median(js):.3f}   IQR "
                  f"[{np.percentile(js,25):.3f}, {np.percentile(js,75):.3f}]")
            print("      compare against the ideal-geometry value from")
            print("      perturbation.py --selftest; a large gap means the")
            print("      ideal-geometry slopes do NOT transfer.")

    # ---- verdict ----
    print("\n" + "=" * 78)
    print("  SEED VERDICT")
    print("=" * 78)
    v2 = results.get("v2_hybrid")
    if not v2 or not np.isfinite(v2["ceiling"]):
        print("  INCONCLUSIVE -- too few scorable variants or the iid/systematic")
        print("  noise split is unestimated. Add matched crystals across forms")
        print("  before using a GO/NO-GO verdict.")
        return results
    if v2["n"] < 20 or not v2["counts_uniform"]:
        print("  INCONCLUSIVE -- the decision rule needs at least 20 biological")
        print("  variants and exact pair-specific crystal counts. The numbers")
        print("  above use median counts when replication varies by pair.")
        return results
    print(f"  v2 sigma_total {v2['sigma']:.2f} deg   tau {v2['tau']:.2f} deg   "
          f"ceiling {v2['ceiling']:.3f}   n={v2['n']} variants")
    if np.isfinite(v2["sigma_syst"]) and v2["sigma_syst"] > v2["sigma_iid"]:
        print("  NOTE systematic > iid: more crystals per variant will NOT help.")
    if v2["ceiling"] >= 0.75 and v2["tau"] > v2["sigma"]:
        print("\n  PROVISIONAL GO -- simulated headroom exists. Review the")
        print("  real-coordinate Jacobian and crystal-form matching before")
        print("  building the full RCSB/SIFTS miner.")
    elif v2["ceiling"] >= 0.65:
        print("\n  PROVISIONAL MARGINAL -- some headroom but thin. Widen the seed before")
        print("  committing to the miner; consider a higher-resolution subset.")
    else:
        print("\n  PROVISIONAL NO-GO -- under the assumed response slope and")
        print("  independent form offsets, the measurement does not resolve")
        print("  the effect. Review the real-coordinate Jacobian and form")
        print("  matching before treating this as a final decision.")
    return results


# --------------------------------------------------------------------------- #
# estimator self-test -- runs with no data
# --------------------------------------------------------------------------- #
def selftest():
    rng = np.random.default_rng(0)
    print("=" * 78)
    print("  seed_test.py estimator self-test (no data required)")
    print("=" * 78)

    print("\n  1. variance_components: recover a known iid / systematic split")
    print(f"      {'true iid':>9} {'true syst':>10} {'est iid':>9} {'est syst':>9}")
    for s_iid, s_syst in ((1.0, 0.0), (1.0, 0.5), (1.0, 2.0), (0.5, 1.5)):
        ei, es = [], []
        for _ in range(300):
            groups = {}
            for f in range(6):                       # 6 crystal forms
                off = rng.normal(0, s_syst)
                groups[f] = list(rng.normal(off, s_iid, 4))
            a, b, _, _ = variance_components(groups)
            ei.append(a)
            es.append(b)
        print(f"      {s_iid:>9.2f} {s_syst:>10.2f} {np.nanmean(ei):>9.2f} "
              f"{np.nanmean(es):>9.2f}")
    print("      -> iid recovered tightly; systematic recovered with more spread")
    print("         (6 forms is a small sample for a between-group variance).")

    print("\n  2. deconvolve_tau: recover a known true-effect SD")
    print(f"      {'true tau':>9} {'sigma_d':>8} {'est tau':>9} {'bias':>8}")
    for tau in (0.0, 1.0, 3.0, 8.0):
        for sd in (1.0, 3.0):
            est = []
            for _ in range(400):
                d = rng.normal(0, tau, 40) + rng.normal(0, sd, 40)
                est.append(deconvolve_tau(d, sd))
            e = np.nanmean(est)
            print(f"      {tau:>9.2f} {sd:>8.2f} {e:>9.2f} {e-tau:>+8.2f}")
    print("      -> unbiased when tau > sigma_d; floors at 0 and reads slightly")
    print("         HIGH when tau << sigma_d (a max(.,0) on a noisy difference).")
    print("         So a small non-zero tau estimate is NOT evidence of signal.")

    print("\n  3. calibrated_threshold: does it deliver the target FPR?")
    print(f"      {'n_wt':>5} {'n_mut':>6} {'k':>7} {'realised FPR':>13}")
    for n_wt, n_mut in ((1, 1), (5, 2), (20, 5)):
        thr, _ = calibrated_threshold(1.0, n_wt, n_mut)
        wt = np.median(rng.normal(0, 1, (100000, n_wt)), axis=1)
        mut = np.median(rng.normal(0, 1, (100000, n_mut)), axis=1)
        print(f"      {n_wt:>5} {n_mut:>6} {thr:>7.3f} "
              f"{np.mean(np.abs(mut-wt) > thr):>13.1%}")

    print("\n  4. oracle_ceiling: monotone in tau/sigma (tau=0 is undefined -> nan,")
    print("     since a perfect predictor of a zero effect has nothing to rank)")
    print(f"      {'tau':>6} {'sigma':>7} {'ceiling':>9}")
    for tau in (0.0, 0.5, 1.0, 3.0, 10.0):
        print(f"      {tau:>6.1f} {1.0:>7.1f} "
              f"{oracle_ceiling(tau, 1.0, 1.0, 5, 2):>9.3f}")
    print("\n  The ESTIMATORS are validated. What is NOT validated is that they")
    print("  will be applied to real coordinates that behave like the assumptions")
    print("  -- that is what running this on a real seed set is for.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--pdb-dir")
    ap.add_argument("--variants")
    ap.add_argument("--fpr", type=float, default=TARGET_FPR)
    a = ap.parse_args()
    if a.selftest:
        selftest()
    elif a.pdb_dir and a.variants:
        run(a.pdb_dir, a.variants, a.fpr)
    else:
        ap.print_help()
