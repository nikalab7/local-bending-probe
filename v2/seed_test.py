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

    pdb_id,chain,variant,resnum,wt_aa,mut_aa,crystal_form,protein_id
    2LZM,A,WT,,,,P3221,T4L
    1L63,A,WT,,,,P3221,T4L
    1L90,A,L99A,99,L,A,P3221,T4L

  * `variant` is the biological variant label; all rows sharing it are redundant
    crystals of the same thing and get aggregated to ONE observation.
  * `WT` is the reference variant.
  * `resnum/wt_aa/mut_aa` describe the single substitution (blank for WT).
  * `crystal_form` is any grouping you trust to carry systematic effects
    (space group, crystal form, deposition batch). Replicates within at least
    one form and at least two forms at a site are needed to estimate both
    components. Otherwise the verdict is inconclusive.
  * The harness is intentionally single-protein. A manifest with multiple
    `protein_id` values is rejected; a full multi-protein benchmark needs a
    protein/site hierarchy and a separate aggregation layer.
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


def metric_value(metric, P):
    """Return value plus QC branch; abstentions stay explicit."""
    value = metric(P)
    branch = metric.branch(P) if hasattr(metric, "branch") else "legacy"
    return (float(value) if np.isfinite(value) else float("nan"), branch)


def bootstrap_tau_ci(deltas, sigma_deltas, slopes, groups, n_boot=2000,
                     rng=None, level=0.90):
    """Cluster-bootstrap tau over protein/site groups, respecting shared WT refs."""
    rng = rng or np.random.default_rng(11)
    keys = sorted(set(groups))
    if len(keys) < 3:
        return float("nan"), float("nan")
    grouped = {key: np.flatnonzero(np.asarray(groups) == key) for key in keys}
    estimates = []
    for _ in range(n_boot):
        chosen = rng.choice(keys, size=len(keys), replace=True)
        idx = np.concatenate([grouped[key] for key in chosen])
        estimates.append(deconvolve_tau(np.asarray(deltas)[idx],
                                        np.asarray(sigma_deltas)[idx],
                                        np.asarray(slopes)[idx]))
    q = (1 - level) / 2
    return float(np.nanquantile(estimates, q)), float(np.nanquantile(estimates, 1 - q))


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
    sigma_arr = np.asarray(sigma_delta, float)
    if (len(d) < 3 or not np.isfinite(sigma_arr).all()
            or not np.isfinite(np.asarray(slope, float)).all()
            or np.any(np.asarray(slope, float) <= 1e-9)):
        return float("nan")
    s = np.asarray(slope, float)
    if not np.isfinite(s).all() or np.any(s <= 1e-9):
        return float("nan")
    noise_var = sigma_arr ** 2
    signal_gain = float(np.mean(s ** 2))
    excess = d.var(ddof=1) - float(np.mean(noise_var))
    return float(np.sqrt(max(excess, 0.0) / signal_gain))


def _null_deltas(sigma_iid, n_wt, n_mut, T, rng, sigma_syst=0.0,
                 wt_forms=None, mut_forms=None):
    """Null deltas with shared offsets for matching crystal forms.

    `sigma_syst` is a form-level random effect, not automatically an independent
    WT and mutant offset. When both sides contain the same form label, its offset
    is shared and cancels in the paired difference. This is the conservative model
    for matched forms and avoids inventing systematic noise twice.
    """
    wt_forms = list(wt_forms or ["__wt_unknown__"] * max(n_wt, 1))
    mut_forms = list(mut_forms or ["__mut_unknown__"] * max(n_mut, 1))
    if len(wt_forms) != max(n_wt, 1) or len(mut_forms) != max(n_mut, 1):
        raise ValueError("form labels must match WT/mutant crystal counts")
    labels = sorted(set(wt_forms + mut_forms))
    offsets = {label: rng.normal(0, sigma_syst, T) for label in labels}
    wt = np.column_stack([rng.normal(0, sigma_iid, T) + offsets[label]
                          for label in wt_forms])
    mut = np.column_stack([rng.normal(0, sigma_iid, T) + offsets[label]
                           for label in mut_forms])
    return np.median(mut, axis=1) - np.median(wt, axis=1)


def calibrated_threshold(sigma_m, n_wt, n_mut, fpr=TARGET_FPR, T=200_000,
                         rng=None, sigma_syst=0.0, wt_forms=None,
                         mut_forms=None):
    """|delta| cut giving `fpr`; systematic offsets survive aggregation."""
    rng = rng or np.random.default_rng(0)
    if not np.isfinite(sigma_m) or not np.isfinite(sigma_syst):
        return float("nan"), float("nan")
    d = _null_deltas(sigma_m, n_wt, n_mut, T, rng, sigma_syst,
                     wt_forms=wt_forms, mut_forms=mut_forms)
    return float(np.quantile(np.abs(d), 1 - fpr)), float(np.sqrt(np.mean(d ** 2)))


def oracle_ceiling(tau, slope, sigma_m, n_wt, n_mut, fpr=TARGET_FPR,
                   T=80_000, rng=None, sigma_syst=0.0, wt_forms=None,
                   mut_forms=None, pair_specs=None):
    """Oracle AUC: labels from noisy aggregated measurement, score = |true|."""
    from sklearn.metrics import roc_auc_score
    rng = rng or np.random.default_rng(7)
    slope_arr = np.asarray(slope, float)
    if (not np.isfinite(tau) or not np.isfinite(sigma_m)
            or not np.isfinite(sigma_syst) or not np.isfinite(slope_arr).all()
            or np.any(slope_arr <= 1e-9) or tau <= 0):
        return float("nan")
    t = rng.normal(0, tau, T)
    if pair_specs:
        pair_index = rng.integers(0, len(pair_specs), T)
        agg = np.empty(T, float)
        thresholds = np.empty(T, float)
        for i, spec in enumerate(pair_specs):
            mask = pair_index == i
            if not mask.any():
                continue
            nw, nm, wf, mf = spec
            agg[mask] = _null_deltas(sigma_m, nw, nm, int(mask.sum()), rng,
                                     sigma_syst, wf, mf)
            null = _null_deltas(sigma_m, nw, nm, 50_000, rng, sigma_syst, wf, mf)
            thresholds[mask] = np.quantile(np.abs(null), 1 - fpr)
    else:
        agg = _null_deltas(sigma_m, n_wt, n_mut, T, rng, sigma_syst,
                           wt_forms=wt_forms, mut_forms=mut_forms)
        thr, _ = calibrated_threshold(sigma_m, n_wt, n_mut, fpr, T=T, rng=rng,
                                      sigma_syst=sigma_syst,
                                      wt_forms=wt_forms, mut_forms=mut_forms)
        thresholds = np.full(T, thr)
    slope_draw = slope_arr if slope_arr.ndim == 0 else slope_arr[
        rng.integers(0, len(slope_arr), T)]
    mov = np.abs(slope_draw * t + agg) > thresholds
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

    protein_ids = {r.get("protein_id", "").strip() for r in rows}
    nonempty_proteins = {p for p in protein_ids if p}
    if len(nonempty_proteins) > 1:
        sys.exit("seed_test.py is single-protein only; split the manifest by protein_id")
    if len(nonempty_proteins) == 1 and "" in protein_ids:
        sys.exit("protein_id must be present on every row when supplied")

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
    print(f"  protein_id: {next(iter(nonempty_proteins), '<single-protein seed>')}")
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

    row_by_id = {id(r): r for r in rows}
    results = {}
    for mname, m in metrics:
        span = m.span
        print("\n" + "-" * 78)
        print(f"  METRIC: {mname}  (span {span})")
        print("-" * 78)

        metric_attempted = metric_finite = 0
        branch_counts = Counter()
        for centre, row_id in valid_at:
            r = row_by_id[row_id]
            res = loaded.get((r["pdb_id"], r["chain"]))
            P = window(res, centre, span) if res else None
            if P is None:
                continue
            metric_attempted += 1
            value, branch = metric_value(m, P)
            branch_counts[branch] += 1
            metric_finite += int(np.isfinite(value))
        coverage = metric_finite / metric_attempted if metric_attempted else 0.0
        print(f"  metric coverage: {metric_finite}/{metric_attempted} ({coverage:.1%})")
        print(f"  branches: {dict(branch_counts)}")

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

        # 2. WT -> mutant delta, one biological variant = one observation.
        # Keep the actual form labels and branch labels so the null model and
        # abstention coverage are pair-specific.
        wt_ref = {}
        for centre in site_nums:
            records = []
            for r in by_variant["WT"]:
                if (centre, id(r)) not in valid_at:
                    continue
                res = loaded.get((r["pdb_id"], r["chain"]))
                if not res:
                    continue
                P = window(res, centre, span)
                if P is not None:
                    v, branch = metric_value(m, P)
                    if np.isfinite(v):
                        records.append((v, r.get("crystal_form") or "unknown", branch))
            if records:
                wt_ref[centre] = records
        deltas, sigma_deltas, thresholds = [], [], []
        n_wt_per, n_mut_per, delta_groups, pair_specs = [], [], [], []
        slope_by_site, cross_branch_rejects = {}, 0

        # Measure the real-coordinate Jacobian before deconvolving tau. The
        # ideal slope is only a fallback for the explicit ideal gate, never for
        # a real seed verdict.
        jacobian_by_site = defaultdict(list)
        for centre in site_nums:
            for r in by_variant["WT"]:
                if (centre, id(r)) not in valid_at:
                    continue
                res = loaded.get((r["pdb_id"], r["chain"]))
                P = window(res, centre, span) if res else None
                if P is None:
                    continue
                j = jacobian(m, P)["slope"]
                if np.isfinite(j):
                    jacobian_by_site[centre].append(float(j))
        slope_pool = np.concatenate([np.asarray(v, float)
                                     for v in jacobian_by_site.values()]) \
            if jacobian_by_site else np.array([], float)
        slope_by_site = {k: float(np.median(v))
                         for k, v in jacobian_by_site.items() if v}

        for var, vrows in by_variant.items():
            if var == "WT" or not vrows[0].get("resnum"):
                continue
            centre = int(vrows[0]["resnum"])
            if centre not in wt_ref or centre not in slope_by_site:
                continue
            mut_records = []
            for r in vrows:
                if (centre, id(r)) not in valid_at:
                    continue
                res = loaded.get((r["pdb_id"], r["chain"]))
                if not res:
                    continue
                P = window(res, centre, span)
                if P is not None:
                    v, branch = metric_value(m, P)
                    if np.isfinite(v):
                        mut_records.append((v, r.get("crystal_form") or "unknown", branch))
            if not mut_records:
                continue
            wt_records = wt_ref[centre]
            branches = {x[2] for x in wt_records + mut_records}
            if mname == "v2_hybrid" and len(branches) != 1:
                cross_branch_rejects += 1
                continue
            wt_values = np.asarray([x[0] for x in wt_records])
            mut_values = np.asarray([x[0] for x in mut_records])
            wt_forms = [x[1] for x in wt_records]
            mut_forms = [x[1] for x in mut_records]
            threshold, sd_delta = calibrated_threshold(
                sigma_iid, len(wt_values), len(mut_values), fpr,
                sigma_syst=sigma_syst, wt_forms=wt_forms, mut_forms=mut_forms)
            delta = float(np.median(mut_values) - np.median(wt_values))
            deltas.append(delta)
            sigma_deltas.append(sd_delta)
            thresholds.append(threshold)
            n_wt_per.append(len(wt_values))
            n_mut_per.append(len(mut_values))
            pair_specs.append((len(wt_values), len(mut_values), wt_forms, mut_forms))
            delta_groups.append(f"{centre}")
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

        # 3-5. Pair-specific noise and real-coordinate slope distributions.
        sigma_deltas = np.asarray(sigma_deltas, float)
        thresholds = np.asarray(thresholds, float)
        slopes_for_delta = np.asarray([slope_by_site[int(g)] for g in delta_groups], float)
        tau = deconvolve_tau(deltas, sigma_deltas, slopes_for_delta)
        tau_ci = bootstrap_tau_ci(deltas, sigma_deltas, slopes_for_delta,
                                  delta_groups)
        n_wt = int(np.median(n_wt_per)) if n_wt_per else 1
        n_mut = int(np.median(n_mut_per)) if n_mut_per else 1
        print(f"\n  [3] tau (true-effect SD, deconvolved)  = {tau:.3f} deg"
              f"   [real Jacobian slope median {np.median(slopes_for_delta):.3f}]")
        print(f"      cluster-bootstrap {90:.0f}% CI = [{tau_ci[0]:.3f}, {tau_ci[1]:.3f}]")
        print(f"  [4] pair-specific mover threshold ({fpr:.0%} FPR, median counts"
              f" {n_wt} WT / {n_mut} mut) = {np.median(thresholds):.3f} deg")
        movers = int((np.abs(deltas) > thresholds).sum())
        print(f"      movers by that cut : {movers}/{len(deltas)} "
              f"({movers/len(deltas):.0%})")
        ceil = oracle_ceiling(tau, slope_pool, sigma_iid, n_wt, n_mut, fpr,
                              sigma_syst=sigma_syst, pair_specs=pair_specs)
        print(f"  [5] ORACLE AUC CEILING = {ceil:.3f}")
        if mname == "v2_hybrid":
            print(f"      cross-branch pairs rejected: {cross_branch_rejects}")
        results[mname] = dict(sigma_iid=sigma_iid, sigma_syst=sigma_syst,
                              sigma=sigma_m, n=len(deltas), tau=tau,
                              tau_ci=tau_ci, ceiling=ceil, movers=movers,
                              slope_pool=slope_pool,
                              counts_uniform=(len(set(zip(n_wt_per, n_mut_per))) == 1),
                              cross_branch_rejects=cross_branch_rejects)

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
    if (v2["n"] < 20 or not v2["counts_uniform"]
            or not np.isfinite(v2["tau_ci"][0])
            or v2["tau_ci"][0] <= 0):
        print("  INCONCLUSIVE -- the decision rule needs at least 20 biological")
        print("  variants, a positive cluster-bootstrap tau lower bound, and exact")
        print("  pair-specific crystal counts. The numbers above are descriptive")
        print("  when any of those conditions is not met.")
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
