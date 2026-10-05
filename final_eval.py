"""
proteinX -- the frozen final evaluation (see PROTOCOL.md; do not edit after
PROTOCOL.md is committed).

  python final_eval.py dev
      Dev data only (pairs.load_proteins). Leave-family-out OOF scores
      (10 folds of >= 30%-identity families), every analysis below.
      -> results/final_dev.json

  python final_eval.py lockbox
      Train the frozen model on ALL dev rows, score the lockbox
      (pairs.load_lockbox) once. Refuses to run if PROTOCOL.md or any frozen
      source file has uncommitted changes, or if results/final_lockbox.json
      already exists (one run only).
      -> results/final_lockbox.json

Frozen pipeline
  label     pairs.METRIC = "multi" (max |z| over bend, phi, psi, CA torsion b),
            NCS averaging, threshold MULTI_T = 2.61 (dev-null calibrated),
            mover = |z| > Z_MOVER after rescaling
  features  FEATURES (WT structure + substitution only; no mutant-side input)
  model     logistic regression C = 0.3, balanced classes, median impute,
            standardize, confidence weights (delta_model.label_weights)

Analyses (both modes)
  headline      AUC + family-bootstrap 95% CI; paired family-bootstrap dAUC
                vs the SS-only and burial-only baselines (same model class)
  ceiling       0.787 = split-half oracle AUC of the dev label (diagnostics)
  within        within-protein AUC (+ CI)
  T4L           with and without T4L
  controls      positive: same pipeline/splits predicting helix (non-SS
                features); negative: labels shuffled within families
  power         Kish effective number of families; minimum detectable AUC
                and dAUC (alpha 0.05 two-sided, power 0.8)
  diagnostics   artifact-suspect exclusions (temperature mismatch, resolution
                gap, mutant altlocs, lattice contact) -- never model inputs
  subsets       confident labels, SS classes (secondary)
Groups with fewer than MIN_GROUP rows or MIN_CLASS per class are not
reported as results ("n too small").
"""
from __future__ import annotations
import argparse
import json
import os
import pickle
import subprocess
import numpy as np

import pairs
import delta_model as dm
from stats_utils import cluster_auc_ci, within_auc_ci, fast_auc

# ------------------------------ frozen settings -----------------------------
FEATURES = dm.SITE + dm.SUBST
BASELINES = {"ss_only": dm.SS, "burial_only": ["hse_up", "n_ca10"]}
POSITIVE_FEATURES = [c for c in dm.SITE if c not in dm.SS] + dm.SUBST
CEILING = 0.787                     # dev split-half oracle AUC (WORKLOG v6)
LEVEL = 0.95
N_BOOT = 2000
N_SHUFFLE_DEV = 20                  # retrained negative controls
N_SHUFFLE_LOCKBOX = 1000            # frozen-model negative controls
MIN_GROUP, MIN_CLASS = 30, 10
MIN_FAMILIES = 10                   # lockbox size rule (PROTOCOL.md section 1)
Z_A, Z_B = 1.96, 0.8416             # alpha 0.05 two-sided, power 0.8
FROZEN_FILES = ("PROTOCOL.md", "final_eval.py", "pairs.py", "delta_model.py",
                "structure_features.py", "plm_features.py", "stats_utils.py",
                "mine_pairs.py", "manifests/mined_ms4.json", "manifests/lockbox.json")
CACHE = os.path.join("results", "_cache")
OUT = {"dev": os.path.join("results", "final_dev.json"),
       "lockbox": os.path.join("results", "final_lockbox.json")}


def freeze_label():
    pairs.METRIC = "multi"
    pairs.NCS_AVERAGE = True
    assert pairs.MULTI_T == 2.61 and pairs.Z_MOVER == 2.0
    assert pairs.MULTI_METRICS == ("bend", "phi", "psi", "ca_tor_b")


# ------------------------------ data ----------------------------------------
def load(which):
    """(rows, null rows, features) for 'dev' or 'lockbox', cached on disk."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f"{which}_rows.pkl")
    if os.path.exists(path):
        with open(path, "rb") as fh:
            return pickle.load(fh)
    freeze_label()
    import diagnostics
    diagnostics.memoize_parsing()
    null = []
    data = pairs.load_proteins(null_out=null) if which == "dev" else pairs.load_lockbox(null_out=null)
    rows = [x for rs, _ in data.values() for x in rs]
    F = dm.featurize(rows)
    out = (rows, null, F)
    with open(path, "wb") as fh:
        pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return out


def arrays(rows):
    z = np.array([x["z"] for x in rows], float)
    return dict(z=z, y=np.abs(z) > pairs.Z_MOVER, w=dm.label_weights(z, pairs.Z_MOVER),
                fam=np.array([str(x["family"]) for x in rows]),
                prot=np.array([str(x["protein"]) for x in rows]))


# ------------------------------ statistics ----------------------------------
def auc_ci(y, s, fam):
    a, lo, hi, nf = cluster_auc_ci(y, s, fam, n_boot=N_BOOT, level=LEVEL)
    return dict(auc=a, lo=lo, hi=hi, n=int(len(y)), movers=int(np.sum(y)), families=nf)


def _boot(fam, seed):
    """Index arrays of family-bootstrap resamples."""
    groups = [np.flatnonzero(fam == f) for f in np.unique(fam)]
    rng = np.random.default_rng(seed)
    for _ in range(N_BOOT):
        yield np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])


def _boot_delta(y, sa, sb, fam, seed):
    return np.array([fast_auc(y[i], sb[i]) - fast_auc(y[i], sa[i])
                     for i in _boot(fam, seed) if 0 < y[i].sum() < len(i)])


def paired(y, s_base, s_model, fam):
    """dAUC = AUC(model) - AUC(base), paired family bootstrap.

    p_two_sided = 2 x the smaller tail of the bootstrap dAUC distribution at 0.
    """
    d = _boot_delta(y, s_base, s_model, fam, 0)
    a = (1 - LEVEL) / 2 * 100
    lo, hi = np.percentile(d, [a, 100 - a])
    return dict(d_auc=fast_auc(y, s_model) - fast_auc(y, s_base), lo=float(lo), hi=float(hi),
                p_two_sided=float(min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))))


def boot_sd(y, s, fam, seed=1):
    """Family-bootstrap SD of the AUC (its standard error under clustering)."""
    return float(np.std([fast_auc(y[i], s[i]) for i in _boot(fam, seed) if 0 < y[i].sum() < len(i)]))


def boot_sd_delta(y, sa, sb, fam, seed=1):
    return float(np.std(_boot_delta(y, sa, sb, fam, seed)))


def power(y, s, s_base, fam):
    """Minimum detectable AUC / dAUC at this sample's family structure.

    n_eff = Kish effective number of families, (sum n_f)^2 / sum n_f^2.
    SE = family-bootstrap SD of the AUC (clustering included); the design
    effect compares it with the independent-rows Hanley-McNeil SE at 0.5.
    """
    n_f = np.array([np.sum(fam == f) for f in np.unique(fam)], float)
    n1, n0 = int(y.sum()), int((~y).sum())
    se_iid = float(np.sqrt((n1 + n0 + 1) / (12.0 * n1 * n0)))
    se = boot_sd(y, s, fam)
    se_d = boot_sd_delta(y, s_base, s, fam)
    return dict(families=int(len(n_f)), n_eff_families=float(n_f.sum() ** 2 / (n_f ** 2).sum()),
                se_auc=se, se_auc_iid=se_iid, design_effect=float((se / se_iid) ** 2),
                min_detectable_auc=float(0.5 + (Z_A + Z_B) * se),
                se_d_auc_vs_ss=se_d, min_detectable_d_auc=float((Z_A + Z_B) * se_d))


def reportable(y):
    y = np.asarray(y, bool)
    return len(y) >= MIN_GROUP and min(y.sum(), (~y).sum()) >= MIN_CLASS


def group(y, s, fam, mask, name):
    if not reportable(y[mask]):
        return dict(name=name, n=int(mask.sum()), movers=int(y[mask].sum()),
                    note=f"n too small (< {MIN_GROUP} rows or < {MIN_CLASS} per class): not reported")
    return dict(name=name, **auc_ci(y[mask], s[mask], fam[mask]))


# ------------------------------ model ---------------------------------------
def X_of(F, cols):
    return dm.matrix(F, cols)


def fit(X, y, w=None):
    return dm.fit_model("logreg", X, y, w)


def oof(X, y, folds, w=None):
    return dm.oof_predict(X, y, folds, "logreg", w)


def shuffle_within(fam, rng):
    """Permutation index that shuffles rows only within each family."""
    idx = np.arange(len(fam))
    for f in np.unique(fam):
        m = np.flatnonzero(fam == f)
        idx[m] = rng.permutation(m)
    return idx


def artifact_masks(rows, F):
    tw = np.array([x.get("temp_wt", np.nan) for x in rows], float)
    tm = np.array([x.get("temp_mut", np.nan) for x in rows], float)
    rw = np.array([x.get("res_wt", np.nan) for x in rows], float)
    rm = np.array([x.get("res_mut", np.nan) for x in rows], float)
    am = np.array([x.get("altloc_mut", 0) or 0 for x in rows], float)
    ls = X_of(F, ["lattice_site"])[:, 0]
    lw = X_of(F, ["lattice_window"])[:, 0]
    with np.errstate(invalid="ignore"):
        m = {"temperature_mismatch_gt50K": np.abs(tm - tw) > 50,
             "mutant_resolution_worse_gt0.5A": (rm - rw) > 0.5,
             "mutant_altloc_in_window": am > 0,
             "mutated_residue_in_lattice_contact": ls > 0}
    m["any"] = np.any(np.vstack(list(m.values())), axis=0)
    with np.errstate(invalid="ignore"):
        m["window_in_lattice_contact"] = lw > 0      # broad (~60% of dev rows): separate
    return m


def diagnostics_block(y, s, fam, rows, F):
    out = {}
    for k, m in artifact_masks(rows, F).items():
        out[k] = dict(n_suspect=int(m.sum()), mover_rate_suspect=float(y[m].mean()) if m.any() else None,
                      mover_rate_rest=float(y[~m].mean()) if (~m).any() else None,
                      excluded=group(y, s, fam, ~m, f"excluding {k}"))
    return out


def subsets(y, s, fam, z, F):
    conf = (np.abs(z) > dm.CONFIDENT_HI) | (np.abs(z) < dm.CONFIDENT_LO)
    ss = X_of(F, dm.SS)
    return [group(y, s, fam, conf, "confident labels (|z|>3 or <1)"),
            group(y, s, fam, ss[:, 0] == 1, "helix"),
            group(y, s, fam, ss[:, 1] == 1, "strand"),
            group(y, s, fam, ss[:, 2] == 1, "loop")]


def null_fp(null):
    zn = np.abs([x["z"] for x in null])
    return dict(n=int(len(zn)), fp=float((zn > pairs.Z_MOVER).mean()) if len(zn) else None)


def headline(y, scores, fam, prot):
    """Model and baselines: AUC + CI, paired dAUC, within-protein AUC."""
    out = {"model": auc_ci(y, scores["model"], fam), "ceiling_oracle_auc_dev": CEILING}
    for b in BASELINES:
        out[b] = auc_ci(y, scores[b], fam)
        out[f"model_vs_{b}"] = paired(y, scores[b], scores["model"], fam)
    a, lo, hi = within_auc_ci(y, scores["model"], prot, fam, n_boot=N_BOOT, level=LEVEL)
    out["within_protein"] = dict(auc=a, lo=lo, hi=hi)
    for b in BASELINES:
        a, lo, hi = within_auc_ci(y, scores[b], prot, fam, n_boot=N_BOOT, level=LEVEL)
        out[f"within_protein_{b}"] = dict(auc=a, lo=lo, hi=hi)
    return out


# ------------------------------ modes ---------------------------------------
def run_dev():
    rows, null, F = load("dev")
    A = arrays(rows); y, w, fam, prot, z = A["y"], A["w"], A["fam"], A["prot"], A["z"]
    folds = dm.site_folds(fam, dm.FAMILY_FOLDS, 0)
    X = X_of(F, FEATURES)
    scores = {"model": oof(X, y, folds, w)}
    for b, cols in BASELINES.items():
        scores[b] = oof(X_of(F, cols), y, folds, w)
    R = dict(mode="dev", n=len(y), movers=int(y.sum()), null=null_fp(null),
             headline=headline(y, scores, fam, prot))
    # T4L: (a) OOF scores restricted to non-T4L rows; (b) T4L removed entirely
    t4l = prot == "T4L"
    nt = ~t4l
    f2 = dm.site_folds(fam[nt], dm.FAMILY_FOLDS, 0)
    s2 = oof(X[nt], y[nt], f2, w[nt])
    R["t4l"] = dict(n_t4l=int(t4l.sum()),
                    t4l_rows_only=group(y, scores["model"], fam, t4l, "T4L rows (OOF)"),
                    without_t4l_rows=group(y, scores["model"], fam, nt, "non-T4L rows (OOF)"),
                    without_t4l_retrained=auc_ci(y[nt], s2, fam[nt]))
    # positive control: helix from non-SS features, same splits
    yh = X_of(F, ["ss_H"])[:, 0] == 1
    ph = oof(X_of(F, POSITIVE_FEATURES), yh, folds)
    R["positive_control_helix"] = auc_ci(yh, ph, fam)
    # negative control: labels (and their weights) shuffled within families
    rng = np.random.default_rng(0); neg = []
    for _ in range(N_SHUFFLE_DEV):
        i = shuffle_within(fam, rng)
        neg.append(fast_auc(y[i], oof(X, y[i], folds, w[i])))
    R["negative_control_shuffled_within_family"] = dict(
        reps=N_SHUFFLE_DEV, mean_auc=float(np.mean(neg)), sd=float(np.std(neg)),
        min=float(np.min(neg)), max=float(np.max(neg)))
    R["power"] = power(y, scores["model"], scores["ss_only"], fam)
    R["artifact_diagnostics"] = diagnostics_block(y, scores["model"], fam, rows, F)
    R["subsets"] = subsets(y, scores["model"], fam, z, F)
    return R


def check_frozen():
    dirty = subprocess.run(["git", "status", "--porcelain", "--", *FROZEN_FILES],
                           capture_output=True, text=True).stdout.strip()
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", "PROTOCOL.md"],
                             capture_output=True).returncode == 0
    if dirty or not tracked:
        raise SystemExit(f"refusing: PROTOCOL.md not committed or frozen files changed:\n{dirty}")
    if os.path.exists(OUT["lockbox"]):
        raise SystemExit(f"refusing: {OUT['lockbox']} exists -- the lockbox is evaluated once")
    return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()


def run_lockbox():
    commit = check_frozen()
    rows_d, _, F_d = load("dev")
    rows_l, null_l, F_l = load("lockbox")
    D, L = arrays(rows_d), arrays(rows_l)
    assert not set(D["fam"]) & set(L["fam"]), "family crosses dev/lockbox"
    assert not set(D["prot"]) & set(L["prot"]), "protein crosses dev/lockbox"
    y, fam, prot, z = L["y"], L["fam"], L["prot"], L["z"]
    X_d, X_l = X_of(F_d, FEATURES), X_of(F_l, FEATURES)
    scores = {"model": fit(X_d, D["y"], D["w"]).predict_proba(X_l)[:, 1]}
    for b, cols in BASELINES.items():
        scores[b] = fit(X_of(F_d, cols), D["y"], D["w"]).predict_proba(X_of(F_l, cols))[:, 1]
    R = dict(mode="lockbox", commit=commit, n=len(y), movers=int(y.sum()),
             proteins=len(set(prot)), families=len(set(fam)), null=null_fp(null_l),
             exploratory=not (reportable(y) and len(set(fam)) >= MIN_FAMILIES))
    R["headline"] = headline(y, scores, fam, prot) if 0 < y.sum() < len(y) else None
    nt = D["prot"] != "T4L"
    s_not4l = fit(X_d[nt], D["y"][nt], D["w"][nt]).predict_proba(X_l)[:, 1]
    R["t4l"] = dict(trained_without_t4l=auc_ci(y, s_not4l, fam),
                    with_vs_without=paired(y, s_not4l, scores["model"], fam))
    # positive control: helix, model trained on dev, scored on lockbox
    yh_d = X_of(F_d, ["ss_H"])[:, 0] == 1; yh_l = X_of(F_l, ["ss_H"])[:, 0] == 1
    ph = fit(X_of(F_d, POSITIVE_FEATURES), yh_d).predict_proba(X_of(F_l, POSITIVE_FEATURES))[:, 1]
    R["positive_control_helix"] = auc_ci(yh_l, ph, fam)
    # negative control: lockbox labels shuffled within families, frozen scores
    rng = np.random.default_rng(0)
    neg = [fast_auc(y[shuffle_within(fam, rng)], scores["model"]) for _ in range(N_SHUFFLE_LOCKBOX)]
    R["negative_control_shuffled_within_family"] = dict(
        reps=N_SHUFFLE_LOCKBOX, mean_auc=float(np.mean(neg)), sd=float(np.std(neg)),
        q025=float(np.percentile(neg, 2.5)), q975=float(np.percentile(neg, 97.5)))
    R["power"] = power(y, scores["model"], scores["ss_only"], fam)
    R["artifact_diagnostics"] = diagnostics_block(y, scores["model"], fam, rows_l, F_l)
    R["subsets"] = subsets(y, scores["model"], fam, z, F_l)
    return R


def show(R):
    h = R.get("headline")
    print(f"[{R['mode']}] n={R['n']} movers={R['movers']} null FP={R['null']['fp']}")
    if h:
        m = h["model"]
        print(f"  model AUC {m['auc']:.3f} [{m['lo']:.3f},{m['hi']:.3f}] ({m['families']} families), "
              f"ceiling {h['ceiling_oracle_auc_dev']}")
        for b in BASELINES:
            d = h[f"model_vs_{b}"]
            print(f"  {b:12s} AUC {h[b]['auc']:.3f} [{h[b]['lo']:.3f},{h[b]['hi']:.3f}]  "
                  f"model - {b}: {d['d_auc']:+.3f} [{d['lo']:+.3f},{d['hi']:+.3f}] p={d['p_two_sided']:.3g}")
        wp = h["within_protein"]
        print(f"  within-protein AUC {wp['auc']:.3f} [{wp['lo']:.3f},{wp['hi']:.3f}]")
    if R.get("exploratory"):
        print("  LOCKBOX TOO SMALL (< %d rows, < %d per class or < %d families): all results exploratory"
              % (MIN_GROUP, MIN_CLASS, MIN_FAMILIES))
    print(f"  T4L: {json.dumps(R['t4l'], default=float)[:400]}")
    pc = R["positive_control_helix"]; nc = R["negative_control_shuffled_within_family"]
    print(f"  positive control (helix) AUC {pc['auc']:.3f} [{pc['lo']:.3f},{pc['hi']:.3f}]")
    print(f"  negative control (shuffled within family) AUC {nc['mean_auc']:.3f} +- {nc['sd']:.3f}")
    p = R["power"]
    print(f"  power: {p['families']} families, n_eff {p['n_eff_families']:.1f}, SE {p['se_auc']:.3f} "
          f"(design effect {p['design_effect']:.1f}); min detectable AUC {p['min_detectable_auc']:.3f}, "
          f"min detectable dAUC vs SS {p['min_detectable_d_auc']:.3f}")
    for k, v in R["artifact_diagnostics"].items():
        e = v["excluded"]
        tail = f"AUC {e['auc']:.3f} [{e['lo']:.3f},{e['hi']:.3f}]" if "auc" in e else e["note"]
        print(f"  excl. {k:32s} suspect n={v['n_suspect']:4d}  -> {tail}")
    for g in R["subsets"]:
        tail = f"AUC {g['auc']:.3f} [{g['lo']:.3f},{g['hi']:.3f}]" if "auc" in g else g["note"]
        print(f"  subset {g['name']:32s} n={g['n']:4d} {tail}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("dev", "lockbox"))
    a = ap.parse_args()
    R = run_dev() if a.mode == "dev" else run_lockbox()
    os.makedirs("results", exist_ok=True)
    with open(OUT[a.mode], "w") as fh:
        json.dump(R, fh, indent=1, default=float)
    show(R)
    print(f"wrote {OUT[a.mode]}")


if __name__ == "__main__":
    main()
