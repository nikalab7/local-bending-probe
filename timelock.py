"""
proteinX -- the pre-registered time-based lockbox (TIMELOCK.md).

The frozen pipeline of PROTOCOL.md is applied to PDB entries released after
CUTOFF. Nothing in the label, features, model or analysis changes; this file
only selects new data, gates unblinding on power, and runs the evaluation once.

  python timelock.py freeze-dev
      (pre-registration only) writes results/frozen_dev_matrix.npz: the dev
      training rows (z, family, protein, SITE + SUBST features) the frozen
      model is trained on, so a later run needs no dev rebuild.
  python timelock.py feasibility
      historical rate of newly eligible families -> results/timelock_feasibility.json
  python timelock.py build
      stratum A: proteins not in dev or the first lockbox, no >= 30% hit to any
                 of them, with >= 1 single-substitution crystal released after
                 CUTOFF (WT crystals of any date)                 -> manifests/timelock_A.json
      stratum B: dev / first-lockbox proteins, single-substitution crystals
                 released after CUTOFF of mutations not in either manifest
                                                                   -> manifests/timelock_B.json
  python timelock.py gate --stratum A|B
      BLINDED: builds the labels but reports only counts (rows, families,
      Kish n_eff, WT-vs-WT null FP) and the predicted minimum detectable AUC.
      Pass if predicted MDA <= GATE_MDA. Never prints mover labels or scores.
  python timelock.py eval --stratum A|B
      Once per stratum, only after a passed gate, in the pinned environment,
      with frozen files identical to the pre-registration commit.
"""
from __future__ import annotations
import argparse
import collections
import concurrent.futures as cf
import hashlib
import json
import os
import subprocess
import numpy as np

CUTOFF = "2026-10-05"            # entries released strictly after this date
GATE_MDA = 0.62                  # unblind only if predicted min. detectable AUC <= this
ASSUMED_MOVER_RATE = 0.30        # for the blinded power prediction (dev 0.36, lockbox 0.27)
DESIGN_EFFECT = {"A": 1.1, "B": 0.9}   # measured: first lockbox (A-like), dev (B-like)
NULL_FP_WARN = 0.04              # above this, the n_wt >= 5 secondary analysis is reported
SECONDARY_MIN_WT = 5
DEV_MATRIX = os.path.join("results", "frozen_dev_matrix.npz")
TL_DIR = "timelock_pdb"
FROZEN = ("PROTOCOL.md", "TIMELOCK.md", "final_eval.py", "pairs.py", "delta_model.py",
          "structure_features.py", "plm_features.py", "stats_utils.py", "mine_pairs.py",
          "timelock.py", "requirements.lock", "manifests/mined_ms4.json",
          "manifests/lockbox.json", DEV_MATRIX)
CORE_PACKAGES = ("numpy", "scipy", "scikit-learn", "gemmi")
GQL_DATE = """query($ids:[String!]!){ entries(entry_ids:$ids){
  rcsb_id rcsb_accession_info{ initial_release_date } } }"""


def manifest_path(stratum):
    return os.path.join("manifests", f"timelock_{stratum}.json")


# ------------------------------ pre-registration ----------------------------
def freeze_dev():
    import final_eval as fe
    import delta_model as dm
    rows, _, F = fe.load("dev")
    rows_l, _, _ = fe.load("lockbox")
    cols = dm.SITE + dm.SUBST
    np.savez_compressed(DEV_MATRIX, z=np.array([x["z"] for x in rows], float),
                        family=np.array([str(x["family"]) for x in rows]),
                        protein=np.array([str(x["protein"]) for x in rows]),
                        cols=np.array(cols), X=dm.matrix(F, cols),
                        known_mutations=np.array(sorted({mutation_id(x) for x in rows + rows_l})))
    print(f"wrote {DEV_MATRIX} ({len(rows)} rows) sha256 {sha256(DEV_MATRIX)}")


def mutation_id(x):
    return f"{x['protein']}:{x['wt']}{x['r']}{x['mut']}"


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def release_dates(pdb_ids):
    import mine_pairs as mp
    ids, out = sorted(set(pdb_ids)), {}
    for i in range(0, len(ids), 300):
        r = mp._post(mp.GRAPHQL, {"query": GQL_DATE, "variables": {"ids": ids[i:i + 300]}})
        for e in r["data"]["entries"] or []:
            if e:
                out[e["rcsb_id"]] = e["rcsb_accession_info"]["initial_release_date"][:10]
    return out


def eligible_date(entries, dates):
    """Earliest date a protein had a form with >= 3 WT crystals and >= 1 single mutant."""
    forms = collections.defaultdict(lambda: dict(wt=[], single=[]))
    for p, e in entries.items():
        forms[e["form"]][e["kind"]].append(dates[p])
    d = [max(sorted(v["wt"])[2], min(v["single"])) for v in forms.values()
         if len(v["wt"]) >= 3 and v["single"]]
    return min(d) if d else None


def feasibility():
    """New families per year (first-eligible year of each dev / lockbox family)."""
    import final_eval as fe
    mans = {t: json.load(open(p)) for t, p in (("dev", "manifests/mined_ms4.json"),
                                              ("lockbox", "manifests/lockbox.json"))}
    dates = release_dates([p for m in mans.values() for v in m["proteins"].values() for p in v["entries"]])
    first = {}
    for t, m in mans.items():
        for v in m["proteins"].values():
            d = eligible_date(v["entries"], dates)
            if d:
                k = (t, str(v["family"]))
                first[k] = min(first.get(k, "9999"), d)
    per_year = collections.Counter(d[:4] for d in first.values())
    recent = [per_year.get(str(y), 0) for y in range(2016, 2026)]
    rows_l, _, _ = fe.load("lockbox")
    lab_fam = len({str(x["family"]) for x in rows_l})
    lb_fam = len({v["family"] for v in mans["lockbox"]["proteins"].values()})
    yield_fam = lab_fam / lb_fam
    rows_per_fam = len(rows_l) / lab_fam
    rate_post = float(np.mean(recent)) * yield_fam

    def families_needed(target_auc, deff=DESIGN_EFFECT["A"], q=ASSUMED_MOVER_RATE):
        se = (target_auc - 0.5) / 2.8016
        n = deff / (12 * q * (1 - q) * se ** 2)          # Hanley-McNeil at AUC 0.5, x design effect
        return n, n / rows_per_fam
    need = {str(a): dict(zip(("rows", "families"), map(float, families_needed(a)))) for a in (0.62, 0.66)}
    out = dict(cutoff=CUTOFF, families_per_first_eligible_year=dict(sorted(per_year.items())),
               mean_new_families_per_year_2016_2025=float(np.mean(recent)),
               post_qc_family_yield_first_lockbox=yield_fam, rows_per_family_first_lockbox=rows_per_fam,
               expected_post_qc_new_families_per_year=rate_post, needed=need,
               years_needed={k: v["families"] / rate_post for k, v in need.items()})
    with open(os.path.join("results", "timelock_feasibility.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "families_per_first_eligible_year"}, indent=1))


# ------------------------------ build ---------------------------------------
def _screen(acc):
    import mine_pairs as mp
    recs = [r for r in map(mp._flatten, mp.entity_meta(mp.entity_ids(acc))) if r]
    return mp.prescreen([r for r in recs if r["resolution"] <= 2.5])


def _entries(ps, keep_single):
    entries = {}
    for f, v in ps["forms"].items():
        singles = [r for r in v["single"] if keep_single(r)]
        if not singles:
            continue
        for kind, rs in (("wt", v["wt"]), ("single", singles)):
            for r in rs:
                entries[r["pdb"]] = dict(chain=r["chain"], form=f, kind=kind, resolution=r["resolution"])
    return entries


def build(workers=8):
    import mine_pairs as mp
    old = {t: json.load(open(p)) for t, p in (("dev", "manifests/mined_ms4.json"),
                                             ("lockbox", "manifests/lockbox.json"))}
    known_acc = {a for m in old.values() for a in m["proteins"]} | set(old["dev"].get("covered", {}))
    known_seq = {v["reference"] for m in old.values() for v in m["proteins"].values()}
    cands = [a for a, _ in mp.candidate_accessions(1, 20000)]

    def screen(acc):
        try:
            return acc, _screen(acc)
        except Exception:
            return acc, None
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        ps = {a: p for a, p in ex.map(screen, cands) if p is not None}
    pdbs = {r["pdb"] for p in ps.values() for v in p["forms"].values() for r in v["single"] + v["wt"]}
    dates = release_dates(pdbs)
    new = lambda r: dates.get(r["pdb"], "0000") > CUTOFF

    # stratum B: new single-mutant crystals of known proteins; rows of mutations
    # already labelled in dev or the first lockbox are dropped when labels are built
    def seen_mutations(acc):
        for m in old.values():
            v = m["proteins"].get(acc)
            if v:
                return v
        return None
    B = {}
    for acc, p in ps.items():
        v = seen_mutations(acc)
        if v is None or p["reference"] != v["reference"]:
            continue
        e = _entries(p, new)
        if any(x["kind"] == "single" for x in e.values()):
            B[acc] = dict(family=str(v["family"]), reference=p["reference"], entries=e)
    # stratum A: new proteins, not homologous to anything known
    A0 = {a: p for a, p in ps.items() if a not in known_acc and p["reference"] not in known_seq}
    A0 = {a: p for a, p in A0.items() if any(x["kind"] == "single" for x in _entries(p, new).values())}

    def homolog(acc):
        try:
            for m in mp.entity_meta(mp.family_hits(A0[acc]["reference"])):
                accs = {r["database_accession"] for r in
                        (m["rcsb_polymer_entity_container_identifiers"]["reference_sequence_identifiers"] or [])
                        if r["database_name"] == "UniProt"}
                seq = (m.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can")
                if accs & known_acc or seq in known_seq:
                    return acc, True
            return acc, False
        except Exception:
            return acc, True                      # cannot verify -> exclude
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        hom = dict(ex.map(homolog, sorted(A0)))
    A1 = {a: p for a, p in A0.items() if not hom[a]}
    fam = mp.families({a: p["reference"] for a, p in A1.items()}, workers) if A1 else {}
    A = {a: dict(family="T" + fam[a], reference=p["reference"], entries=_entries(p, new))
         for a, p in A1.items()}
    for s, man in (("A", A), ("B", B)):
        with open(manifest_path(s), "w") as fh:
            json.dump(dict(cutoff=CUTOFF, stratum=s, proteins=man), fh, indent=1, sort_keys=True)
        print(f"stratum {s}: {len(man)} proteins, "
              f"{sum(x['kind'] == 'single' for v in man.values() for x in v['entries'].values())} "
              f"new single-mutant entries -> {manifest_path(s)}")
    os.makedirs(TL_DIR, exist_ok=True)
    allp = sorted({p for man in (A, B) for v in man.values() for p in v["entries"]})
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(lambda p: mp.fetch_gz(p, TL_DIR), allp))
    with open(os.path.join("manifests", "timelock_files.json"), "w") as fh:
        json.dump({p: sha256(os.path.join(TL_DIR, f"{p}.pdb.gz")) for p in allp
                   if os.path.exists(os.path.join(TL_DIR, f"{p}.pdb.gz"))}, fh, indent=0, sort_keys=True)


# ------------------------------ gate / eval ---------------------------------
def _point_lockbox_at(stratum):
    """Route final_eval's lockbox loader to a timelock stratum (no code changes).

    Stratum B (known proteins): rows of mutations already labelled in dev or the
    first lockbox are dropped, and family / protein IDs get a "B:" prefix so the
    frozen disjointness assertions hold; the family bootstrap still resamples
    whole families. The overlap with training proteins is the design of B.
    """
    import mine_pairs as mp
    import final_eval as fe
    mp.LOCKBOX_MANIFEST, mp.LOCKBOX_DIR = manifest_path(stratum), TL_DIR
    fe.CACHE = os.path.join("results", "_cache", f"timelock_{stratum}")
    fe.OUT["lockbox"] = os.path.join("results", f"timelock_{stratum}.json")
    real_load = fe.load
    known = set(np.load(DEV_MATRIX)["known_mutations"].tolist())

    def load(which):
        if which == "dev" and os.path.exists(DEV_MATRIX):
            return _frozen_dev()
        rows, null, F = real_load(which)
        if which == "lockbox" and stratum == "B":
            keep = [i for i, x in enumerate(rows) if mutation_id(x) not in known]
            rows = [dict(rows[i], family="B:" + str(rows[i]["family"]),
                         protein="B:" + str(rows[i]["protein"])) for i in keep]
            F = [F[i] for i in keep]
        return rows, null, F
    fe.load = load
    return fe


def _frozen_dev():
    """final_eval.load('dev') replacement: the committed frozen dev matrix."""
    D = np.load(DEV_MATRIX, allow_pickle=False)
    cols = [str(c) for c in D["cols"]]
    rows = [dict(z=float(z), family=str(f), protein=str(p))
            for z, f, p in zip(D["z"], D["family"], D["protein"])]
    F = [dict(zip(cols, map(float, x))) for x in D["X"]]
    return rows, [], F


def predicted_mda(n, stratum, q=ASSUMED_MOVER_RATE):
    if n < 2:
        return float("inf"), float("inf")
    n1, n0 = q * n, (1 - q) * n
    se = np.sqrt(DESIGN_EFFECT[stratum] * (n + 1) / (12 * n1 * n0))
    return float(0.5 + 2.8016 * se), float(se)


def gate(stratum):
    fe = _point_lockbox_at(stratum)
    rows, null, _ = fe.load("lockbox")
    fam = np.array([str(x["family"]) for x in rows])
    n_f = np.array([np.sum(fam == f) for f in np.unique(fam)], float)
    mda, se = predicted_mda(len(rows), stratum)
    zn = np.abs([x["z"] for x in null])
    nwt = np.array([x["n_wt"] for x in rows])
    out = dict(stratum=stratum, cutoff=CUTOFF, n_rows=len(rows), families=int(len(n_f)),
               n_eff_families=float(n_f.sum() ** 2 / (n_f ** 2).sum()) if len(n_f) else 0.0,
               rows_n_wt_ge_5=int((nwt >= SECONDARY_MIN_WT).sum()),
               null_n=int(len(zn)), null_fp=float((zn > 2.0).mean()) if len(zn) else None,
               predicted_se=se, predicted_mda=mda, gate_mda=GATE_MDA, passed=bool(mda <= GATE_MDA))
    with open(os.path.join("results", f"timelock_gate_{stratum}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))          # counts only: no labels, no scores


def prereg_commit():
    """The pre-registration commit: the one that added TIMELOCK.md."""
    c = subprocess.run(["git", "log", "--diff-filter=A", "--format=%H", "--", "TIMELOCK.md"],
                       capture_output=True, text=True).stdout.split()
    if not c:
        raise SystemExit("refusing: TIMELOCK.md is not committed (no pre-registration commit)")
    return c[-1]


def check_pinned():
    import importlib.metadata as md
    pin = {l.split("==")[0].lower(): l.split("==")[1].strip() for l in open("requirements.lock")
           if "==" in l}
    bad = [f"{p} {md.version(p)} != {pin[p]}" for p in CORE_PACKAGES if md.version(p) != pin[p]]
    c = prereg_commit()
    diff = subprocess.run(["git", "diff", "--name-only", c, "--", *FROZEN],
                          capture_output=True, text=True).stdout.split()
    if bad or diff:
        raise SystemExit(f"refusing: environment {bad} / files changed since {c[:8]}: {diff}")
    return c


def evaluate(stratum):
    commit = check_pinned()
    g = json.load(open(os.path.join("results", f"timelock_gate_{stratum}.json")))
    if not g["passed"]:
        raise SystemExit(f"refusing: gate not passed (predicted MDA {g['predicted_mda']:.3f} > {GATE_MDA})")
    fe = _point_lockbox_at(stratum)
    if os.path.exists(fe.OUT["lockbox"]):
        raise SystemExit(f"refusing: {fe.OUT['lockbox']} exists -- each stratum is evaluated once")
    fe.check_frozen = lambda: commit
    R = fe.run_lockbox()
    R.update(stratum=stratum, prereg_commit=commit, cutoff=CUTOFF, gate=g)
    if g["null_fp"] is not None and g["null_fp"] > NULL_FP_WARN:
        R["secondary_n_wt_ge_5"] = _secondary(fe, stratum)
    with open(fe.OUT["lockbox"], "w") as fh:
        json.dump(R, fh, indent=1, default=float)
    fe.show(R)


def _secondary(fe, stratum):
    """Pre-registered secondary: rows whose WT floor has >= 5 crystals."""
    rows_d, _, F_d = _frozen_dev()
    rows_l, _, F_l = fe.load("lockbox")
    keep = [i for i, x in enumerate(rows_l) if x["n_wt"] >= SECONDARY_MIN_WT]
    D, L = fe.arrays(rows_d), fe.arrays([rows_l[i] for i in keep])
    m = fe.fit(fe.X_of(F_d, fe.FEATURES), D["y"], D["w"])
    s = m.predict_proba(fe.X_of([F_l[i] for i in keep], fe.FEATURES))[:, 1]
    return fe.group(L["y"], s, L["fam"], np.ones(len(keep), bool), f"n_wt >= {SECONDARY_MIN_WT}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("freeze-dev", "feasibility", "build", "gate", "eval"))
    ap.add_argument("--stratum", choices=("A", "B"))
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    if a.cmd in ("gate", "eval") and not a.stratum:
        ap.error("--stratum is required")
    {"freeze-dev": freeze_dev, "feasibility": feasibility,
     "build": lambda: build(a.workers), "gate": lambda: gate(a.stratum),
     "eval": lambda: evaluate(a.stratum)}[a.cmd]()


if __name__ == "__main__":
    main()
