"""
proteinX -- is the mutant-vs-WT excess experimental? (LABTEST.md)

Exploratory; design, definitions and decision rules were committed in
LABTEST.md before any metadata was fetched. The frozen label is unchanged;
only NULL_HELDOUT is raised (more WT-WT pseudo-mutants) and the reference WT
sets are captured by wrapping pairs.null_rows at run time.

  python labtest.py labels      -> results/_cache/labtest_{dev,lockbox}.pkl
  python labtest.py metadata    -> results/lab_metadata.json.gz
  python labtest.py analyse     -> results/labtest.json
  python labtest.py redo        -> results/labtest_redo.json
"""
from __future__ import annotations
import argparse
import concurrent.futures as cf
import gzip
import json
import os
import pickle
import re
import urllib.request
import numpy as np

CACHE = os.path.join("results", "_cache")
META = os.path.join("results", "lab_metadata.json.gz")
NULL_HELDOUT_TEST = 40
PDB_DIRS = ("t4l_pdb", "val_pdb", "mined_pdb", "lockbox_pdb")
REDO_DIR, HYBRID_DIR = "redo_pdb", "redo_hybrid"
REDO_N, REDO_MAX_ENTRIES, REDO_MIN_COVER, SEED = 40, 150, 0.9, 0
N_BOOT = 2000
CONSORTIUM = re.compile(r"project|consortium|cent(er|re)|initiative|genomics|institute|laboratory|university", re.I)
PROGRAMS = ("REFMAC", "PHENIX", "BUSTER", "CNS", "X-PLOR", "XPLOR", "SHELX", "TNT", "PROLSQ", "RESTRAIN")
GQL = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
  audit_author{ name } software{ name classification }
  rcsb_accession_info{ deposit_date } pdbx_database_status{ status_code_sf } } }"""


# ------------------------------ labels --------------------------------------
def build_labels(which, null_heldout=NULL_HELDOUT_TEST):
    import pairs
    import diagnostics
    import final_eval as fe
    fe.freeze_label()
    pairs.NULL_HELDOUT = null_heldout
    diagnostics.memoize_parsing()
    refsets, cryst = {}, {}
    orig = pairs.null_rows

    def wrapped(rows, structs, wt_by_form, cons, prior):
        if rows:
            prot = rows[0]["protein"]
            for f, w in wt_by_form.items():
                refsets[(prot, f)] = sorted(w)
        for p, st in structs.items():
            cryst[p] = (st.get("temperature"), st.get("resolution"))
        return orig(rows, structs, wt_by_form, cons, prior)
    pairs.null_rows = wrapped
    null = []
    data = pairs.load_proteins(null_out=null) if which == "dev" else pairs.load_lockbox(null_out=null)
    pairs.null_rows = orig
    rows = [x for rs, _ in data.values() for x in rs]
    return rows, null, refsets, cryst


def labels():
    import final_eval as fe
    os.makedirs(CACHE, exist_ok=True)
    for which in ("dev", "lockbox"):
        rows, null, refsets, cryst = build_labels(which)
        frozen, _, _ = fe.load(which)
        key = lambda xs: sorted((x["protein"], x["r"], x["wt"], x["mut"], round(x["z"], 6)) for x in xs)
        if key(rows) != key(frozen):
            raise SystemExit(f"{which}: real rows differ from the frozen-label rows -- stopping")
        slim = lambda x: {k: v for k, v in x.items() if k not in ("scaffold_res", "scaffold_helix", "scaffold_sheet")}
        with open(os.path.join(CACHE, f"labtest_{which}.pkl"), "wb") as fh:
            pickle.dump(dict(rows=[slim(x) for x in rows], null=[slim(x) for x in null],
                             refsets=refsets, cryst=cryst), fh, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"{which}: {len(rows)} rows (identical to frozen), {len(null)} null pseudo-mutants, "
              f"{len(refsets)} reference sets")


# ------------------------------ metadata ------------------------------------
def all_entry_ids():
    ids = set()
    for d in PDB_DIRS:
        if os.path.isdir(d):
            ids |= {f.split(".")[0].upper() for f in os.listdir(d)}
    return sorted(ids)


def metadata(workers=6):
    import mine_pairs as mp
    ids = all_entry_ids()
    chunks = [ids[i:i + 200] for i in range(0, len(ids), 200)]

    def fetch(ch):
        r = mp._post(mp.GRAPHQL, {"query": GQL, "variables": {"ids": ch}})
        return [e for e in (r.get("data") or {}).get("entries") or [] if e]
    out = {}
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        for es in ex.map(fetch, chunks):
            for e in es:
                out[e["rcsb_id"]] = dict(
                    authors=[a["name"] for a in e.get("audit_author") or []],
                    refinement=[s["name"] for s in e.get("software") or []
                                if (s.get("classification") or "").lower() == "refinement" and s.get("name")],
                    deposit=(e.get("rcsb_accession_info") or {}).get("deposit_date"),
                    sf=(e.get("pdbx_database_status") or {}).get("status_code_sf"))
    with gzip.open(META, "wt") as fh:
        json.dump(out, fh)
    print(f"metadata for {len(out)}/{len(ids)} entries -> {META}")


def load_meta():
    with gzip.open(META, "rt") as fh:
        return json.load(fh)


def author_keys(names):
    keys = set()
    for n in names:
        if "," not in n or CONSORTIUM.search(n):
            continue
        last, first = [t.strip() for t in n.split(",", 1)]
        if last:
            keys.add(f"{last.lower()},{first[:1].lower()}")
    return keys


def program(names):
    if not names:
        return None
    n = names[0].upper()
    for p in PROGRAMS:
        if n.startswith(p):
            return "XPLOR" if p == "X-PLOR" else p
    return "OTHER:" + n.split()[0]


def year(d):
    return int(d[:4]) if d else None


# ------------------------------ relations -----------------------------------
def relation(e_ids, R_ids, meta, cryst):
    """Lab / program / year / condition relation of entity e to reference set R."""
    e_auth = set().union(*[author_keys(meta.get(p, {}).get("authors", [])) for p in e_ids])
    R = [p for p in R_ids if p in meta]
    if not R or not e_auth:
        return None
    f_lab = float(np.mean([bool(e_auth & author_keys(meta[p]["authors"])) for p in R]))
    e_prog = {program(meta.get(p, {}).get("refinement", [])) for p in e_ids} - {None}
    R_prog = [program(meta[p]["refinement"]) for p in R]
    f_prog = float(np.mean([q in e_prog for q in R_prog])) if e_prog and any(R_prog) else np.nan
    ey = [year(meta.get(p, {}).get("deposit")) for p in e_ids]
    Ry = [year(meta[p]["deposit"]) for p in R]
    ey = [v for v in ey if v]; Ry = [v for v in Ry if v]
    gap = abs(np.median(ey) - np.median(Ry)) if ey and Ry else np.nan

    def med(ids, k):
        v = [cryst.get(p, (None, None))[k] for p in ids]
        v = [x for x in v if x is not None and np.isfinite(x)]
        return float(np.median(v)) if v else np.nan
    dres = med(e_ids, 1) - med(R, 1)
    dT = abs(med(e_ids, 0) - med(R, 0))
    return dict(f_lab=f_lab, f_prog=f_prog, year_gap=float(gap), dres=float(dres), dT=float(dT))


def lab_class(f):
    return "same" if f >= 0.5 else ("cross" if f == 0 else "mixed")


def annotate(D, meta):
    out = []
    for kind, xs in (("real", D["rows"]), ("null", D["null"])):
        for x in xs:
            R = D["refsets"].get((x["protein"], x["form"]))
            if R is None:
                continue
            e = list(x["pdbs"])
            rel = relation(e, [p for p in R if p not in e], meta, D["cryst"])
            if rel is None:
                continue
            out.append(dict(kind=kind, protein=str(x["protein"]), family=str(x["family"]), form=x["form"],
                            y=bool(abs(x["z"]) > 2.0), wt=x["wt"], mut=x["mut"], r=x["r"], s=x["s"],
                            heldout=x.get("heldout"), **rel))
    return out


# ------------------------------ statistics ----------------------------------
def _boot_idx(fam, seed=SEED):
    fams = np.unique(fam); by = {f: np.flatnonzero(fam == f) for f in fams}
    rng = np.random.default_rng(seed)
    for _ in range(N_BOOT):
        yield np.concatenate([by[f] for f in rng.choice(fams, len(fams))])


def strat_diff(y, a, b, strata):
    """CMH-type risk difference rate(a) - rate(b) over strata containing both."""
    num = den = 0.0
    for s in np.unique(strata[a | b]):
        m = strata == s
        na, nb = (a & m).sum(), (b & m).sum()
        if na and nb:
            w = na * nb / (na + nb)
            num += w * (y[a & m].mean() - y[b & m].mean()); den += w
    return num / den if den else np.nan


def strat_ci(y, a, b, strata, fam):
    est = strat_diff(y, a, b, strata)
    bs = [strat_diff(y[i], a[i], b[i], strata[i]) for i in _boot_idx(fam)]
    bs = [v for v in bs if np.isfinite(v)]
    n_strata = sum(1 for s in np.unique(strata) if (a & (strata == s)).any() and (b & (strata == s)).any())
    keep = np.isin(strata, [s for s in np.unique(strata)
                            if (a & (strata == s)).any() and (b & (strata == s)).any()])
    return dict(diff=float(est), ci=[float(v) for v in np.percentile(bs, [2.5, 97.5])] if bs else None,
                n_strata=int(n_strata), n_a=int((a & keep).sum()), n_b=int((b & keep).sum()),
                families=int(len(set(fam[keep]))))


def rate_ci(y, m, fam):
    if not m.any():
        return dict(n=0)
    bs = [y[i][m[i]].mean() for i in _boot_idx(fam) if m[i].any()]
    return dict(n=int(m.sum()), families=int(len(set(fam[m]))), rate=float(y[m].mean()),
                ci=[float(v) for v in np.percentile(bs, [2.5, 97.5])])


def excess_ci(y, real, null, m, fam):
    """(real mover rate - null mover rate) within mask m, family bootstrap."""
    a, b = real & m, null & m
    if not a.any() or not b.any():
        return dict(n_real=int(a.sum()), n_null=int(b.sum()))
    bs = [y[i][a[i]].mean() - y[i][b[i]].mean() for i in _boot_idx(fam) if a[i].any() and b[i].any()]
    return dict(n_real=int(a.sum()), n_null=int(b.sum()), real=float(y[a].mean()), null=float(y[b].mean()),
                excess=float(y[a].mean() - y[b].mean()), ci=[float(v) for v in np.percentile(bs, [2.5, 97.5])])


def decide(est, cross_rate):
    if est["ci"] is None or cross_rate is None:
        return "inconclusive (no estimate)"
    lo, hi = est["ci"]
    if est["diff"] >= 0.05 and lo > 0 and cross_rate >= 0.10:
        return "experimental excess SUPPORTED"
    if hi < 0.05 and cross_rate < 0.10:
        return "experimental excess REJECTED"
    return "INCONCLUSIVE"


# ------------------------------ analysis ------------------------------------
def analyse():
    meta = load_meta()
    A = []
    for which in ("dev", "lockbox"):
        with open(os.path.join(CACHE, f"labtest_{which}.pkl"), "rb") as fh:
            A += annotate(pickle.load(fh), meta)
    g = lambda k: np.array([a[k] for a in A])
    y, kind, fam = g("y"), g("kind"), g("family")
    strata = np.array([f"{a['protein']}|{a['form']}" for a in A])
    f_lab, f_prog, gap, dres, dT = g("f_lab"), g("f_prog").astype(float), g("year_gap").astype(float), \
        g("dres").astype(float), g("dT").astype(float)
    with np.errstate(invalid="ignore"):
        matched = (np.abs(dres) <= 0.3) & (dT <= 50)
    cross, same, mixed = f_lab == 0, f_lab >= 0.5, (f_lab > 0) & (f_lab < 0.5)
    real, null = kind == "real", kind == "null"
    R = dict(n_real=int(real.sum()), n_null=int(null.sum()), primary={}, secondary={})
    # ---- primary: WT-WT null
    P = R["primary"]
    for nm, m in (("cross", cross), ("same", same), ("mixed", mixed)):
        P[f"rate_{nm}"] = rate_ci(y, null & m, fam)
        P[f"rate_{nm}_matched"] = rate_ci(y, null & m & matched, fam)
    P["stratified_cross_minus_same"] = strat_ci(y, null & cross, null & same, strata, fam)
    P["stratified_cross_minus_same_matched"] = strat_ci(y, null & cross & matched, null & same & matched, strata, fam)
    cr = P["rate_cross_matched"].get("rate")
    P["decision"] = decide(P["stratified_cross_minus_same_matched"], cr)
    # ---- secondary: mutant-WT rows
    S = R["secondary"]
    for nm, m in (("cross", cross), ("same", same), ("mixed", mixed)):
        S[f"rate_{nm}"] = rate_ci(y, real & m, fam)
    S["stratified_cross_minus_same"] = strat_ci(y, real & cross, real & same, strata, fam)
    S["stratified_cross_minus_same_matched"] = strat_ci(y, real & cross & matched, real & same & matched, strata, fam)
    with np.errstate(invalid="ignore"):
        dprog, sprog = f_prog == 0, f_prog >= 0.5
    S["program_different_minus_same"] = strat_ci(y, real & dprog, real & sprog, strata, fam)
    S["program_different_minus_same_null"] = strat_ci(y, null & dprog, null & sprog, strata, fam)
    S["year_gap_bins"] = {}
    for nm, (lo, hi) in (("0-2", (0, 2)), ("3-9", (3, 9)), (">=10", (10, 1e9))):
        with np.errstate(invalid="ignore"):
            m = (gap >= lo) & (gap <= hi)
        S["year_gap_bins"][nm] = dict(real=rate_ci(y, real & m, fam), null=rate_ci(y, null & m, fam))
    S["mutational_excess_same_lab"] = excess_ci(y, real, null, same, fam)
    S["mutational_excess_cross_lab"] = excess_ci(y, real, null, cross, fam)
    S["mutational_excess_same_lab_matched"] = excess_ci(y, real, null, same & matched, fam)
    S["mutational_excess_cross_lab_matched"] = excess_ci(y, real, null, cross & matched, fam)
    with open(os.path.join("results", "labtest.json"), "w") as fh:
        json.dump(R, fh, indent=1, default=float)
    show(R)


def _fmt_rate(d):
    return f"{d['rate']:.3f} [{d['ci'][0]:.3f},{d['ci'][1]:.3f}] n={d['n']} fam={d['families']}" if d.get("n") else "n=0"


def _fmt_diff(d):
    ci = f"[{d['ci'][0]:+.3f},{d['ci'][1]:+.3f}]" if d["ci"] else "[--]"
    return f"{d['diff']:+.3f} {ci} strata={d['n_strata']} n={d['n_a']}/{d['n_b']} fam={d['families']}"


def show(R):
    P, S = R["primary"], R["secondary"]
    print(f"real rows {R['n_real']}, null pseudo-mutants {R['n_null']}")
    print("PRIMARY (WT-WT null)")
    for k in ("cross", "same", "mixed"):
        print(f"  {k:6s} {_fmt_rate(P['rate_' + k])} | matched {_fmt_rate(P['rate_' + k + '_matched'])}")
    print(f"  within-protein cross - same          {_fmt_diff(P['stratified_cross_minus_same'])}")
    print(f"  within-protein cross - same, matched {_fmt_diff(P['stratified_cross_minus_same_matched'])}")
    print(f"  DECISION: {P['decision']}")
    print("SECONDARY (mutant-WT rows)")
    for k in ("cross", "same", "mixed"):
        print(f"  {k:6s} {_fmt_rate(S['rate_' + k])}")
    for k in ("stratified_cross_minus_same", "stratified_cross_minus_same_matched",
              "program_different_minus_same", "program_different_minus_same_null"):
        print(f"  {k:40s} {_fmt_diff(S[k])}")
    for b, v in S["year_gap_bins"].items():
        print(f"  year gap {b:5s} real {_fmt_rate(v['real'])} | null {_fmt_rate(v['null'])}")
    for k in ("same_lab", "cross_lab", "same_lab_matched", "cross_lab_matched"):
        e = S["mutational_excess_" + k]
        if "excess" in e:
            print(f"  mutational excess {k:18s} {e['real']:.3f} - {e['null']:.3f} = {e['excess']:+.3f} "
                  f"[{e['ci'][0]:+.3f},{e['ci'][1]:+.3f}] (n {e['n_real']}/{e['n_null']})")


# ------------------------------ PDB-REDO ------------------------------------
def redo_url(p):
    p = p.lower()
    return f"https://pdb-redo.eu/db/{p}/{p}_final.pdb"


def redo_available(p):
    try:
        req = urllib.request.Request(redo_url(p), method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status == 200
    except Exception:
        return False


def fetch_redo(p):
    path = os.path.join(REDO_DIR, f"{p.upper()}_final.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    for k in range(3):
        try:
            with urllib.request.urlopen(redo_url(p), timeout=120) as r:
                data = r.read()
            with open(path, "wb") as fh:
                fh.write(data)
            return path
        except Exception:
            pass
    return None


COORD = ("ATOM  ", "HETATM", "ANISOU", "TER", "MODEL ", "ENDMDL", "CONECT", "MASTER", "END")


def hybrid(orig_path, redo_path, out_path):
    """Original header records + PDB-REDO coordinate records."""
    op = gzip.open if orig_path.endswith(".gz") else open
    with op(orig_path, "rt") as fh:
        head = [l for l in fh if not l.startswith(COORD)]
    with open(redo_path) as fh:
        coords = [l for l in fh if l.startswith(("ATOM  ", "HETATM", "ANISOU", "TER", "MODEL ", "ENDMDL"))]
    with gzip.open(out_path, "wt") as fh:
        fh.writelines(head + coords + ["END\n"])
    return out_path


def redo(workers=6):
    import pairs
    import diagnostics
    import final_eval as fe
    import mine_pairs as mp
    fe.freeze_label(); pairs.NULL_HELDOUT = NULL_HELDOUT_TEST; diagnostics.memoize_parsing()
    mans = {}
    for path, d in (("manifests/mined_ms4.json", mp.MINED_DIR), ("manifests/lockbox.json", mp.LOCKBOX_DIR)):
        for acc, v in json.load(open(path))["proteins"].items():
            mans[acc] = (v, d)
    labelled = set()
    for which in ("dev", "lockbox"):
        with open(os.path.join(CACHE, f"labtest_{which}.pkl"), "rb") as fh:
            labelled |= {str(x["protein"]) for x in pickle.load(fh)["rows"]}
    elig = sorted(a for a, (v, _) in mans.items() if a in labelled and a != "T4L" and len(v["entries"]) <= REDO_MAX_ENTRIES)
    order = list(np.random.default_rng(SEED).permutation(elig))
    os.makedirs(REDO_DIR, exist_ok=True); os.makedirs(HYBRID_DIR, exist_ok=True)
    chosen, checked = [], 0
    for acc in order:
        if len(chosen) >= REDO_N:
            break
        ents = sorted(mans[acc][0]["entries"]); checked += 1
        with cf.ThreadPoolExecutor(max_workers=workers) as ex:
            avail = dict(zip(ents, ex.map(redo_available, ents)))
        cover = np.mean(list(avail.values()))
        if cover >= REDO_MIN_COVER:
            chosen.append((acc, avail))
    print(f"PDB-REDO subset: {len(chosen)} proteins (checked {checked} of {len(elig)} eligible)")
    rows_o, rows_r, null_o, null_r, n_ent, n_redo = [], [], [], [], 0, 0
    for acc, avail in chosen:
        v, d = mans[acc]
        chain = {p: e["chain"] for p, e in v["entries"].items()}
        orig = {p: os.path.join(d, f"{p}.pdb.gz") for p in v["entries"]}
        orig = {p: f for p, f in orig.items() if os.path.exists(f)}
        hyb = {}
        for p, f in orig.items():
            n_ent += 1
            rp = fetch_redo(p) if avail.get(p) else None
            if rp:
                hyb[p] = hybrid(f, rp, os.path.join(HYBRID_DIR, f"{p}.pdb.gz")); n_redo += 1
            else:
                hyb[p] = f
        for paths, rows_out, null_out in ((orig, rows_o, null_o), (hyb, rows_r, null_r)):
            nl = []
            rs, _ = pairs.build_pairs(paths, chain, acc, min_cons=3, family=v["family"], null_out=nl)
            rows_out += rs; null_out += nl
    key = lambda x: (str(x["protein"]), x["wt"], x["r"], x["mut"])
    ko = {key(x): x for x in rows_o}; kr = {key(x): x for x in rows_r}
    common = sorted(set(ko) & set(kr))
    nk = lambda x: (str(x["protein"]), x["s"], x["heldout"])
    no = {nk(x): x for x in null_o}; nr = {nk(x): x for x in null_r}
    ncommon = sorted(set(no) & set(nr))
    yo = np.array([abs(ko[k]["z"]) > 2 for k in common]); yr = np.array([abs(kr[k]["z"]) > 2 for k in common])
    fam = np.array([str(ko[k]["family"]) for k in common])
    nyo = np.array([abs(no[k]["z"]) > 2 for k in ncommon]); nyr = np.array([abs(nr[k]["z"]) > 2 for k in ncommon])
    nfam = np.array([str(no[k]["family"]) for k in ncommon])

    def paired(a, b, f):
        bs = [b[i].mean() - a[i].mean() for i in _boot_idx(f)]
        return dict(n=int(len(a)), families=int(len(set(f))), original=float(a.mean()), redo=float(b.mean()),
                    diff=float(b.mean() - a.mean()), ci=[float(v) for v in np.percentile(bs, [2.5, 97.5])])
    out = dict(proteins=len(chosen), entries=n_ent, entries_redo=n_redo,
               rows_original=len(rows_o), rows_redo=len(rows_r), rows_paired=len(common),
               mover_rate=paired(yo, yr, fam), null_fp=paired(nyo, nyr, nfam))
    # excess (real - null) before and after, family bootstrap on the union of families
    allf = np.concatenate([fam, nfam]); n1 = len(common)
    def exc(ix, real_y, null_y):
        i_r, i_n = ix[ix < n1], ix[ix >= n1] - n1
        return (real_y[i_r].mean() if len(i_r) else np.nan) - (null_y[i_n].mean() if len(i_n) else np.nan)
    ix_all = np.arange(len(allf))
    eo, er = exc(ix_all, yo, nyo), exc(ix_all, yr, nyr)
    bs = [exc(i, yr, nyr) - exc(i, yo, nyo) for i in _boot_idx(allf)]
    bs = [v for v in bs if np.isfinite(v)]
    out["excess"] = dict(original=float(eo), redo=float(er), change=float(er - eo),
                         ci=[float(v) for v in np.percentile(bs, [2.5, 97.5])])
    d = out["mover_rate"]
    if d["diff"] <= -0.05 and d["ci"][1] < 0 and out["excess"]["change"] < 0:
        out["decision"] = "re-refinement SUPPORTS the lab/refinement explanation"
    elif d["ci"][0] > -0.05:
        out["decision"] = "re-refinement does NOT support it"
    else:
        out["decision"] = "INCONCLUSIVE"
    with open(os.path.join("results", "labtest_redo.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("labels", "metadata", "analyse", "redo"))
    a = ap.parse_args()
    {"labels": labels, "metadata": metadata, "analyse": analyse, "redo": redo}[a.cmd]()


if __name__ == "__main__":
    main()
