"""
proteinX -- CLEAN WT/mutant pair labels (implements SPEC sections 2.2-2.4).

The gate scripts built labels one crystal at a time, compared every mutant
against WT crystals of ANY form, never read resolution or ligands, and set the
mover threshold from a raw per-window MAD of as few as 3 crystals. On T4L that
let one variant (L99A, ~60 ligand-soak crystals) supply 25% of the validation
rows. This module fixes the labels:

  1. RESOLUTION   both sides <= MAX_RES (X-ray only; no CRYST1/resolution -> out).
  2. CRYSTAL FORM a mutant is compared ONLY to WT crystals of the same form
                  (space group + unit cell within CELL_TOL), so crystal packing
                  and hinge state are held constant.
  3. LIGANDS      a mutant crystal whose het groups near the window differ from
                  what WT crystals of that form typically carry is excluded.
  4. AGGREGATION  one row per MUTATION: median bending over its crystals (the
                  form with the most crystals is used), not one row per PDB entry.
  5. NOISE FLOOR  per-window WT sigma is shrunk toward a prior (few-crystal MADs
                  are unreliable). The default prior is heteroscedastic: within
                  each form, log(sigma) is regressed on the window's WT B-factor
                  (z within chain), because flexible windows are also noisier.
                  Shrinking every window toward ONE pooled sigma pulls high-B
                  windows' sigma down, inflates their z, and manufactures a
                  B-factor -> "mover" association. prior="pooled" keeps the old
                  behaviour for comparison. The threshold uses the standard
                  error of a median-vs-median difference:
                      SE = 1.2533 * sigma * sqrt(1/m + 1/n)
                  mover := |delta| / SE > Z_MOVER   (~5% false-positive rate).

Window: the canonical 5-C-alpha window centred on the mutated residue
(start s = r - 2), measured with bending_metric.bending_angle.
"""
from __future__ import annotations
import os
from collections import Counter, defaultdict
import numpy as np
from bending_metric import bending_angle, is_continuous

MAX_RES = 2.5          # Angstrom, SPEC 2.3
CELL_TOL = 0.02        # relative tolerance on a, b, c
ANGLE_TOL = 2.0        # degrees on alpha, beta, gamma
MIN_WT = 3             # WT crystals of the same form needed for a floor
SHRINK_K = 4           # pseudo-crystals of prior sigma added to each window
PRIOR = "bfactor"      # sigma prior: "bfactor" (log-linear in window B) | "pooled"
PRIOR_MIN_WINDOWS = 10 # fewer usable windows in a form -> fall back to pooled
Z_MOVER = 2.0
LIG_CUTOFF = 8.0       # Angstrom, het atom to any window C-alpha
LIG_TYPICAL = 0.5      # het present in >= this fraction of form-WT = "normal"
MED_EFF = 1.2533       # sd(median) / sd(mean) for normal data
WATER = {"HOH", "DOD", "WAT", "H2O"}

THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'}


# ------------------------------- parsing ------------------------------------
def parse_structure(path):
    """One pass over a PDB file.

    C-alpha selection mirrors feasibility_t4l.parse_ca exactly (model 1, ATOM
    records, standard residues, no insertion codes, highest occupancy) so the
    bending values are identical to the rest of the project.

    Returns dict(resolution, spacegroup, cell, helix, sheet, hets, chains) where
    chains = {chain: {resSeq: (aa, ca_xyz, ca_bfactor)}} and hets is a list of
    (resname, xyz) for every non-water HETATM atom.
    """
    out = dict(resolution=None, spacegroup=None, cell=None,
               helix=set(), sheet=set(), hets=[])
    best = {}
    with open(path) as fh:
        for line in fh:
            rec = line[:6]
            if rec.startswith("ENDMDL"):
                break
            if rec == "CRYST1":
                try:
                    out["cell"] = tuple(float(line[a:b]) for a, b in
                                        ((6, 15), (15, 24), (24, 33),
                                         (33, 40), (40, 47), (47, 54)))
                    out["spacegroup"] = line[55:66].strip() or None
                except ValueError:
                    pass
            elif rec == "REMARK" and line[6:10].strip() == "2" \
                    and "RESOLUTION." in line:
                tok = line.split("RESOLUTION.")[1].split()
                try:
                    out["resolution"] = float(tok[0])
                except (ValueError, IndexError):
                    pass
            elif rec == "HELIX ":
                try:
                    ch, a, b = line[19], int(line[21:25]), int(line[33:37])
                    out["helix"].update((ch, r) for r in range(a, b + 1))
                except ValueError:
                    pass
            elif rec == "SHEET ":
                try:
                    ch, a, b = line[21], int(line[22:26]), int(line[33:37])
                    out["sheet"].update((ch, r) for r in range(a, b + 1))
                except ValueError:
                    pass
            elif rec == "HETATM":
                name = line[17:20].strip()
                if name in WATER or name in THREE2ONE:
                    continue
                try:
                    out["hets"].append((name, np.array(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])])))
                except ValueError:
                    pass
            elif rec == "ATOM  ":
                if line[12:16].strip() != "CA":
                    continue
                aa = THREE2ONE.get(line[17:20].strip())
                if aa is None or line[26] != " ":
                    continue
                try:
                    rs = int(line[22:26])
                    xyz = np.array([float(line[30:38]), float(line[38:46]),
                                    float(line[46:54])])
                    occ = float(line[54:60]) if line[54:60].strip() else 1.0
                    bf = float(line[60:66]) if line[60:66].strip() else np.nan
                except ValueError:
                    continue
                key = (line[21], rs)
                if key not in best or occ > best[key][0]:
                    best[key] = (occ, aa, xyz, bf)
    chains = {}
    for (ch, rs), (_, aa, xyz, bf) in best.items():
        chains.setdefault(ch, {})[rs] = (aa, xyz, bf)
    out["chains"] = chains
    return out


def ss_of(st, ch, r):
    if (ch, r) in st["helix"]:
        return "H"
    if (ch, r) in st["sheet"]:
        return "E"
    return "L"


def chain_bz(res):
    """{resSeq: C-alpha B-factor z-scored within the chain} (NaN-safe)."""
    keys = sorted(res)
    bf = np.array([res[k][2] for k in keys], float)
    fin = np.isfinite(bf)
    if fin.sum() < 3 or bf[fin].std() == 0:
        return {k: np.nan for k in keys}
    z = (bf - bf[fin].mean()) / bf[fin].std()
    return dict(zip(keys, z))


def fit_sigma_prior(sig, bw, n, pooled):
    """Per-window prior sigma from WT B-factor: log(sigma) = a + b * bw.

    sig, bw, n : per-window raw sigma, WT window B z-score, WT crystal count.
    Slope by weighted least squares (weight n: more crystals, less noisy
    log-sigma); intercept set so the median residual is 0, which keeps the
    prior on the same scale as the pooled median. The prior is clipped to the
    5-95th percentile of observed sigma so it never extrapolates. Too few
    windows or no B variation -> constant pooled prior, slope 0.
    Returns (prior array, slope).
    """
    sig, bw, n = (np.asarray(v, float) for v in (sig, bw, n))
    use = (sig > 0) & np.isfinite(bw)
    if use.sum() < PRIOR_MIN_WINDOWS or np.std(bw[use]) < 1e-6:
        return np.full(len(sig), pooled), 0.0
    ls = np.log(sig[use])
    b = float(np.polyfit(bw[use], ls, 1, w=np.sqrt(n[use]))[0])
    a = float(np.median(ls - b * bw[use]))
    lo, hi = np.percentile(sig[use], [5, 95])
    bwf = np.where(np.isfinite(bw), bw, np.nanmedian(bw[use]))
    return np.clip(np.exp(a + b * bwf), lo, hi), b


def window_bend(res, s):
    """Bending of window s..s+4 or None (missing residue / chain break)."""
    win = [s + k for k in range(5)]
    if any(w not in res for w in win):
        return None
    pts = np.array([res[w][1] for w in win])
    if not is_continuous(pts):
        return None
    return bending_angle(pts)


# ----------------------------- crystal forms --------------------------------
def assign_forms(structs):
    """{pid: form_id}; form = space group + cell within tolerance (greedy)."""
    refs = defaultdict(list)          # spacegroup -> [reference cells]
    form = {}
    for pid in sorted(structs):
        st = structs[pid]
        sg, cell = st["spacegroup"], st["cell"]
        if sg is None or cell is None:
            continue
        for k, ref in enumerate(refs[sg]):
            if all(abs(cell[i] - ref[i]) <= CELL_TOL * ref[i] for i in range(3)) and \
               all(abs(cell[i] - ref[i]) <= ANGLE_TOL for i in range(3, 6)):
                form[pid] = f"{sg}#{k}"
                break
        else:
            refs[sg].append(cell)
            form[pid] = f"{sg}#{len(refs[sg]) - 1}"
    return form


def near_hets(st, res, s):
    """Het residue names with an atom within LIG_CUTOFF of window s..s+4."""
    ca = np.array([res[w][1] for w in range(s, s + 5) if w in res])
    if not len(ca) or not st["hets"]:
        return frozenset()
    names = set()
    for name, xyz in st["hets"]:
        if np.min(np.linalg.norm(ca - xyz, axis=1)) < LIG_CUTOFF:
            names.add(name)
    return frozenset(names)


# ------------------------------ pair builder --------------------------------
def build_pairs(pdb_paths, pick_chain, protein, min_cons=5, prior=PRIOR):
    """Clean one-row-per-mutation labels for one protein.

    pdb_paths  : {pid: path}
    pick_chain : chains -> (chain_id, residues) or (None, None)
    prior      : noise-floor prior, "bfactor" or "pooled" (see module doc)
    Returns (rows, diag). Each row holds the label (delta, se, z, mover) plus
    what delta_model needs to featurise the site (scaffold residues etc.).
    """
    diag = Counter()
    structs = {}
    for pid, path in pdb_paths.items():
        try:
            st = parse_structure(path)
        except Exception:
            diag["unparseable"] += 1
            continue
        ch, res = pick_chain({c: {r: (v[0], v[1]) for r, v in rs.items()}
                              for c, rs in st["chains"].items()})
        if ch is None:
            diag["no_matching_chain"] += 1
            continue
        st["chid"], st["res"] = ch, st["chains"][ch]
        structs[pid] = st
    diag["structures"] = len(structs)

    # consensus sequence -> WT vs single mutants (SPEC 2.2, N = 1)
    counts = defaultdict(Counter)
    for st in structs.values():
        for rs, v in st["res"].items():
            counts[rs][v[0]] += 1
    cons = {rs: c.most_common(1)[0][0] for rs, c in counts.items()
            if sum(c.values()) >= min_cons}
    kind = {}
    for pid, st in structs.items():
        m = tuple(sorted((rs, cons[rs], v[0]) for rs, v in st["res"].items()
                         if rs in cons and v[0] != cons[rs]))
        kind[pid] = "wt" if not m else (m[0] if len(m) == 1 else None)
    diag["wt_any_res"] = sum(k == "wt" for k in kind.values())
    diag["single_any_res"] = sum(isinstance(k, tuple) for k in kind.values())

    # resolution + crystal form
    ok = {p for p, st in structs.items()
          if st["resolution"] is not None and st["resolution"] <= MAX_RES}
    diag["dropped_resolution"] = len(structs) - len(ok)
    form = assign_forms({p: structs[p] for p in ok})
    ok = {p for p in ok if p in form}
    wt_by_form = defaultdict(list)
    mut_by_key = defaultdict(list)     # (form, mutation) -> [pid]
    for p in ok:
        if kind[p] == "wt":
            wt_by_form[form[p]].append(p)
        elif isinstance(kind[p], tuple):
            mut_by_key[(form[p], kind[p])].append(p)
    diag["forms_with_floor"] = sum(len(v) >= MIN_WT for v in wt_by_form.values())

    # per-form WT window statistics + sigma prior
    if prior not in ("bfactor", "pooled"):
        raise ValueError(f"unknown prior {prior!r}")
    wt_stats = {}
    for f, wts in wt_by_form.items():
        if len(wts) < MIN_WT:
            continue
        bz = [chain_bz(structs[p]["res"]) for p in wts]
        per = {}
        for s in range(min(cons), max(cons) - 3):
            vals, bws = [], []
            for p, z in zip(wts, bz):
                v = window_bend(structs[p]["res"], s)
                if v is not None:
                    vals.append(v)
                    zw = np.array([z[w] for w in range(s, s + 5)])
                    bws.append(float(zw[np.isfinite(zw)].mean()) if np.isfinite(zw).any()
                               else np.nan)
            if len(vals) >= MIN_WT:
                vals = np.array(vals)
                med = float(np.median(vals))
                bws = np.array(bws)
                bw = float(np.median(bws[np.isfinite(bws)])) if np.isfinite(bws).any() else np.nan
                per[s] = (med, float(1.4826 * np.median(np.abs(vals - med))), len(vals), bw)
        pos = [v[1] for v in per.values() if v[1] > 0]
        if not pos:
            continue
        pooled = float(np.median(pos))
        keys = sorted(per)
        if prior == "bfactor":
            pr, slope = fit_sigma_prior([per[s][1] for s in keys], [per[s][3] for s in keys],
                                        [per[s][2] for s in keys], pooled)
        else:
            pr, slope = np.full(len(keys), pooled), 0.0
        diag[f"prior_slope[{f}]"] = round(slope, 3)
        wt_stats[f] = dict(per=per, pooled=pooled, wts=wts,
                           prior=dict(zip(keys, map(float, pr))))

    # one candidate per (form, mutation); then keep the best-covered form
    cand = defaultdict(list)
    for (f, mut), pids in mut_by_key.items():
        r, wtaa, mutaa = mut
        s = r - 2
        if f not in wt_stats or s not in wt_stats[f]["per"]:
            diag["variant_form_no_wt_floor"] += 1
            continue
        if any(w not in cons for w in range(s, s + 5)):
            diag["variant_window_outside_consensus"] += 1
            continue
        W = wt_stats[f]
        wt_med, sig_hat, n, bw_wt = W["per"][s]
        # ligand state that is normal for this form at this window
        wt_lig = Counter()
        for p in W["wts"]:
            wt_lig.update(near_hets(structs[p], structs[p]["res"], s))
        typical = {k for k, c in wt_lig.items() if c >= LIG_TYPICAL * len(W["wts"])}
        bends, used, lig_drop = [], [], 0
        for p in pids:
            if near_hets(structs[p], structs[p]["res"], s) != typical:
                lig_drop += 1
                continue
            b = window_bend(structs[p]["res"], s)
            if b is not None:
                bends.append(b); used.append(p)
        diag["mutant_crystals_ligand_mismatch"] += lig_drop
        if not bends:
            diag["variant_form_no_clean_crystal"] += 1
            continue
        m = len(bends)
        sig_prior = W["prior"][s]
        sig = float(np.sqrt((n * sig_hat ** 2 + SHRINK_K * sig_prior ** 2) / (n + SHRINK_K)))
        se = MED_EFF * sig * np.sqrt(1.0 / m + 1.0 / n)
        delta = float(np.median(bends)) - wt_med
        # best-resolution WT of this form that actually resolves the window
        scaffold = min((p for p in W["wts"] if window_bend(structs[p]["res"], s) is not None),
                       key=lambda p: (structs[p]["resolution"], -len(structs[p]["res"])))
        cand[mut].append(dict(
            protein=protein, r=r, wt=wtaa, mut=mutaa, form=f, s=s,
            n_wt=n, n_mut=m, wt_bend=wt_med, sigma=sig, sigma_raw=sig_hat,
            sigma_prior=sig_prior, b_window_wt=bw_wt,
            delta=delta, se=float(se), z=float(delta / se),
            mover=bool(abs(delta) > Z_MOVER * se),
            res_mut=float(np.median([structs[p]["resolution"] for p in used])),
            pdbs=sorted(used), scaffold=scaffold,
            scaffold_res=structs[scaffold]["res"],
            scaffold_ss=ss_of(structs[scaffold], structs[scaffold]["chid"], r),
            wt_seq="".join(cons[w] for w in range(s, s + 5))))
    rows = [max(v, key=lambda x: (x["n_mut"], x["n_wt"])) for v in cand.values()]
    rows.sort(key=lambda x: (x["r"], x["mut"]))
    diag["mutations"] = len(rows)
    diag["movers"] = sum(x["mover"] for x in rows)
    return rows, dict(diag)


# --------------------------------- CLI --------------------------------------
def load_proteins(prior=PRIOR):
    """Clean pairs for T4L + the Gate 4 proteins, from the gate caches."""
    from feasibility_t4l import PDB_DIR, pick_t4l_chain
    from stats_utils import pinned_ids

    def t4l_pick(chains):
        best, sc = (None, None), 1e9
        for ch, res in chains.items():
            if pick_t4l_chain({ch: res}) is not None and abs(len(res) - 164) < sc:
                best, sc = (ch, res), abs(len(res) - 164)
        return best

    out = {}
    t4l = {f[:-4]: os.path.join(PDB_DIR, f) for f in os.listdir(PDB_DIR)
           if f.endswith(".pdb")} if os.path.isdir(PDB_DIR) else {}
    if t4l:
        out["T4L"] = build_pairs(t4l, t4l_pick, "T4L", min_cons=10, prior=prior)
    try:
        from powered_loop_gate import PROTEINS, VAL_DIR, pick_chain, fetch_ids, MAX_PER
    except Exception:
        return out
    for up, name, hint in PROTEINS:
        ids = pinned_ids(f"val_{up}_{name}", lambda: fetch_ids(up))[:MAX_PER]
        paths = {p: os.path.join(VAL_DIR, f"{p}.pdb") for p in ids
                 if os.path.exists(os.path.join(VAL_DIR, f"{p}.pdb"))}
        if paths:
            out[name] = build_pairs(paths, lambda c, h=hint: pick_chain(c, h), name,
                                    prior=prior)
    return out


def main():
    import argparse
    import csv
    ap = argparse.ArgumentParser()
    ap.add_argument("--prior", choices=("bfactor", "pooled"), default=PRIOR,
                    help="noise-floor sigma prior (default: %(default)s)")
    args = ap.parse_args()
    data = load_proteins(args.prior)
    if not data:
        print("no cached PDB files -- run feasibility_t4l.py / powered_loop_gate.py first")
        return
    fields = ["protein", "r", "wt", "mut", "form", "n_wt", "n_mut", "wt_bend",
              "sigma_raw", "sigma_prior", "b_window_wt", "sigma", "delta", "se", "z", "mover", "res_mut",
              "scaffold_ss", "scaffold", "pdbs"]
    with open("pairs_clean.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for name, (rows, diag) in data.items():
            print(f"\n== {name} ==")
            for k, v in diag.items():
                print(f"  {k:36s} {v}")
            for x in rows:
                w.writerow({k: (" ".join(x[k]) if k == "pdbs" else x[k]) for k in fields})
    print("\nwrote pairs_clean.csv")


if __name__ == "__main__":
    main()
