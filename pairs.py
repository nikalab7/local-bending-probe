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
                  Common crystallization additives (ADDITIVES: SO4, GOL, BME,
                  EDO, ...; metals excluded) are ignored in that comparison:
                  they alone dropped ~150 clean mutations, and ignoring them
                  leaves the WT-vs-WT null calibration unchanged (3.1 vs 3.2%).
  4. AGGREGATION  one row per MUTATION: median bending over its crystals (the
                  form with the most crystals is used), not one row per PDB entry.
  5. NOISE FLOOR  per-window WT sigma = sample SD / c4(n) (unbiased). The first
                  version used 1.4826 * MAD with no finite-sample correction,
                  which underestimates sigma by a third at n = 3: 17-20% of
                  WT-vs-WT pseudo-mutants were "movers" in forms with 3-4 WT
                  crystals. Bias-corrected MAD and Qn are available
                  (SIGMA_EST) but calibrate worse on the null (6-7% and 6%
                  false positives vs 4% for SD at |z| > 2): WT crystals of one
                  form are a mixture (most agree, a few sit in another
                  conformation), a robust scale ignores the minority, and a
                  mutant crystal in that minority state then looks like a
                  mover. SD sees the minority, so its errors are conservative.
                  Sigma is then shrunk toward a prior. The default prior is heteroscedastic: within each
                  form, log(sigma) is regressed on the window's WT B-factor
                  (z within chain), because flexible windows are also noisier.
                  Shrinking every window toward ONE pooled sigma pulls high-B
                  windows' sigma down, inflates their z, and manufactures a
                  B-factor -> "mover" association. prior="pooled" keeps that
                  behaviour for comparison. The threshold uses the standard
                  error of a median-vs-median difference:
                      SE = sigma * sqrt(e_m^2 / m + e_n^2 / n)
                  with e_k = sd(median of k) * sqrt(k) / sigma (1 for k <= 2,
                  -> 1.2533 for large k). delta / SE is mapped from a Student-t
                  with n - 1 + SHRINK_K degrees of freedom onto the normal scale
                  (Z_CALIB = "t": sigma from few crystals is itself uncertain),
                  and mover := |z| > Z_MOVER. pairs.null_rows / diagnostics.py
                  check the false-positive rate empirically: 2.4% overall and
                  <= 4.9% in every n_wt bin (nominal 4.6%).
  6. WT LIGANDS   the ligand rule applies to the WT reference too: a WT crystal
                  whose het groups near a window differ from the form's typical
                  set is left out of that window's median and sigma, as long as
                  >= MIN_WT clean crystals remain (WT_LIG_FILTER = "prefer").

Window: the canonical 5-C-alpha window centred on the mutated residue
(start s = r - 2), measured with bending_metric.bending_angle.
"""
from __future__ import annotations
import gzip
import math
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
# z calibration for the uncertainty of sigma: "t" maps delta/SE through a
# Student-t with nu = n - 1 + SHRINK_K (the prior acts as SHRINK_K pseudo-
# crystals) onto the normal scale, so |z| > Z_MOVER means the same tail
# probability for a 3-crystal and a 30-crystal floor; "none" keeps delta/SE.
Z_CALIB = "t"
LIG_CUTOFF = 8.0       # Angstrom, het atom to any window C-alpha
LIG_TYPICAL = 0.5      # het present in >= this fraction of form-WT = "normal"
MED_EFF = 1.2533       # sd(median) / sd(mean) for normal data, large samples
# ligand rule for the WT reference, per window: "prefer" = use WT crystals in the
# form's typical ligand state when >= MIN_WT remain, else all; "strict" = typical
# state only (drops windows); "off" = all WT crystals
WT_LIG_FILTER = "prefer"
# finite-sample bias correction for 1.4826 * MAD (Croux & Rousseeuw 1992);
# n > 9: n / (n - 0.8). Checked by simulation in tests.
MAD_BN = {2: 1.196, 3: 1.495, 4: 1.363, 5: 1.206, 6: 1.200, 7: 1.140, 8: 1.129,
          9: 1.107}
# Qn finite-sample factors (Croux & Rousseeuw 1992); n > 9: n/(n+1.4) odd, n/(n+3.8) even
QN_DN = {2: 0.399, 3: 0.994, 4: 0.512, 5: 0.844, 6: 0.611, 7: 0.857, 8: 0.669,
         9: 0.872}
SIGMA_EST = "sd"       # per-window WT sigma: "sd" | "qn" | "mad" (see module doc)
# sd(median of k normals) * sqrt(k), by simulation (4e5 draws); k > 15 uses
# MED_EFF * sqrt(k / (k + 1.45)), which matches simulation to < 0.5%.
MED_EFF_K = {1: 1.0, 2: 1.0, 3: 1.162, 4: 1.094, 5: 1.196, 6: 1.137, 7: 1.215,
             8: 1.160, 9: 1.217, 10: 1.175, 11: 1.230, 12: 1.189, 13: 1.232,
             14: 1.195, 15: 1.238}
NULL_HELDOUT = 10      # WT crystals per form held out as pseudo-mutants (null control)
WATER = {"HOH", "DOD", "WAT", "H2O"}
# Crystallization additives / buffer components: ignored by the ligand-state
# rule when IGNORE_ADDITIVES. Metals are NOT listed (structural / catalytic
# Mg, Ca, Zn, Fe, Mn sites are real ligand states).
ADDITIVES = frozenset("""SO4 PO4 CL NA K BR IOD NO3 SCN NH4 GOL EDO PEG PGE PG4 1PE P6G
    PE4 PE5 2PE 12P 15P ACT ACY FMT BME HED DMS MPD MRD TRS EPE MES IMD IPA EOH MOH
    MLI TAR BU3 CAC CXS PGO DIO BTB B3P""".split())
IGNORE_ADDITIVES = True

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
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
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


def mad_sigma(vals):
    """Bias-corrected robust sigma: 1.4826 * MAD * b_n (0 for fewer than 2 values)."""
    vals = np.asarray(vals, float)
    n = len(vals)
    if n < 2:
        return 0.0
    s = 1.4826 * np.median(np.abs(vals - np.median(vals)))
    return float(s * MAD_BN.get(n, n / (n - 0.8)))


def qn_sigma(vals):
    """Rousseeuw-Croux Qn scale (82% efficient, robust), finite-sample corrected."""
    vals = np.asarray(vals, float)
    n = len(vals)
    if n < 2:
        return 0.0
    h = n // 2 + 1
    i, j = np.triu_indices(n, 1)
    d = np.sort(np.abs(vals[i] - vals[j]))
    dn = QN_DN.get(n, n / (n + 1.4) if n % 2 else n / (n + 3.8))
    return float(2.2219 * d[h * (h - 1) // 2 - 1] * dn)


def sd_sigma(vals):
    """Sample SD / c4(n): unbiased for normal data and the most efficient, but not
    robust -- an outlying WT crystal inflates sigma (fewer movers, never more)."""
    vals = np.asarray(vals, float)
    n = len(vals)
    if n < 2:
        return 0.0
    c4 = math.sqrt(2.0 / (n - 1)) * math.exp(math.lgamma(n / 2) - math.lgamma((n - 1) / 2))
    return float(vals.std(ddof=1) / c4)


def sigma_hat(vals):
    return {"mad": mad_sigma, "qn": qn_sigma, "sd": sd_sigma}[SIGMA_EST](vals)


def med_eff(k):
    """sd(median of k iid normals) / (sigma / sqrt(k))."""
    return MED_EFF_K.get(k) or float(MED_EFF * np.sqrt(k / (k + 1.45)))


def fit_sigma_prior(sig, bw, n, pooled):
    """Per-window prior sigma from WT B-factor: log(sigma) = a + b * bw.

    sig, bw, n : per-window raw sigma, WT window B z-score, WT crystal count.
    Slope by weighted least squares on log-sigma (weight n: more crystals,
    less noisy log-sigma; windows with sigma = 0 cannot enter the log fit).
    The intercept is set from the MEAN of sigma * exp(-b * bw) over all
    windows: sigma is unbiased after the MAD correction, but its distribution
    is right-skewed at small n, so a median-based intercept would put the
    prior below the true noise level. The prior is clipped to the 5-95th
    percentile of observed sigma so it never extrapolates. Too few windows or
    no B variation -> constant pooled prior, slope 0.
    Returns (prior array, slope).
    """
    sig, bw, n = (np.asarray(v, float) for v in (sig, bw, n))
    use = (sig > 0) & np.isfinite(bw)
    if use.sum() < PRIOR_MIN_WINDOWS or np.std(bw[use]) < 1e-6:
        return np.full(len(sig), pooled), 0.0
    ls = np.log(sig[use])
    b = float(np.polyfit(bw[use], ls, 1, w=np.sqrt(n[use]))[0])
    fin = np.isfinite(bw)
    a = float(np.log(np.mean(sig[fin] * np.exp(-b * bw[fin]))))
    lo, hi = np.percentile(sig[use], [5, 95])
    bwf = np.where(fin, bw, np.nanmedian(bw[use]))
    return np.clip(np.exp(a + b * bwf), lo, hi), b


def bend_of(st, s):
    """window_bend of the selected chain, cached on the structure dict."""
    cache = st.setdefault("_bend", {})
    key = (st["chid"], s)
    if key not in cache:
        cache[key] = window_bend(st["res"], s)
    return cache[key]


def bz_of(st):
    cache = st.setdefault("_bz", {})
    if st["chid"] not in cache:
        cache[st["chid"]] = chain_bz(st["res"])
    return cache[st["chid"]]


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


def het_map(st, res):
    """{resSeq: het names with an atom within LIG_CUTOFF of that C-alpha}.

    Computed once per (structure, chain) and cached on the structure dict.
    """
    cache = st.setdefault("_het_map", {})
    key = id(res)
    if key not in cache:
        out = {}
        if st["hets"] and res:
            keys = sorted(res)
            ca = np.array([res[k][1] for k in keys])
            names = [n for n, _ in st["hets"]]
            xyz = np.array([x for _, x in st["hets"]])
            close = np.linalg.norm(xyz[:, None, :] - ca[None, :, :], axis=2) < LIG_CUTOFF
            for i, j in zip(*np.nonzero(close)):
                out.setdefault(keys[j], set()).add(names[i])
        cache[key] = {k: frozenset(v) for k, v in out.items()}
    return cache[key]


def near_hets(st, res, s):
    """Het residue names with an atom within LIG_CUTOFF of window s..s+4
    (crystallization additives left out when IGNORE_ADDITIVES)."""
    m = het_map(st, res)
    names = set()
    for w in range(s, s + 5):
        if w in res:
            names |= m.get(w, frozenset())
    if IGNORE_ADDITIVES:
        names -= ADDITIVES
    return frozenset(names)


def form_stats(structs, wts, cons, prior):
    """WT window statistics of one crystal form, with the sigma prior.

    Returns dict(per={s: (median bend, raw sigma, n, WT window B z)}, pooled,
    prior={s: prior sigma}, slope, wts) or None (too few WT crystals).
    """
    if len(wts) < MIN_WT:
        return None
    bz = [bz_of(structs[p]) for p in wts]
    per = {}
    mode = {True: "strict", False: "off"}.get(WT_LIG_FILTER, WT_LIG_FILTER)
    for s in range(min(cons), max(cons) - 3):
        typical = typical_hets(structs, wts, s) if mode != "off" else None
        obs = []                                   # (bend, window B z, clean?)
        for p, z in zip(wts, bz):
            v = bend_of(structs[p], s)
            if v is None:
                continue
            zw = np.array([z[w] for w in range(s, s + 5)])
            clean = typical is None or near_hets(structs[p], structs[p]["res"], s) == typical
            obs.append((v, float(zw[np.isfinite(zw)].mean()) if np.isfinite(zw).any()
                        else np.nan, clean))
        clean_obs = [o for o in obs if o[2]]
        if mode == "strict" or (mode == "prefer" and len(clean_obs) >= MIN_WT):
            obs = clean_obs
        vals = [o[0] for o in obs]
        bws = [o[1] for o in obs]
        if len(vals) >= MIN_WT:
            vals = np.array(vals)
            med = float(np.median(vals))
            bws = np.array(bws)
            bw = float(np.median(bws[np.isfinite(bws)])) if np.isfinite(bws).any() else np.nan
            per[s] = (med, sigma_hat(vals), len(vals), bw)
    if not per or not any(v[1] > 0 for v in per.values()):
        return None
    pooled = float(np.mean([v[1] for v in per.values()]))
    keys = sorted(per)
    if prior == "bfactor":
        pr, slope = fit_sigma_prior([per[s][1] for s in keys], [per[s][3] for s in keys],
                                    [per[s][2] for s in keys], pooled)
    else:
        pr, slope = np.full(len(keys), pooled), 0.0
    return dict(per=per, pooled=pooled, wts=wts, slope=slope,
                prior=dict(zip(keys, map(float, pr))))


def typical_hets(structs, wts, s):
    """Het names near window s that at least LIG_TYPICAL of the WT crystals carry."""
    c = Counter()
    for p in wts:
        c.update(near_hets(structs[p], structs[p]["res"], s))
    return {k for k, v in c.items() if v >= LIG_TYPICAL * len(wts)}


def calibrated_z(t_stat, n):
    """delta/SE -> normal-scale z with the same two-sided tail under t(n - 1 + K)."""
    if Z_CALIB != "t":
        return float(t_stat)
    from scipy.stats import norm, t as student_t
    p = max(2.0 * student_t.sf(abs(t_stat), n - 1 + SHRINK_K), 1e-300)
    return float(np.sign(t_stat) * norm.isf(p / 2.0))


def label(W, s, bends):
    """(delta, se, z, sigma, sigma_prior) of mutant bends vs form stats W at window s."""
    wt_med, sig_hat, n, _ = W["per"][s]
    sig_prior = W["prior"][s]
    sig = float(np.sqrt((n * sig_hat ** 2 + SHRINK_K * sig_prior ** 2) / (n + SHRINK_K)))
    m = len(bends)
    se = float(sig * np.sqrt(med_eff(m) ** 2 / m + med_eff(n) ** 2 / n))
    delta = float(np.median(bends)) - wt_med
    return delta, se, calibrated_z(delta / se, n), sig, sig_prior


# ------------------------------ pair builder --------------------------------
def build_pairs(pdb_paths, pick_chain, protein, min_cons=5, prior=PRIOR,
                family=None, null_out=None):
    """Clean one-row-per-mutation labels for one protein.

    pdb_paths  : {pid: path}
    pick_chain : chains -> (chain_id, residues) or (None, None)
               or a dict {pid: chain_id} (chains fixed in advance, e.g. by SIFTS)
    prior      : noise-floor prior, "bfactor" or "pooled" (see module doc)
    family     : sequence-family label stored on every row (default: protein)
    null_out   : if a list, WT-vs-WT pseudo-mutant rows are appended to it
                 (see null_rows) -- the negative control for the site model
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
        if isinstance(pick_chain, dict):
            ch = pick_chain.get(pid)
            ch = ch if ch in st["chains"] else None
        else:
            ch, _ = pick_chain({c: {r: (v[0], v[1]) for r, v in rs.items()}
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
        W = form_stats(structs, wts, cons, prior)
        if W is not None:
            diag[f"prior_slope[{f}]"] = round(W["slope"], 3)
            wt_stats[f] = W

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
        typical = typical_hets(structs, W["wts"], s)
        bends, used, lig_drop = [], [], 0
        for p in pids:
            if near_hets(structs[p], structs[p]["res"], s) != typical:
                lig_drop += 1
                continue
            b = bend_of(structs[p], s)
            if b is not None:
                bends.append(b); used.append(p)
        diag["mutant_crystals_ligand_mismatch"] += lig_drop
        if not bends:
            diag["variant_form_no_clean_crystal"] += 1
            continue
        m = len(bends)
        delta, se, z, sig, sig_prior = label(W, s, bends)
        # best-resolution WT of this form that resolves the window, preferring
        # crystals in the form's typical ligand state there
        cands = [p for p in W["wts"] if bend_of(structs[p], s) is not None]
        clean_wt = [p for p in cands if near_hets(structs[p], structs[p]["res"], s) == typical]
        scaffold = min(clean_wt or cands,
                       key=lambda p: (structs[p]["resolution"], -len(structs[p]["res"])))
        cand[mut].append(dict(
            protein=protein, family=family or protein, r=r, wt=wtaa, mut=mutaa, form=f, s=s,
            n_wt=n, n_mut=m, wt_bend=wt_med, sigma=sig, sigma_raw=sig_hat,
            sigma_prior=sig_prior, b_window_wt=bw_wt,
            delta=delta, se=se, z=float(z),
            mover=bool(abs(z) > Z_MOVER),
            res_mut=float(np.median([structs[p]["resolution"] for p in used])),
            pdbs=sorted(used), scaffold=scaffold,
            scaffold_res=structs[scaffold]["res"],
            scaffold_path=pdb_paths[scaffold], scaffold_chain=structs[scaffold]["chid"],
            scaffold_helix=structs[scaffold]["helix"], scaffold_sheet=structs[scaffold]["sheet"],
            scaffold_ss=ss_of(structs[scaffold], structs[scaffold]["chid"], r),
            wt_seq="".join(cons[w] for w in range(s, s + 5))))
    rows = [max(v, key=lambda x: (x["n_mut"], x["n_wt"])) for v in cand.values()]
    rows.sort(key=lambda x: (x["r"], x["mut"]))
    diag["mutations"] = len(rows)
    diag["movers"] = sum(x["mover"] for x in rows)
    if null_out is not None:
        nr = null_rows(rows, structs, wt_by_form, cons, prior)
        diag["null_pseudo_pairs"] = len(nr)
        diag["null_pseudo_movers"] = sum(x["mover"] for x in nr)
        null_out.extend(nr)
    return rows, dict(diag)


def null_rows(rows, structs, wt_by_form, cons, prior):
    """WT-vs-WT pseudo-mutants at the windows of the real rows.

    For each form, up to NULL_HELDOUT WT crystals are held out one at a time;
    the held-out crystal plays a single-crystal "mutant" against statistics
    (median, sigma, prior) recomputed WITHOUT it, through the same label()
    code. No mutation exists, so every pseudo-mover is noise. If site
    features predict pseudo-movers, the site signal on real data can be
    noise structure; if they do not, it cannot.
    """
    out, done, cache = [], set(), {}
    for x in rows:
        f, s = x["form"], x["s"]
        wts = sorted(wt_by_form[f])
        held = [p for p in sorted(wts, key=lambda p: (structs[p]["resolution"], p))][:NULL_HELDOUT]
        for h in held:
            if (f, s, h) in done:
                continue
            done.add((f, s, h))
            rest = [p for p in wts if p != h]
            if (f, h) not in cache:
                cache[(f, h)] = form_stats(structs, rest, cons, prior)
            W = cache[(f, h)]
            if W is None or s not in W["per"]:
                continue
            if near_hets(structs[h], structs[h]["res"], s) != typical_hets(structs, rest, s):
                continue
            b = bend_of(structs[h], s)
            if b is None:
                continue
            delta, se, z, sig, sig_prior = label(W, s, [b])
            out.append({**x, "n_wt": W["per"][s][2], "n_mut": 1, "delta": delta,
                        "se": se, "z": float(z), "sigma": sig, "sigma_prior": sig_prior,
                        "mover": bool(abs(z) > Z_MOVER), "pdbs": [h], "heldout": h})
    return out


# --------------------------------- CLI --------------------------------------
def load_proteins(prior=PRIOR, null_out=None, mined=True):
    """Clean pairs for T4L + the Gate 4 proteins (gate caches) and, if
    mine_pairs.py has been run, every mined protein (manifests/mined.json).

    Returns {protein: (rows, diag)}. Rows carry a sequence-family label
    (>= 30% identity clusters from the miner; the protein name without it).
    """
    from feasibility_t4l import PDB_DIR, pick_t4l_chain
    from stats_utils import pinned_ids
    man = None
    if mined:
        from mine_pairs import load_manifest, MINED_DIR
        man = load_manifest()
    fam = {v["name"]: v["family"] for v in man["covered"].values()} if man else {}

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
        out["T4L"] = build_pairs(t4l, t4l_pick, "T4L", min_cons=10, prior=prior,
                                 family=fam.get("T4L"), null_out=null_out)
    try:
        from powered_loop_gate import PROTEINS, VAL_DIR, pick_chain, fetch_ids, MAX_PER
    except Exception:
        PROTEINS = []
    for up, name, hint in PROTEINS:
        ids = pinned_ids(f"val_{up}_{name}", lambda: fetch_ids(up))[:MAX_PER]
        paths = {p: os.path.join(VAL_DIR, f"{p}.pdb") for p in ids
                 if os.path.exists(os.path.join(VAL_DIR, f"{p}.pdb"))}
        if paths:
            out[name] = build_pairs(paths, lambda c, h=hint: pick_chain(c, h), name,
                                    prior=prior, family=fam.get(name), null_out=null_out)
    for acc, v in sorted((man or {}).get("proteins", {}).items()):
        paths = {p: os.path.join(MINED_DIR, f"{p}.pdb.gz") for p in v["entries"]}
        paths = {p: f for p, f in paths.items() if os.path.exists(f)}
        if paths:
            out[acc] = build_pairs(paths, {p: e["chain"] for p, e in v["entries"].items()},
                                   acc, min_cons=3, prior=prior, family=v["family"],
                                   null_out=null_out)
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
    fields = ["protein", "family", "r", "wt", "mut", "form", "n_wt", "n_mut", "wt_bend",
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
