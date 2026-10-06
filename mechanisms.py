"""
proteinX -- known-mechanism label check and backbone geometry (MECHANISMS.md).

Exploratory, after the lockbox run; the hypotheses, groups and thresholds were
committed in MECHANISMS.md before this was run. Data: dev + first-lockbox rows
(final_eval caches), frozen label unchanged.

  python mechanisms.py   -> results/mechanisms.json
"""
from __future__ import annotations
import json
import os
import numpy as np

import delta_model as dm
import final_eval as fe
import structure_features as sf
from posthoc_dev import rate_ci

PDB_DIRS = ("lockbox_pdb", "mined_pdb", "t4l_pdb", "val_pdb")
O_FLIP, O_FLIP_SENS, DIH_FLIP = 2.0, 1.5, 60.0
N_NULL, SEED = 1500, 0
MIN_GROUP, MIN_CLASS = fe.MIN_GROUP, fe.MIN_CLASS
_INDEX, _ATOMS, _CHAINS = {}, {}, {}


# ------------------------------ structures ----------------------------------
def path_of(pdb):
    if not _INDEX:
        for d in PDB_DIRS:
            if os.path.isdir(d):
                for f in os.listdir(d):
                    key = f.split(".")[0].upper()
                    _INDEX.setdefault(key, os.path.join(d, f))
    return _INDEX.get(pdb.upper())


def chains_of(path):
    if path not in _CHAINS:
        import gzip
        opener = gzip.open if path.endswith(".gz") else open
        ch = set()
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("ENDMDL"):
                    break
                if line.startswith("ATOM  "):
                    ch.add(line[21])
        _CHAINS[path] = sorted(ch)
    return _CHAINS[path]


def atoms(path, chain):
    k = (path, chain)
    if k not in _ATOMS:
        _ATOMS[k] = sf.parse_atoms(path, chain)[0]
    return _ATOMS[k]


def at(res, r, a):
    return res.get(r, (None, {}))[1].get(a)


def kabsch(P, Q):
    """Rotation R and translation t minimizing |(P R + t) - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    R = U @ D @ Vt
    return R, qc - pc @ R


def wrap(a):
    return (a + 180.0) % 360.0 - 180.0


def window_geometry(wt, mt, s):
    """C-alpha superposition of window s..s+4 and per-plane changes."""
    w = range(s, s + 5)
    P = [at(mt, i, "CA") for i in w]; Q = [at(wt, i, "CA") for i in w]
    if any(x is None for x in P + Q):
        return None
    P, Q = np.array(P), np.array(Q)
    R, t = kabsch(P, Q)
    rmsd = float(np.sqrt(((P @ R + t - Q) ** 2).sum(1).mean()))
    planes = []
    for i in range(s, s + 4):
        om, ow = at(mt, i, "O"), at(wt, i, "O")
        od = float(np.linalg.norm(om @ R + t - ow)) if om is not None and ow is not None else np.nan
        pm, pw = sf.phi_psi(mt, i)[1], sf.phi_psi(wt, i)[1]
        fm, fw = sf.phi_psi(mt, i + 1)[0], sf.phi_psi(wt, i + 1)[0]
        planes.append(dict(i=i, o_disp=od, dpsi=float(wrap(pm - pw)), dphi=float(wrap(fm - fw))))
    dih = [abs(v) for p in planes for v in (p["dpsi"], p["dphi"]) if np.isfinite(v)]
    return dict(ca_rmsd=rmsd, planes=planes, max_dihedral=max(dih) if dih else np.nan,
                max_o=float(np.nanmax([p["o_disp"] for p in planes])) if planes else np.nan)


def is_flip(p, o_thr=O_FLIP):
    return (np.isfinite(p["o_disp"]) and p["o_disp"] > o_thr
            and abs(p["dpsi"]) > DIH_FLIP and abs(p["dphi"]) > DIH_FLIP)


def row_geometry(x, expect_aa):
    """Median geometry over the row's mutant (or held-out WT) crystals vs the scaffold."""
    wt = atoms(x["scaffold_path"], x["scaffold_chain"])
    out = []
    for pdb in x["pdbs"]:
        path = path_of(pdb)
        if path is None:
            continue
        best = None
        for ch in chains_of(path):
            mt = atoms(path, ch)
            if mt.get(x["r"], ("?",))[0] != expect_aa:
                continue
            g = window_geometry(wt, mt, x["s"])
            if g and (best is None or g["ca_rmsd"] < best["ca_rmsd"]):
                best = g
        if best:
            out.append(best)
    if not out:
        return None
    flips = [any(is_flip(p) for p in g["planes"]) for g in out]
    flips_s = [any(is_flip(p, O_FLIP_SENS) for p in g["planes"]) for g in out]
    flip_planes = [(p["dpsi"], p["dphi"]) for g in out for p in g["planes"] if is_flip(p)]
    return dict(ca_rmsd=float(np.median([g["ca_rmsd"] for g in out])),
                max_o=float(np.median([g["max_o"] for g in out])),
                max_dihedral=float(np.median([g["max_dihedral"] for g in out])),
                flip=float(np.median(flips)) >= 0.5, flip_sens=float(np.median(flips_s)) >= 0.5,
                flip_planes=flip_planes)


# ------------------------------ mechanism groups ----------------------------
def scaffold_facts(x):
    wt = atoms(x["scaffold_path"], x["scaffold_chain"])
    r = x["r"]
    phi = sf.phi_psi(wt, r)[0]
    sg = at(wt, r, "SG")
    ss_bond = bool(sg is not None and any(
        v[0] == "C" and k != r and v[1].get("SG") is not None and np.linalg.norm(v[1]["SG"] - sg) < 2.5
        for k, v in wt.items()))
    a = [at(wt, r - 1, "CA"), at(wt, r - 1, "C"), at(wt, r, "N"), at(wt, r, "CA")]
    omega = sf.dihedral(*a) if all(v is not None for v in a) else np.nan
    helix = x["scaffold_helix"]; ch = x["scaffold_chain"]
    interior = all((ch, k) in helix for k in (r, r - 1, r - 2, r - 3))
    return dict(phi=phi, disulfide=ss_bond, cis=bool(np.isfinite(omega) and abs(omega) < 30),
                helix_interior=interior)


def groups(rows, facts, F):
    wt = np.array([x["wt"] for x in rows]); mu = np.array([x["mut"] for x in rows])
    phi = np.array([f["phi"] for f in facts], float)
    hse = dm.matrix(F, ["hse_up"])[:, 0]; blo = dm.matrix(F, ["blosum62"])[:, 0]
    gp = np.isin(wt, ["G", "P"]) | np.isin(mu, ["G", "P"])
    with np.errstate(invalid="ignore"):
        return {
            "M1 X->Pro inside a helix": (mu == "P") & np.array([f["helix_interior"] for f in facts]),
            "M2 Gly at positive phi -> X": (wt == "G") & (phi > 0) & (mu != "G"),
            "M3 X->Gly": mu == "G",
            "M4 Pro->X": wt == "P",
            "M5 disulfide Cys->X": (wt == "C") & np.array([f["disulfide"] for f in facts]),
            "M6 cis peptide at r": np.array([f["cis"] for f in facts]),
            "M7 conservative, exposed, no G/P (control)": (blo >= 1) & ~gp & (hse < np.nanmedian(hse)),
        }


def reportable(y):
    return len(y) >= MIN_GROUP and min(y.sum(), (~y).sum()) >= MIN_CLASS


def compare(v, fam, mA, mB, seed=0):
    """Rates of indicator v in groups A and B, family-bootstrap CI on each and on A - B."""
    v = np.asarray(v, float)
    idx = np.flatnonzero(mA | mB)
    fams = np.unique(fam[idx]); by = {f: idx[fam[idx] == f] for f in fams}
    rng = np.random.default_rng(seed); bs = []
    for _ in range(fe.N_BOOT):
        i = np.concatenate([by[f] for f in rng.choice(fams, len(fams))])
        a, b = v[i][mA[i]], v[i][mB[i]]
        if len(a) and len(b):
            bs.append((a.mean(), b.mean()))
    bs = np.array(bs); q = lambda c: [float(x) for x in np.percentile(c, [2.5, 97.5])]
    return dict(n_a=int(mA.sum()), n_b=int(mB.sum()), rate_a=float(v[mA].mean()), rate_b=float(v[mB].mean()),
                ci_a=q(bs[:, 0]), ci_b=q(bs[:, 1]), diff=float(v[mA].mean() - v[mB].mean()),
                diff_ci=q(bs[:, 0] - bs[:, 1]))


# ------------------------------ main ----------------------------------------
def main():
    rows_d, null_d, F_d = fe.load("dev")
    rows_l, null_l, F_l = fe.load("lockbox")
    rows, F = rows_d + rows_l, F_d + F_l
    source = np.array(["dev"] * len(rows_d) + ["lockbox"] * len(rows_l))
    fam = np.array([str(x["family"]) for x in rows])
    z = np.abs([x["z"] for x in rows]); y = z > 2
    yb = np.abs([x["z_bend"] for x in rows]) > 2
    facts = [scaffold_facts(x) for x in rows]
    G = groups(rows, facts, F)
    R = dict(note="exploratory; pre-specified in MECHANISMS.md", n=len(rows),
             n_dev=int((source == "dev").sum()), n_lockbox=int((source == "lockbox").sum()),
             part1={}, misses={}, part2={})
    print(f"rows {len(rows)} (dev {R['n_dev']}, lockbox {R['n_lockbox']}), movers {y.mean():.3f}")
    # ---- part 1
    for name, m in G.items():
        rec = dict(n=int(m.sum()), families=int(len(set(fam[m]))),
                   reportable=bool(reportable(y[m])))
        if m.sum() and (~m).sum():
            rec["combined"] = rate_ci(y, fam, m); rec["bend"] = rate_ci(yb, fam, m)
        R["part1"][name] = rec
        if "combined" in rec:
            c, b = rec["combined"], rec["bend"]
            print(f"  {name:44s} n={rec['n']:4d} fam={rec['families']:3d} combined {c['rate']:.2f} "
                  f"[{c['rate_ci'][0]:.2f},{c['rate_ci'][1]:.2f}] vs {c['rest_rate']:.2f} diff {c['difference']:+.2f} "
                  f"[{c['difference_ci'][0]:+.2f},{c['difference_ci'][1]:+.2f}] | bend {b['rate']:.2f} vs {b['rest_rate']:.2f}"
                  f"{'' if rec['reportable'] else '  (descriptive only)'}")
    # ---- geometry
    geo = [row_geometry(x, x["mut"]) for x in rows]
    rng = np.random.default_rng(SEED)
    null = null_d + null_l
    pick = rng.choice(len(null), min(N_NULL, len(null)), replace=False)
    ngeo = [row_geometry(null[i], null[i]["wt"]) for i in pick]
    nfam = np.array([str(null[i]["family"]) for i in pick])
    ok = np.array([g is not None for g in geo]); nok = np.array([g is not None for g in ngeo])
    print(f"geometry: rows {ok.sum()}/{len(rows)}, null {nok.sum()}/{len(ngeo)}")
    # misses of strong mechanisms
    for name in ("M1 X->Pro inside a helix", "M2 Gly at positive phi -> X",
                 "M5 disulfide Cys->X", "M6 cis peptide at r"):
        idx = np.flatnonzero(G[name] & ~y)
        R["misses"][name] = [dict(protein=rows[i]["protein"], mutation=f"{rows[i]['wt']}{rows[i]['r']}{rows[i]['mut']}",
                                  n_wt=rows[i]["n_wt"], n_mut=rows[i]["n_mut"], sigma=rows[i]["sigma"],
                                  z={k: round(rows[i][k], 2) for k in ("z_bend", "z_phi", "z_psi", "z_ca_tor_b")},
                                  max_dihedral_change=None if geo[i] is None else round(geo[i]["max_dihedral"], 1),
                                  ca_rmsd=None if geo[i] is None else round(geo[i]["ca_rmsd"], 2))
                             for i in idx]
    # ---- part 2: combine real rows and null rows into one indicator array
    flip = np.array([g["flip"] if g else np.nan for g in geo] + [g["flip"] if g else np.nan for g in ngeo], float)
    flip_s = np.array([g["flip_sens"] if g else np.nan for g in geo] + [g["flip_sens"] if g else np.nan for g in ngeo], float)
    allfam = np.concatenate([fam, np.array(["null:" + f for f in nfam])])
    n_r = len(rows)
    real = np.zeros(len(flip), bool); real[:n_r] = True
    okall = np.isfinite(flip)
    yy = np.zeros(len(flip), bool); yy[:n_r] = y
    yyb = np.zeros(len(flip), bool); yyb[:n_r] = yb
    dih_only = real & yy & ~yyb & okall
    bend_mov = real & yyb & okall
    nonmov = real & ~yy & okall
    nullm = ~real & okall
    R["part2"]["counts"] = dict(dihedral_only_movers=int(dih_only.sum()), bend_movers=int(bend_mov.sum()),
                                non_movers=int(nonmov.sum()), null=int(nullm.sum()))
    for lab, v in (("flip", flip), ("flip_sensitivity_1.5A", flip_s)):
        R["part2"][f"G1 {lab}"] = {f"dihedral-only vs {k}": compare(v, allfam, dih_only, m)
                                   for k, m in (("bend movers", bend_mov), ("non-movers", nonmov), ("null", nullm))}
    for k, c in R["part2"]["G1 flip"].items():
        print(f"  G1 flip rate {k:32s} {c['rate_a']:.3f} [{c['ci_a'][0]:.3f},{c['ci_a'][1]:.3f}] vs "
              f"{c['rate_b']:.3f} [{c['ci_b'][0]:.3f},{c['ci_b'][1]:.3f}] diff {c['diff']:+.3f} "
              f"[{c['diff_ci'][0]:+.3f},{c['diff_ci'][1]:+.3f}]  (n {c['n_a']} vs {c['n_b']})")
    # G2 crankshaft compensation in flip planes of real movers / null
    def comp(gs):
        pts = np.array([p for g in gs if g for p in g["flip_planes"]], float)
        if len(pts) < 3:
            return dict(n_planes=int(len(pts)))
        slope = float(np.polyfit(pts[:, 0], pts[:, 1], 1)[0])
        return dict(n_planes=int(len(pts)), corr=float(np.corrcoef(pts[:, 0], pts[:, 1])[0, 1]), slope=slope,
                    median_abs_sum=float(np.median(np.abs(wrap(pts[:, 0] + pts[:, 1])))))
    R["part2"]["G2 compensation"] = dict(movers=comp([geo[i] for i in np.flatnonzero(y & ok)]),
                                         null=comp(ngeo))
    print(f"  G2 {R['part2']['G2 compensation']}")
    # G3 Gly involvement among real rows
    gly = np.zeros(len(flip), bool)
    gly[:n_r] = np.array([x["wt"] == "G" or x["mut"] == "G" for x in rows])
    R["part2"]["G3 flip rate, Gly involved vs not"] = compare(flip, allfam, real & gly & okall, real & ~gly & okall)
    c = R["part2"]["G3 flip rate, Gly involved vs not"]
    print(f"  G3 flip rate Gly {c['rate_a']:.3f} vs no Gly {c['rate_b']:.3f} diff {c['diff']:+.3f} "
          f"[{c['diff_ci'][0]:+.3f},{c['diff_ci'][1]:+.3f}] (n {c['n_a']} vs {c['n_b']})")
    # G4 C-alpha RMSD
    rm = np.array([g["ca_rmsd"] if g else np.nan for g in geo] + [g["ca_rmsd"] if g else np.nan for g in ngeo])
    q = lambda m: [float(np.nanpercentile(rm[m], p)) for p in (25, 50, 75)]
    R["part2"]["G4 ca_rmsd quartiles"] = {"dihedral-only movers": q(dih_only), "bend movers": q(bend_mov),
                                          "non-movers": q(nonmov), "null": q(nullm)}
    print(f"  G4 C-alpha RMSD quartiles {R['part2']['G4 ca_rmsd quartiles']}")
    with open(os.path.join("results", "mechanisms.json"), "w") as fh:
        json.dump(R, fh, indent=1, default=float)
    print("wrote results/mechanisms.json")


if __name__ == "__main__":
    main()
