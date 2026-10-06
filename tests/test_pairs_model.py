"""
Offline end-to-end tests for pairs.py + delta_model.py on SYNTHETIC PDB files.

Run:  python -m pytest tests/   (or: python tests/test_pairs_model.py)

The synthetic set is built so each label-cleaning rule has something to catch:
  * true movers         : 6 deg hinge at C-alpha r in the mutant
  * crystal-form effect : form B bends residue 40 in WT AND mutant; a mutant
                          seen only in form B must NOT be called a mover
  * ligand soaks        : 10 of 12 crystals of one variant carry a ligand next
                          to the window -> only the 2 apo crystals are used
  * low resolution      : a 3.0 A mutant crystal is dropped
"""
import os
import sys
import tempfile
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import pairs                                   # noqa: E402
import delta_model                             # noqa: E402

# The synthetic structures carry C-alphas only, so the label-cleaning tests run on
# the canonical single "bend" metric; the frozen combined label is unit-tested below.
pairs.METRIC = "bend"
pairs.NCS_AVERAGE = False

AA3 = {v: k for k, v in pairs.THREE2ONE.items()}
N_RES = 60
SEQ = "".join(np.random.default_rng(7).choice(list("AVLIFMSTNQDEKRG"), N_RES))
FORMS = {"A": ("P 21 21 21", (50.0, 60.0, 70.0, 90.0, 90.0, 90.0)),
         "B": ("P 1 21 1", (40.0, 50.0, 60.0, 90.0, 100.0, 90.0))}
MOVER_SITES = [12, 20, 28]
LIGAND_SITE = 33
FORM_B_SITE = 40            # form B carries an 8 deg hinge here in every crystal


def ideal_helix():
    t = np.radians(100.0) * np.arange(1, N_RES + 1)
    return np.stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(1, N_RES + 1)], 1)


def hinge(ca, r, deg):
    """Rotate residues after r about the local tangent through C-alpha r.

    The axis is perpendicular to the helix axis and the radial vector, so the
    effect on the bend is the same at every site. Bond lengths are unchanged.
    """
    o = ca[r - 1]
    ax = np.cross([0.0, 0.0, 1.0], np.r_[o[:2], 0.0]); ax /= np.linalg.norm(ax)
    t = np.radians(deg)
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    R = np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K        # Rodrigues
    out = ca.copy()
    out[r:] = (ca[r:] - o) @ R.T + o
    return out


def write_pdb(path, seq, ca, form, resolution, ligand_at=None, bfac=None):
    sg, cell = FORMS[form]
    lines = [f"REMARK   2 RESOLUTION.    {resolution:4.2f} ANGSTROMS.",
             "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} {:<11s}{:4d}".format(*cell, sg, 4)]
    h = list(" " * 40)
    h[0:5] = "HELIX"; h[19] = "A"; h[21:25] = f"{1:4d}"; h[33:37] = f"{N_RES:4d}"
    lines.append("".join(h))
    for i, (aa, xyz) in enumerate(zip(seq, ca), start=1):
        b = 20.0 if bfac is None else bfac[i - 1]
        lines.append(f"ATOM  {i:5d}  CA  {AA3[aa]} A{i:4d}    "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{b:6.2f}           C")
    lines.append(f"HETATM{9000:5d}  O   HOH A 900    {99.0:8.3f}{99.0:8.3f}{99.0:8.3f}  1.00 30.00")
    if ligand_at is not None:
        x = ca[ligand_at - 1] + np.array([3.0, 0.0, 0.0])
        lines.append(f"HETATM{9001:5d}  C1  BNZ A 901    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00 30.00")
    lines.append("END")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def make_dataset(d):
    rng = np.random.default_rng(0)
    base = ideal_helix()
    paths = {}
    k = 0

    def crystal(form, seq, move=None, res=1.8, lig=None):
        nonlocal k
        ca = base + rng.normal(0, 0.05, base.shape)
        if form == "B":                                   # packing-induced bend
            ca = hinge(ca, FORM_B_SITE, 8.0)
        if move is not None:                               # 6 deg hinge at C-alpha r
            ca = hinge(ca, move, 6.0)
        pid = f"S{k:03d}"; k += 1
        write_pdb(os.path.join(d, pid + ".pdb"), seq, ca, form, res, lig)
        paths[pid] = os.path.join(d, pid + ".pdb")

    def mutate(r, aa):
        return SEQ[:r - 1] + aa + SEQ[r:]

    for _ in range(8):
        crystal("A", SEQ)
    for _ in range(4):
        crystal("B", SEQ)
    for r in MOVER_SITES:                                  # true movers
        for _ in range(2):
            crystal("A", mutate(r, "P" if SEQ[r - 1] != "P" else "G"), move=r)
    for r in (16, 24):                                     # null mutants
        crystal("A", mutate(r, "W" if SEQ[r - 1] != "W" else "Y"))
    m33 = mutate(LIGAND_SITE, "A" if SEQ[LIGAND_SITE - 1] != "A" else "G")
    for _ in range(10):                                    # ligand soaks
        crystal("A", m33, move=LIGAND_SITE, lig=LIGAND_SITE)
    for _ in range(2):                                     # apo
        crystal("A", m33)
    crystal("B", mutate(FORM_B_SITE, "W" if SEQ[FORM_B_SITE - 1] != "W" else "Y"))
    crystal("A", mutate(45, "W" if SEQ[44] != "W" else "Y"), move=45, res=3.0)  # low-res
    return paths


def build():
    d = tempfile.mkdtemp()
    paths = make_dataset(d)
    return pairs.build_pairs(paths, lambda c: ("A", c["A"]), "synthetic", min_cons=5)


def test_parser_roundtrip():
    d = tempfile.mkdtemp(); p = os.path.join(d, "x.pdb")
    write_pdb(p, SEQ, ideal_helix(), "B", 2.1, ligand_at=10)
    st = pairs.parse_structure(p)
    assert st["resolution"] == 2.1
    assert st["spacegroup"] == "P 1 21 1"
    assert np.allclose(st["cell"], FORMS["B"][1])
    assert [n for n, _ in st["hets"]] == ["BNZ"]           # water ignored
    assert ("A", 30) in st["helix"]
    assert len(st["chains"]["A"]) == N_RES
    from feasibility_t4l import parse_ca                    # same C-alpha selection
    ca_old = parse_ca(p)["A"]
    assert all(np.allclose(ca_old[r][1], st["chains"]["A"][r][1]) for r in ca_old)


def test_label_cleaning():
    rows, diag = build()
    by = {x["r"]: x for x in rows}
    for r in MOVER_SITES:
        assert by[r]["mover"], (r, by[r]["z"])
        assert by[r]["n_mut"] == 2                           # aggregated, one row
    for r in (16, 24):
        assert not by[r]["mover"], (r, by[r]["z"])
    # crystal form: compared to form-B WT, so the form shift is not a "mover"
    assert by[FORM_B_SITE]["form"].startswith("P 1 21 1")
    assert not by[FORM_B_SITE]["mover"], by[FORM_B_SITE]["z"]
    # ligand soaks excluded; the variant is ONE row built from the 2 apo crystals
    assert by[LIGAND_SITE]["n_mut"] == 2
    assert not by[LIGAND_SITE]["mover"]
    assert diag["mutant_crystals_ligand_mismatch"] == 10
    # low-resolution crystal dropped
    assert 45 not in by and diag["dropped_resolution"] == 1
    assert len(rows) == len({(x["r"], x["mut"]) for x in rows})


def test_scaffold_resolves_window():
    """The best-resolution WT lacks residues 10-14; it must not be the scaffold there."""
    d = tempfile.mkdtemp(); paths = make_dataset(d)
    ca = ideal_helix()
    keep = [i for i in range(N_RES) if not 10 <= i + 1 <= 14]
    p = os.path.join(d, "GAP.pdb")
    write_pdb(p, "".join(SEQ[i] for i in keep), ca[keep], "A", 0.9)
    # renumber: write_pdb numbers sequentially, so patch residue numbers back
    lines = open(p).read().splitlines()
    it = iter(i + 1 for i in keep)
    lines = [l[:22] + f"{next(it):4d}" + l[26:] if l.startswith("ATOM") else l for l in lines]
    open(p, "w").write("\n".join(lines) + "\n")
    paths["GAP"] = p
    rows, _ = pairs.build_pairs(paths, lambda c: ("A", c["A"]), "synthetic", min_cons=5)
    x = [x for x in rows if x["r"] == 12][0]
    assert x["scaffold"] != "GAP"
    delta_model.featurize([x])


def test_sigma_prior_fit():
    rng = np.random.default_rng(1)
    bw = rng.normal(size=200)
    sig = np.exp(0.0 + 0.5 * bw + 0.2 * rng.normal(size=200))
    pr, slope = pairs.fit_sigma_prior(sig, bw, np.full(200, 8), pooled=1.0)
    assert abs(slope - 0.5) < 0.1
    assert pr[np.argmax(bw)] > pr[np.argmin(bw)]
    # no B variation (or too few windows) -> constant pooled prior
    pr, slope = pairs.fit_sigma_prior(sig, np.zeros(200), np.full(200, 8), pooled=1.3)
    assert slope == 0.0 and np.all(pr == 1.3)
    pr, slope = pairs.fit_sigma_prior(sig[:5], bw[:5], np.full(5, 8), pooled=1.3)
    assert slope == 0.0 and np.all(pr == 1.3)


def test_bfactor_prior_tracks_flexible_region():
    """WT crystals noisier where B is high -> the prior is larger there, and a
    null mutant in the flexible region is less likely to look like a mover."""
    d = tempfile.mkdtemp(); rng = np.random.default_rng(5)
    base = ideal_helix()
    bfac = np.where((np.arange(1, N_RES + 1) >= 40) & (np.arange(1, N_RES + 1) <= 55), 60.0, 15.0)
    noise = np.where(bfac > 30, 0.25, 0.03)[:, None]
    paths = {}
    for k in range(6):
        ca = base + rng.normal(0, 1, base.shape) * noise
        p = os.path.join(d, f"W{k}.pdb"); write_pdb(p, SEQ, ca, "A", 1.8, bfac=bfac)
        paths[f"W{k}"] = p
    for r in (20, 48):
        aa = "W" if SEQ[r - 1] != "W" else "Y"
        ca = base + rng.normal(0, 1, base.shape) * noise
        p = os.path.join(d, f"M{r}.pdb")
        write_pdb(p, SEQ[:r - 1] + aa + SEQ[r:], ca, "A", 1.8, bfac=bfac)
        paths[f"M{r}"] = p
    out = {}
    for prior in ("bfactor", "pooled"):
        rows, diag = pairs.build_pairs(paths, lambda c: ("A", c["A"]), "synthetic",
                                       min_cons=5, prior=prior)
        out[prior] = {x["r"]: x for x in rows}
    bf, po = out["bfactor"], out["pooled"]
    assert bf[48]["sigma_prior"] > 2 * bf[20]["sigma_prior"]
    assert po[48]["sigma_prior"] == po[20]["sigma_prior"]
    assert abs(bf[48]["z"]) < abs(po[48]["z"])            # flexible: less inflated z
    assert bf[48]["b_window_wt"] > bf[20]["b_window_wt"]


def test_null_control_rows():
    """WT-vs-WT pseudo-mutants: held-out crystal excluded from its own floor,
    and (no mutation exists) few pseudo-movers."""
    d = tempfile.mkdtemp(); paths = make_dataset(d)
    null = []
    rows, diag = pairs.build_pairs(paths, lambda c: ("A", c["A"]), "synthetic",
                                   min_cons=5, null_out=null)
    assert null and diag["null_pseudo_pairs"] == len(null)
    real = {(x["form"], x["s"]): x for x in rows}
    for x in null:
        assert x["heldout"] not in real[(x["form"], x["s"])]["pdbs"]
        assert x["n_wt"] == real[(x["form"], x["s"])]["n_wt"] - 1
        assert x["n_mut"] == 1
    assert np.mean([x["mover"] for x in null]) < 0.2
    assert len({(x["form"], x["s"], x["heldout"]) for x in null}) == len(null)


def test_chain_dict_and_gzip():
    import gzip
    import shutil
    d = tempfile.mkdtemp(); paths = make_dataset(d)
    gz = {}
    for pid, p in paths.items():
        with open(p, "rb") as a, gzip.open(p + ".gz", "wb") as b:
            shutil.copyfileobj(a, b)
        gz[pid] = p + ".gz"
    rows_a, _ = pairs.build_pairs(paths, lambda c: ("A", c["A"]), "synthetic", min_cons=5)
    rows_b, diag = pairs.build_pairs(gz, {p: "A" for p in gz}, "synthetic", min_cons=5,
                                     family="FAM")
    assert [(x["r"], x["mut"], round(x["z"], 9)) for x in rows_a] == \
        [(x["r"], x["mut"], round(x["z"], 9)) for x in rows_b]
    assert all(x["family"] == "FAM" for x in rows_b)
    _, diag = pairs.build_pairs(gz, {p: "Z" for p in gz}, "synthetic", min_cons=5)
    assert diag["no_matching_chain"] == len(gz)


def test_prescreen_picks_usable_forms():
    import mine_pairs
    ref = "".join(np.random.default_rng(2).choice(list("ACDEFHIKLMNQRST"), 120))   # no G/P/W/Y
    mut = lambda i, a: ref[:i] + a + ref[i + 1:]
    k = iter(range(1000))

    def rec(seq, form="A", nent=1, res=1.8):
        sg, cell = FORMS[form]
        return dict(entity=f"X{next(k):03d}_1", pdb=f"X{next(k):03d}", chain="A", seq=seq,
                    resolution=res, n_protein_entities=nent, spacegroup=sg, cell=cell,
                    description="demo")
    recs = [rec(ref) for _ in range(4)]                          # form A: 4 WT
    recs += [rec(mut(10, "W")), rec(mut(20, "P")), rec(mut(30, "G"), nent=2)]
    recs += [rec(ref, "B") for _ in range(2)] + [rec(mut(40, "W"), "B")]  # B: 2 WT only
    recs += [rec(mut(50, "W")[:-1] + "Y")]                        # double mutant
    ps = mine_pairs.prescreen(recs)
    assert ps["reference"] == ref and list(ps["forms"]) == ["P 21 21 21#0"]
    f = ps["forms"]["P 21 21 21#0"]
    assert len(f["wt"]) == 4 and len(f["single"]) == 2          # complex + double out
    assert mine_pairs.hamming1("ABC", "ABD") == 2 and mine_pairs.hamming1("ABC", "AXD") is None
    assert mine_pairs.prescreen([rec(ref) for _ in range(5)]) is None   # no mutant


def test_form_matching_matters():
    """Sanity: against form-A WT the form-B mutant WOULD look like a mover."""
    rows, _ = build()
    x = [x for x in rows if x["r"] == FORM_B_SITE][0]
    d = tempfile.mkdtemp(); paths = make_dataset(d)
    st = {p: pairs.parse_structure(v) for p, v in paths.items()}
    wt_a = [pairs.window_bend(s["chains"]["A"], x["s"]) for s in st.values()
            if s["spacegroup"] == "P 21 21 21" and "".join(v[0] for _, v in
            sorted(s["chains"]["A"].items())) == SEQ]
    cross = np.median([pairs.window_bend(st[p]["chains"]["A"], x["s"]) for p in x["pdbs"]]) \
        - np.median(wt_a)
    assert abs(cross) > 2 * abs(x["delta"]) + 1.0


def test_fast_auc_matches_sklearn():
    from sklearn.metrics import roc_auc_score
    from stats_utils import fast_auc
    rng = np.random.default_rng(3)
    for _ in range(20):
        y = rng.random(80) > 0.6
        s = np.round(rng.normal(size=80), 1)                 # force ties
        assert abs(fast_auc(y, s) - roc_auc_score(y, s)) < 1e-12


def test_blosum_symmetric():
    assert all(delta_model.BLOSUM62[(a, b)] == delta_model.BLOSUM62[(b, a)]
               for a in delta_model.AA for b in delta_model.AA)
    assert delta_model.BLOSUM62[("W", "W")] == 11 and delta_model.BLOSUM62[("A", "R")] == -1


def test_featurize_runs():
    rows, _ = build()
    F = delta_model.featurize(rows)
    X = delta_model.matrix(F, delta_model.SITE + delta_model.SUBST)
    assert X.shape == (len(rows), len(delta_model.SITE) + len(delta_model.SUBST))
    assert np.isfinite(X[:, delta_model.SITE.index("n_ca10")]).all()


def _synthetic_features(signal, n_sites=60, per_site=3, seed=0):
    rng = np.random.default_rng(seed)
    F, y, g, prot = [], [], [], []
    for s in range(n_sites):
        site_flex = rng.normal()
        for _ in range(per_site):
            f = {c: rng.normal() for c in delta_model.SITE + delta_model.SUBST
                 + delta_model.CONTEXT + delta_model.INTERACT}
            f["b_site"] = site_flex + 0.3 * rng.normal()
            F.append(f)
            y.append(signal * site_flex + rng.normal() > 0.8)
            g.append(f"P{s % 2}:{s}"); prot.append(f"P{s % 2}")
    return F, np.array(y), np.array(g), prot


def test_cv_finds_planted_site_signal():
    F, y, g, prot = _synthetic_features(signal=2.0)
    res, con = delta_model.evaluate(F, y, g, prot, repeats=2)
    assert res["site|logreg"]["leave_site_out"]["auc"] > 0.75
    assert res["site|logreg"]["leave_site_out"]["lo"] > 0.6


def test_cv_null_is_chance():
    F, y, g, prot = _synthetic_features(signal=0.0)
    res, con = delta_model.evaluate(F, y, g, prot, repeats=2)
    m = res["site+subst|logreg"]["leave_site_out"]
    assert m["lo"] < 0.5 < m["hi"] or abs(m["auc"] - 0.5) < 0.12


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("PASS", name)


# ----------------------------- calibration ----------------------------------
def test_sigma_estimators_unbiased():
    """Finite-sample corrected MAD, Qn and SD/c4 are unbiased for normal data
    (the uncorrected 1.4826*MAD is ~0.67 sigma at n = 3)."""
    rng = np.random.default_rng(11)
    for n in (3, 5, 8, 12):
        v = rng.normal(size=(4000, n))
        for est in (pairs.mad_sigma, pairs.qn_sigma, pairs.sd_sigma):
            m = np.mean([est(row) for row in v])
            assert abs(m - 1.0) < 0.05, (est.__name__, n, m)
    raw = np.mean([1.4826 * np.median(np.abs(r - np.median(r))) for r in rng.normal(size=(4000, 3))])
    assert raw < 0.75


def test_median_efficiency_table():
    rng = np.random.default_rng(12)
    for k in (1, 2, 3, 4, 7, 20, 30):
        sim = np.median(rng.normal(size=(100000, k)), axis=1).std() * np.sqrt(k)
        assert abs(pairs.med_eff(k) - sim) / sim < 0.02, (k, pairs.med_eff(k), sim)


def test_label_se_uses_exact_efficiency():
    W = dict(per={10: (100.0, 1.0, 9, 0.0)}, prior={10: 1.0})
    delta, se, z, sig, _ = pairs.label(W, 10, [102.0])
    assert sig == 1.0 and delta == 2.0
    assert abs(se - np.sqrt(1.0 + pairs.med_eff(9) ** 2 / 9)) < 1e-12


# ------------------------- full-atom features --------------------------------
def _place(a, b, c, bond, angle, torsion):
    """NeRF: atom d with |cd| = bond, angle(b, c, d) = angle, dihedral(a, b, c, d) = torsion."""
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc); n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    ang, tor = np.radians(angle), np.radians(torsion)
    d2 = np.array([-bond * np.cos(ang), bond * np.sin(ang) * np.cos(tor),
                   bond * np.sin(ang) * np.sin(tor)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def _backbone(phis, psis):
    """Ideal backbone N, CA, C (+ O) for the given phi/psi lists, omega = 180."""
    N = np.array([0.0, 0.0, 0.0]); CA = np.array([1.458, 0.0, 0.0])
    C = _place(np.array([0.0, 1.0, 0.0]), N, CA, 1.525, 111.0, -60.0)
    atoms = [{"N": N, "CA": CA, "C": C}]
    for i in range(1, len(phis)):
        p = atoms[-1]
        n = _place(p["N"], p["CA"], p["C"], 1.329, 116.2, psis[i - 1])
        ca = _place(p["CA"], p["C"], n, 1.458, 121.7, 180.0)
        c = _place(p["C"], n, ca, 1.525, 111.0, phis[i])
        atoms.append({"N": n, "CA": ca, "C": c})
    for i in range(len(atoms) - 1):
        atoms[i]["O"] = _place(atoms[i + 1]["N"], atoms[i]["CA"], atoms[i]["C"], 1.231, 120.5, 180.0)
    return atoms


def _write_backbone_pdb(path, atoms, seq):
    import structure_features as sf
    three = {v: k for k, v in sf.THREE2ONE.items()}
    lines, k = [], 1
    for i, (at, aa) in enumerate(zip(atoms, seq), start=1):
        for name, xyz in at.items():
            lines.append(f"ATOM  {k:5d}  {name:<3s} {three[aa]} A{i:4d}    "
                         f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00           {name[0]}")
            k += 1
    open(path, "w").write("\n".join(lines + ["END"]) + "\n")


def test_dihedral_sign_convention():
    import structure_features as sf
    for th in (-120.0, -60.0, 30.0, 150.0):
        t = np.radians(th)
        d = sf.dihedral(np.array([1.0, 0, 0]), np.zeros(3), np.array([0, 0, 1.0]),
                        np.array([np.cos(t), np.sin(t), 1.0]))
        assert abs(d - th) < 1e-9


def test_phi_psi_from_built_backbone():
    import structure_features as sf
    phis = [-57.0, -57.0, -120.0, 60.0, -57.0, -65.0, -57.0]
    psis = [-47.0, -47.0, 130.0, 45.0, -47.0, 140.0, -47.0]
    atoms = _backbone(phis, psis)
    d = tempfile.mkdtemp(); p = os.path.join(d, "bb.pdb")
    _write_backbone_pdb(p, atoms, "AAGGAPA")
    residues, het = sf.parse_atoms(p, "A")
    assert len(residues) == 7 and len(het) == 0
    for r in range(2, 7):
        phi, psi = sf.phi_psi(residues, r)
        assert abs(phi - phis[r - 1]) < 0.5 and abs(psi - psis[r - 1]) < 0.5, (r, phi, psi)
    ctx = sf.context_features(residues, het, 4, "L", None)
    assert ctx["phi_pos"] == 1.0 and ctx["phi_alpha"] == 0.0
    ctx3 = sf.context_features(residues, het, 3, "L", None)
    assert ctx3["phi_beta"] == 1.0


def test_interaction_features_logic():
    import structure_features as sf
    base = dict(phi=70.0, phi_pos=1.0, cb_density=10.0, sc_bb_hbonds=2.0, helix_nterm=0.0)
    f = sf.interaction_features(base, "G", "A", "L")
    assert f["gly_phi_pos_loss"] == 1.0 and f["helix_prop_change"] == 0.0
    assert sf.interaction_features(base, "G", "G", "L")["gly_phi_pos_loss"] == 0.0
    assert sf.interaction_features({**base, "phi_pos": 0.0, "phi": -65.0}, "G", "A", "L")["gly_phi_pos_loss"] == 0.0
    hel = {**base, "phi": -60.0, "phi_pos": 0.0}
    f = sf.interaction_features(hel, "A", "P", "H")
    assert f["pro_in_helix"] == 1.0 and f["pro_phi_strain"] < 0.1
    assert abs(f["helix_prop_change"] - 3.16) < 1e-9
    assert sf.interaction_features({**hel, "phi": -150.0}, "A", "P", "L")["pro_phi_strain"] > 1.0
    assert sf.interaction_features(base, "S", "A", "L")["sc_hb_loss"] == 2.0      # loses both
    assert sf.interaction_features(base, "S", "T", "L")["sc_hb_loss"] == 0.0      # can still H-bond
    assert sf.interaction_features(base, "A", "W", "L")["overpack"] > 0.0
    assert sf.interaction_features(base, "W", "A", "L")["overpack"] == 0.0


def test_within_protein_auc_ignores_base_rates():
    from stats_utils import within_auc, fast_auc
    rng = np.random.default_rng(4)
    strata = np.repeat(["P", "Q"], 200)
    y = np.r_[rng.random(200) < 0.7, rng.random(200) < 0.2]
    score = (strata == "P").astype(float) + 1e-3 * rng.normal(size=400)   # protein identity only
    assert fast_auc(y, score) > 0.7
    assert abs(within_auc(y, score, strata) - 0.5) < 0.1


def test_calibrated_z_t_to_normal():
    """Few WT crystals -> heavier t tails -> the same delta/SE maps to a smaller |z|;
    many crystals -> identity; sign and order preserved."""
    assert abs(pairs.calibrated_z(2.0, 10 ** 6) - 2.0) < 1e-3
    small, big = pairs.calibrated_z(2.5, 3), pairs.calibrated_z(2.5, 30)
    assert small < big < 2.5
    assert pairs.calibrated_z(-2.5, 3) == -small
    zs = [pairs.calibrated_z(v, 5) for v in (0.5, 1.0, 2.0, 3.0, 8.0)]
    assert zs == sorted(zs)
    # tail probability preserved: P(|t_nu| > t) == P(|N| > z)
    from scipy.stats import norm, t as student_t
    nu = 3 - 1 + pairs.SHRINK_K
    assert abs(2 * norm.sf(small) - 2 * student_t.sf(2.5, nu)) < 1e-9


def test_circular_metric_wraps_at_180():
    """WT psi near +179 and a mutant at -179 differ by 2 degrees, not 358."""
    W = dict(per={10: (0.0, 1.0, 9, 0.0, 179.0)}, prior={10: 1.0})
    assert abs(pairs.rel(W, 10, -179.0) - 2.0) < 1e-9
    assert abs(pairs.rel(W, 10, 178.0) + 1.0) < 1e-9
    assert float(pairs.wrap(190.0)) == -170.0


def test_phi_psi_metric_matches_structure_features():
    phis = [-57.0, -57.0, -120.0, 60.0, -57.0, -65.0, -57.0]
    psis = [-47.0, -47.0, 130.0, 45.0, -47.0, 140.0, -47.0]
    d = tempfile.mkdtemp(); p = os.path.join(d, "bb.pdb")
    _write_backbone_pdb(p, _backbone(phis, psis), "AAGGAPA")
    st = pairs.parse_structure(p)
    st["chid"], st["res"] = "A", st["chains"]["A"]
    for r in (3, 4, 5):
        assert abs(pairs.METRICS["phi"][0](st, r - 2) - phis[r - 1]) < 0.5
        assert abs(pairs.METRICS["psi"][0](st, r - 2) - psis[r - 1]) < 0.5


def test_ncs_average():
    """Two identical chains: the averaged metric lies between the two copies."""
    d = tempfile.mkdtemp(); p = os.path.join(d, "ncs.pdb")
    ca = ideal_helix()
    write_pdb(p, SEQ, ca, "A", 1.8)
    ca2 = hinge(ca, 30, 10.0)
    lines = open(p).read().splitlines()
    extra = []
    for i, (aa, xyz) in enumerate(zip(SEQ, ca2), start=1):
        extra.append(f"ATOM  {i:5d}  CA  {AA3[aa]} B{i:4d}    "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{20.0:6.2f}           C")
    open(p, "w").write("\n".join(lines[:-1] + extra + ["END"]) + "\n")
    st = pairs.parse_structure(p)
    st["chid"], st["res"] = "A", st["chains"]["A"]
    a = pairs.window_bend(st["chains"]["A"], 28); b = pairs.window_bend(st["chains"]["B"], 28)
    assert pairs.ncs_copies(st) == ["B"]
    old = pairs.NCS_AVERAGE
    try:
        pairs.NCS_AVERAGE = True
        v = pairs.bend_of(st, 28)
    finally:
        pairs.NCS_AVERAGE = old
    assert abs(v - (a + b) / 2) < 1e-9 and abs(a - b) > 1.0


def test_enm_features_on_helix():
    import structure_features as sf
    ca = ideal_helix()
    res = {i + 1: ("A", ca[i], 20.0) for i in range(N_RES)}
    f = sf.enm_features(res, 30, 28)
    assert all(np.isfinite(v) for v in f.values())
    assert f["enm_msf_site"] > 0 and f["enm_bend_response"] > 0
    # chain ends fluctuate more than the middle in an elastic network
    end = sf.enm_features(res, 3, 1)["enm_msf_site"]
    assert end > f["enm_msf_site"]


def test_lattice_contacts_runs():
    import structure_features as sf
    d = tempfile.mkdtemp(); p = os.path.join(d, "x.pdb")
    write_pdb(p, SEQ, ideal_helix(), "A", 1.8)
    site, win = sf.lattice_contacts(p, "A", 30, 28)
    assert (np.isnan(site) and np.isnan(win)) or (0 <= site <= win)


def test_multi_label_combination(monkeypatch):
    """Combined label = largest |z| over MULTI_METRICS, rescaled so |z| > Z_MOVER
    exactly when max |z_m| > MULTI_T; only mutations labelled by every metric."""
    fake = {"bend": {(10, "A"): 1.0, (12, "G"): -2.0, (14, "W"): 0.5},
            "phi": {(10, "A"): -3.0, (12, "G"): 1.0, (14, "W"): 0.2},
            "psi": {(10, "A"): 0.1, (12, "G"): 2.5, (14, "W"): 2.0},
            "ca_tor_b": {(10, "A"): 0.3, (12, "G"): 0.4}}

    def single(paths, pick, protein, null_out=None, **kw):
        rows = [dict(r=r, mut=m, z=z, protein=protein) for (r, m), z in fake[pairs.METRIC].items()]
        return rows, dict(mutations=len(rows))
    monkeypatch.setattr(pairs, "build_pairs_single", single)
    monkeypatch.setattr(pairs, "METRIC", "multi")
    rows, diag = pairs.build_pairs({}, None, "P")
    by = {(x["r"], x["mut"]): x for x in rows}
    assert set(by) == {(10, "A"), (12, "G")}                    # (14, W) lacks ca_tor_b
    s = pairs.Z_MOVER / pairs.MULTI_T
    assert abs(by[(10, "A")]["z"] - (-3.0 * s)) < 1e-12 and by[(10, "A")]["z_metric"] == "phi"
    assert by[(10, "A")]["mover"] == (3.0 > pairs.MULTI_T)
    assert by[(12, "G")]["mover"] is False and abs(by[(12, "G")]["z"] - 2.5 * s) < 1e-12
    assert pairs.METRIC == "multi"


# ------------------------- protocol guarantees --------------------------------
def test_features_are_wt_only():
    """Rule: no model input may come from the mutant side. Perturb every mutant
    crystal (coordinates, resolution, temperature, altlocs) and require identical
    features for every mutation; only the label may change."""
    d1 = tempfile.mkdtemp(); paths = make_dataset(d1)
    rows_a, _ = pairs.build_pairs(paths, lambda c: ("A", c["A"]), "synthetic", min_cons=5)
    rng = np.random.default_rng(9)
    for pid, p in paths.items():
        st = pairs.parse_structure(p)
        seq = "".join(v[0] for _, v in sorted(st["chains"]["A"].items()))
        if seq == SEQ:
            continue                                      # WT crystal: untouched
        out = []
        for line in open(p).read().splitlines():
            if line.startswith("ATOM"):
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                xyz += rng.normal(0, 0.4, 3)
                line = line[:30] + "".join(f"{v:8.3f}" for v in xyz) + line[54:]
                out.append(line)
                out.append(line[:16] + "B" + line[17:])   # an alternate conformer
                continue
            if line.startswith("REMARK   2 RESOLUTION"):
                line = "REMARK   2 RESOLUTION.    2.40 ANGSTROMS."
                out.append(line)
                out.append("REMARK 200  TEMPERATURE           (KELVIN) : 293")
                continue
            out.append(line)
        open(p, "w").write("\n".join(out) + "\n")
    rows_b, _ = pairs.build_pairs(paths, lambda c: ("A", c["A"]), "synthetic", min_cons=5)
    A = {(x["r"], x["mut"]): x for x in rows_a}; B = {(x["r"], x["mut"]): x for x in rows_b}
    common = sorted(set(A) & set(B))
    assert len(common) >= 4
    cols = sorted({c for v in delta_model.FEATURE_SETS.values() for c in v})
    FA = delta_model.featurize([A[k] for k in common]); FB = delta_model.featurize([B[k] for k in common])
    for fa, fb in zip(FA, FB):
        for c in cols:
            va, vb = fa.get(c, np.nan), fb.get(c, np.nan)
            assert (np.isnan(va) and np.isnan(vb)) or va == vb, c
    assert any(A[k]["z"] != B[k]["z"] for k in common)    # the label did see the change


def test_no_family_crosses_folds():
    rng = np.random.default_rng(1)
    fam = np.array([f"F{rng.integers(0, 40)}" for _ in range(500)])
    for k in (5, 10):
        folds = delta_model.site_folds(fam, k, 0)
        for f in np.unique(fam):
            assert len(set(folds[fam == f])) == 1


def test_lockbox_disjoint_from_dev():
    import json
    import mine_pairs
    if not (os.path.exists(mine_pairs.MANIFEST) and os.path.exists(mine_pairs.LOCKBOX_MANIFEST)):
        import pytest
        pytest.skip("manifests not present")
    dev = json.load(open(mine_pairs.MANIFEST)); lb = json.load(open(mine_pairs.LOCKBOX_MANIFEST))
    dev_acc = set(dev["proteins"]) | set(dev["covered"])
    assert not dev_acc & set(lb["proteins"])
    dev_fam = {v["family"] for v in dev["proteins"].values()} | {v["family"] for v in dev["covered"].values()}
    assert not dev_fam & {v["family"] for v in lb["proteins"].values()}
    dev_seq = {v["reference"] for v in dev["proteins"].values()}
    assert not dev_seq & {v["reference"] for v in lb["proteins"].values()}
    dev_pdb = {p for v in dev["proteins"].values() for p in v["entries"]}
    assert not dev_pdb & {p for v in lb["proteins"].values() for p in v["entries"]}


def test_final_eval_shuffle_within_family_and_power():
    """Negative-control shuffles stay inside families; power output is sane."""
    import final_eval as fe
    rng = np.random.default_rng(0)
    fam = np.repeat(np.array(["a", "b", "c", "d"]), [5, 7, 3, 9])
    idx = fe.shuffle_within(fam, rng)
    assert sorted(idx) == list(range(len(fam)))
    assert (fam[idx] == fam).all()
    fam = np.repeat(np.arange(40).astype(str), 10)
    y = rng.random(400) < 0.3
    s = y + rng.normal(0, 1, 400)
    fe.N_BOOT = 200
    p = fe.power(y, s, rng.normal(0, 1, 400), fam)
    assert abs(p["n_eff_families"] - 40) < 1e-9
    assert 0.5 < p["min_detectable_auc"] < 0.75
    assert not fe.reportable(y[:20]) and fe.reportable(y)


def test_timelock_power_gate_formula():
    """Blinded gate: predicted MDA falls with n, matches the first lockbox's SE, guards n < 2."""
    import timelock as tl
    mda, se = tl.predicted_mda(111, "A")
    assert abs(se - 0.063) < 0.001                 # observed family-bootstrap SE was 0.066
    assert tl.predicted_mda(500, "A")[0] < mda
    assert tl.predicted_mda(0, "A")[0] == float("inf")
    assert tl.mutation_id(dict(protein="P1", wt="L", r=99, mut="A")) == "P1:L99A"


def test_circular_metric_wraps_across_180(monkeypatch):
    """WT angles straddling +-180 give a small sigma; a mutant at -177 gives a small delta."""
    vals = {"W1": 179.0, "W2": -179.0, "W3": 178.0, "W4": -178.0, "W5": 180.0}
    monkeypatch.setattr(pairs, "METRIC", "psi")
    monkeypatch.setattr(pairs, "bend_of", lambda st, s: vals[st["id"]])
    monkeypatch.setattr(pairs, "bz_of", lambda st: {})
    monkeypatch.setattr(pairs, "typical_hets", lambda structs, wts, s: None)
    structs = {p: {"id": p} for p in vals}
    W = pairs.form_stats(structs, sorted(vals), set(range(1, 10)), "pooled")
    med, sig, n, _, ref = W["per"][1]
    assert sig < 3.0 and n == 5 and ref is not None
    delta, se, z, _, _ = pairs.label(W, 1, [pairs.rel(W, 1, -177.0)])
    assert abs(delta) < 5.0
    # without the circular reference the same values would look ~180 deg apart
    assert np.std(list(vals.values())) > 100
