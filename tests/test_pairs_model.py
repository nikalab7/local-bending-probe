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
            f = {c: rng.normal() for c in delta_model.SITE + delta_model.SUBST}
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
