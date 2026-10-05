"""
proteinX -- full-atom WT site context and substitution x context features.

delta_model's original descriptors see only C-alpha atoms, and its
substitution features are context-free (a volume change counts the same in
a helix core and on a surface loop). That makes "the substitution adds
nothing" a weak null: the physically known ways a point mutation bends a
backbone are all substitution x local-structure interactions:

  * Gly in a positive-phi (left-handed) conformation replaced by a non-Gly
  * Pro introduced where phi is far from -65 deg, or inside a helix
  * a larger side chain forced into a packed site (overpacking)
  * a side chain that H-bonds the backbone (capping box, Asx/ST turns)
    replaced by one that cannot
  * helix / strand propensity changes inside that secondary structure

This module parses the WT scaffold with all heavy atoms and computes:

  CONTEXT  ("where", full-atom): phi/psi region of the residue, backbone
           H-bonds, side-chain heavy-atom contacts, C-beta density, side-chain
           to backbone H-bonds, distance to the nearest het group, position
           inside the SS element.
  INTERACT ("what given where"): the mechanisms above, each zero unless the
           substitution and the local structure both fit it.

Geometry is distance-based (no hydrogens): an H-bond is a donor/acceptor
N/O pair within HB_CUT. Good enough for counting, not for energetics.
"""
from __future__ import annotations
import gzip
import numpy as np

HB_CUT = 3.5           # Angstrom, N/O - N/O
CONTACT_CUT = 4.5      # Angstrom, heavy atom - heavy atom
CB_CUT = 10.0          # Angstrom, C-beta density radius
LIG_CAP = 20.0         # Angstrom, distance-to-ligand cap
WATER = {"HOH", "DOD", "WAT", "H2O"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
SC_POLAR = {"OD1", "OD2", "ND2", "OG", "OG1", "NE2", "OE1", "OE2", "NZ", "NH1",
            "NH2", "NE", "OH", "ND1", "NE1", "SG"}
CAN_HBOND_SC = set("DENQSTKRHYWC")          # residues with polar side-chain atoms
CHARGED = set("DEKR")
THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'}
# helix propensity, kcal/mol relative to Ala (Pace & Scholtz 1998)
HELIX_DDG = dict(A=0.00, L=0.21, R=0.21, M=0.24, K=0.26, Q=0.39, E=0.40, I=0.41,
                 W=0.49, S=0.50, Y=0.53, F=0.54, H=0.61, V=0.61, N=0.65, T=0.66,
                 C=0.68, D=0.69, G=1.00, P=3.16)
# beta-strand propensity (Chou & Fasman 1978)
SHEET_P = dict(V=1.70, I=1.60, Y=1.47, F=1.38, W=1.37, L=1.30, C=1.19, T=1.19,
               Q=1.10, M=1.05, R=0.93, N=0.89, H=0.87, A=0.83, S=0.75, G=0.75,
               K=0.74, P=0.55, D=0.54, E=0.37)
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4,
           G=60.1, H=153.2, I=166.7, L=166.7, K=168.6, M=162.9, F=189.9,
           P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)

CONTEXT = ["phi_pos", "phi_alpha", "phi_beta", "bb_hbonds", "sc_contacts",
           "cb_density", "sc_bb_hbonds", "log_lig_dist", "ss_end_dist",
           "helix_nterm", "helix_cterm"]
LATTICE = ["lattice_site", "lattice_window"]   # crystal-packing contacts (gemmi)
INTERACT = ["gly_phi_pos_loss", "pro_phi_strain", "pro_in_helix", "overpack",
            "cavity_cb", "sc_hb_loss", "helix_prop_change", "sheet_prop_change",
            "buried_charge_in"]


# --------------------------------- parsing ----------------------------------
def parse_atoms(path, chain):
    """Heavy atoms of one chain + non-water het atoms (model 1, best altloc).

    Returns (residues, het_xyz): residues = {resSeq: (aa, {atom: xyz})} for
    standard residues without insertion codes; het_xyz = (k, 3) array.
    """
    opener = gzip.open if path.endswith(".gz") else open
    best, het = {}, []
    with opener(path, "rt") as fh:
        for line in fh:
            rec = line[:6]
            if rec.startswith("ENDMDL"):
                break
            if rec == "HETATM":
                if line[17:20].strip() in WATER or line[17:20].strip() in THREE2ONE:
                    continue
                try:
                    het.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    pass
                continue
            if rec != "ATOM  " or line[21] != chain or line[26] != " ":
                continue
            aa = THREE2ONE.get(line[17:20].strip())
            name = line[12:16].strip()
            if aa is None or name.startswith("H") or (line[76:78].strip() in ("H", "D")):
                continue
            try:
                rs = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                occ = float(line[54:60]) if line[54:60].strip() else 1.0
            except ValueError:
                continue
            key = (rs, name)
            if key not in best or occ > best[key][0]:
                best[key] = (occ, aa, xyz)
    residues = {}
    for (rs, name), (_, aa, xyz) in best.items():
        residues.setdefault(rs, (aa, {}))[1][name] = xyz
    return residues, (np.array(het) if het else np.zeros((0, 3)))


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def phi_psi(residues, r):
    """(phi, psi) of residue r in degrees; NaN where atoms are missing."""
    def at(k, a):
        return residues.get(k, (None, {}))[1].get(a)
    phi = psi = np.nan
    a = [at(r - 1, "C"), at(r, "N"), at(r, "CA"), at(r, "C")]
    if all(x is not None for x in a) and np.linalg.norm(a[0] - a[1]) < 2.0:
        phi = dihedral(*a)
    b = [at(r, "N"), at(r, "CA"), at(r, "C"), at(r + 1, "N")]
    if all(x is not None for x in b) and np.linalg.norm(b[2] - b[3]) < 2.0:
        psi = dihedral(*b)
    return phi, psi


# --------------------------------- features ---------------------------------
def context_features(residues, het_xyz, r, ss, ss_span):
    """Full-atom "where" descriptors of residue r.

    ss_span: (start, end) of the SS element containing r, or None for loops.
    """
    aa, atoms = residues[r]
    phi, psi = phi_psi(residues, r)
    names = []
    xyz = []
    for k, (_, at) in residues.items():
        for n, x in at.items():
            names.append((k, n)); xyz.append(x)
    xyz = np.array(xyz)
    own = np.array([k == r for k, _ in names])
    # backbone H-bonds of residue r (N or O to backbone O/N of residues >= 2
    # apart; the i +- 1 peptide N...O contact at ~2.3 A is covalent geometry)
    bb_hb = 0
    for a in ("N", "O"):
        if a in atoms:
            partner = "O" if a == "N" else "N"
            sel = np.array([(n == partner) and abs(k - r) >= 2 for k, n in names])
            if sel.any():
                bb_hb += int((np.linalg.norm(xyz[sel] - atoms[a], axis=1) < HB_CUT).sum())
    # side-chain contacts and side-chain -> backbone H-bonds
    sc = {n: x for n, x in atoms.items() if n not in BACKBONE}
    sc_contacts = sc_bb_hb = 0
    if sc:
        scx = np.array(list(sc.values()))
        d = np.linalg.norm(xyz[~own][:, None, :] - scx[None, :, :], axis=2)
        sc_contacts = int((d.min(axis=1) < CONTACT_CUT).sum())
        polar = np.array([n in SC_POLAR for n in sc])
        if polar.any():
            bbsel = np.array([(n in ("N", "O")) and (k != r) for k, n in names])
            if bbsel.any():
                dp = np.linalg.norm(xyz[bbsel][:, None, :] - scx[polar][None, :, :], axis=2)
                sc_bb_hb = int((dp.min(axis=1) < HB_CUT).sum())
    # C-beta density (C-alpha for Gly)
    def cb(at):
        return at.get("CB", at.get("CA"))
    me = cb(atoms)
    cbs = np.array([cb(at) for k, (_, at) in residues.items() if k != r and cb(at) is not None])
    cb_density = int((np.linalg.norm(cbs - me, axis=1) < CB_CUT).sum()) if me is not None and len(cbs) else 0
    # distance to the nearest het atom
    allr = np.array(list(atoms.values()))
    lig = float(np.min(np.linalg.norm(het_xyz[:, None, :] - allr[None, :, :], axis=2))) \
        if len(het_xyz) else LIG_CAP
    lig = min(lig, LIG_CAP)
    # position inside the SS element
    if ss_span is None:
        ss_end, nterm, cterm = 0.0, 0.0, 0.0
    else:
        a0, a1 = ss_span
        ss_end = float(min(r - a0, a1 - r))
        nterm = float(ss == "H" and r - a0 < 3)
        cterm = float(ss == "H" and a1 - r < 3)
    return dict(
        phi=phi, psi=psi,
        phi_pos=float(np.isfinite(phi) and phi > 0),
        phi_alpha=float(np.isfinite(phi) and phi < 0 and np.isfinite(psi) and -120 < psi < 50),
        phi_beta=float(np.isfinite(phi) and phi < 0 and np.isfinite(psi) and not -120 < psi < 50),
        bb_hbonds=float(bb_hb), sc_contacts=float(sc_contacts), cb_density=float(cb_density),
        sc_bb_hbonds=float(sc_bb_hb), log_lig_dist=float(np.log1p(lig)),
        ss_end_dist=ss_end, helix_nterm=nterm, helix_cterm=cterm)


def interaction_features(ctx, wt, mut, ss):
    """"What given where": zero unless substitution and local structure both fit."""
    burial = ctx["cb_density"] / 20.0
    dv = (VOL[mut] - VOL[wt]) / 100.0
    phi = ctx["phi"]
    pro_dev = 0.0
    if mut == "P" and np.isfinite(phi):
        pro_dev = min(abs(phi - (-65.0)) / 60.0, 2.0)
    return dict(
        gly_phi_pos_loss=float(wt == "G" and mut != "G" and ctx["phi_pos"] > 0),
        pro_phi_strain=pro_dev,
        pro_in_helix=float(mut == "P" and ss == "H" and ctx["helix_nterm"] == 0),
        overpack=max(0.0, dv) * burial,
        cavity_cb=max(0.0, -dv) * burial,
        sc_hb_loss=float(ctx["sc_bb_hbonds"]) * float(mut not in CAN_HBOND_SC),
        helix_prop_change=(HELIX_DDG[mut] - HELIX_DDG[wt]) * float(ss == "H"),
        sheet_prop_change=(SHEET_P[mut] - SHEET_P[wt]) * float(ss == "E"),
        buried_charge_in=float(mut in CHARGED and wt not in CHARGED) * burial)


def ss_span_of(helix, sheet, chain, r):
    """(start, end) of the contiguous HELIX/SHEET run containing (chain, r)."""
    for recs in (helix, sheet):
        if (chain, r) in recs:
            a = r
            while (chain, a - 1) in recs:
                a -= 1
            b = r
            while (chain, b + 1) in recs:
                b += 1
            return a, b
    return None


_CACHE = {}
_LATTICE = {}
_ENM = {}
ENM = ["enm_msf_site", "enm_bend_response", "enm_bend_response_rel"]
ENM_CUTOFF = 15.0      # Angstrom, ANM spring cutoff between C-alphas


def _anm_pinv(ca):
    """Pseudo-inverse of the anisotropic-network Hessian (6 rigid-body modes removed)."""
    n = len(ca)
    d = ca[:, None, :] - ca[None, :, :]
    r2 = (d ** 2).sum(-1)
    contact = (r2 < ENM_CUTOFF ** 2) & ~np.eye(n, dtype=bool)
    H = np.zeros((3 * n, 3 * n))
    for i, j in zip(*np.nonzero(np.triu(contact))):
        e = d[i, j] / np.sqrt(r2[i, j])
        block = -np.outer(e, e)
        H[3*i:3*i+3, 3*j:3*j+3] = block
        H[3*j:3*j+3, 3*i:3*i+3] = block
        H[3*i:3*i+3, 3*i:3*i+3] -= block
        H[3*j:3*j+3, 3*j:3*j+3] -= block
    w, v = np.linalg.eigh(H)
    keep = w > 1e-6 * max(w.max(), 1e-12)
    keep[:6] = False
    return (v[:, keep] / w[keep]) @ v[:, keep].T


def enm_features(res, r, s, n_dir=12):
    """Mechanics of the WT scaffold (C-alpha anisotropic network model).

    enm_msf_site          mean-square fluctuation of residue r, relative to the chain mean
    enm_bend_response     |change of the window bend angle| per unit force applied at r,
                          averaged over n_dir directions (perturbation response)
    enm_bend_response_rel the same, relative to the median over all windows of the chain
    """
    from bending_metric import bending_angle
    key = id(res)
    nums = sorted(res)
    if key not in _ENM:
        ca = np.array([res[k][1] for k in nums])
        _ENM[key] = (_anm_pinv(ca), ca, {k: i for i, k in enumerate(nums)}, {})
    G, ca, idx, resp_cache = _ENM[key]
    if r not in idx or any(w not in idx for w in range(s, s + 5)):
        return {k: np.nan for k in ENM}
    msf = np.array([np.trace(G[3*i:3*i+3, 3*i:3*i+3]) for i in range(len(nums))])
    dirs = np.random.default_rng(0).normal(size=(n_dir, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)

    def response(rr, ss_):
        if (rr, ss_) in resp_cache:
            return resp_cache[(rr, ss_)]
        i = idx[rr]; win = [idx[w] for w in range(ss_, ss_ + 5)]
        b0 = bending_angle(ca[win]); out = []
        for u in dirs:
            dx = (G[:, 3*i:3*i+3] @ u).reshape(-1, 3)
            eps = 0.1 / max(np.abs(dx[win]).max(), 1e-9)       # small, linear-regime step
            out.append(abs(bending_angle(ca[win] + eps * dx[win]) - b0) / eps)
        resp_cache[(rr, ss_)] = float(np.mean(out))
        return resp_cache[(rr, ss_)]
    resp = response(r, s)
    ref = [response(k + 2, k) for k in nums[::5] if k + 2 in idx and all(w in idx for w in range(k, k + 5))]
    return dict(enm_msf_site=float(msf[idx[r]] / msf.mean()),
                enm_bend_response=resp,
                enm_bend_response_rel=float(resp / np.median(ref)) if ref else np.nan)


def lattice_contacts(path, chain, r, s, cutoff=4.0):
    """Crystal-packing contacts of the WT scaffold (needs gemmi; NaN without it).

    Returns (atoms of residue r, atoms of window s..s+4) within `cutoff` of an
    atom of a symmetry mate (another copy in the lattice, image_idx != 0) or of
    another chain of the asymmetric unit. These are packing contacts, not
    contacts inside the molecule.
    """
    try:
        import gemmi
    except ImportError:
        return np.nan, np.nan
    if path not in _LATTICE:
        st = gemmi.read_structure(path)
        st.remove_hydrogens()
        model = st[0]
        ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
        _LATTICE[path] = (model, ns)
    model, ns = _LATTICE[path]
    try:
        ch = model[chain]
    except Exception:
        return np.nan, np.nan
    site = win = 0
    for res in ch:
        if res.het_flag == "H" or not s <= res.seqid.num <= s + 4 or res.seqid.icode != " ":
            continue
        for atom in res:
            for m in ns.find_atoms(atom.pos, "\0", radius=cutoff):
                cra = m.to_cra(model)
                if m.image_idx != 0 or cra.chain.name != chain:
                    if cra.residue.het_flag != "H":          # protein neighbour, not water/ligand
                        win += 1
                        site += res.seqid.num == r
                        break
    return float(site), float(win)


def row_features(x):
    """CONTEXT + INTERACT features for one pairs row (needs scaffold_path/chain)."""
    key = (x["scaffold_path"], x["scaffold_chain"])
    if key not in _CACHE:
        _CACHE[key] = parse_atoms(*key)
    residues, het = _CACHE[key]
    nan = {k: np.nan for k in CONTEXT + INTERACT}
    if x["r"] not in residues or "CA" not in residues[x["r"]][1]:
        return nan
    span = ss_span_of(x["scaffold_helix"], x["scaffold_sheet"], x["scaffold_chain"], x["r"])
    ctx = context_features(residues, het, x["r"], x["scaffold_ss"], span)
    out = {k: ctx[k] for k in CONTEXT}
    out.update(interaction_features(ctx, x["wt"], x["mut"], x["scaffold_ss"]))
    out["lattice_site"], out["lattice_window"] = lattice_contacts(
        x["scaffold_path"], x["scaffold_chain"], x["r"], x["s"])
    out.update(enm_features(x["scaffold_res"], x["r"], x["s"]))
    return out
