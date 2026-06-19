"""
proteinX feasibility GO/NO-GO -- T4 lysozyme (UniProt P00720).

Question: do real point mutations produce a LOCAL Ca-bending change above the
crystal-to-crystal noise floor, often enough to make a benchmark exist?

Design choices that keep this honest:
  * Noise floor is PER-WINDOW, estimated from redundant WT* crystals. T4L has
    multiple crystal forms + a hinge-bending mode; using crystal-to-crystal WT
    scatter as the floor automatically discounts that, so we don't credit
    packing/hinge motion as a "mutation effect".
  * 2-sigma floor is deliberate: differencing noise is ~sqrt(2)*sigma ~ 1.41 sigma,
    so an event above 2*sigma also clears the differencing-noise band -> it is the
    subset a noise-limited ORACLE model could in principle resolve.
  * Baseline to beat = predict-zero, reported on the above-floor subset (where
    zero is a bad predictor), not the full set (where zero wins by true-negatives).
"""
import json, os, sys, urllib.request, concurrent.futures as cf
import numpy as np
from bending_metric import bending_angle, is_continuous

PDB_DIR = "t4l_pdb"; os.makedirs(PDB_DIR, exist_ok=True)
MAX_DL = 800
THREE2ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}


def fetch_id_list():
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": "P00720"}},
         "return_type": "entry",
         "request_options": {"return_all_hits": True}}
    req = urllib.request.Request(
        "https://search.rcsb.org/rcsbsearch/v2/query",
        data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    res = json.load(urllib.request.urlopen(req, timeout=40))
    ids = [r["identifier"] for r in res["result_set"]]
    print(f"RCSB search: {len(ids)} entries map to P00720")
    return ids


def fetch_pdb(pid):
    path = os.path.join(PDB_DIR, f"{pid}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return pid, True
    try:
        data = urllib.request.urlopen(
            f"https://files.rcsb.org/download/{pid}.pdb", timeout=30).read()
        with open(path, "wb") as f:
            f.write(data)
        return pid, True
    except Exception:
        return pid, False


def parse_ca(path):
    """Return {chainID: {resSeq: (one_letter, np.array([x,y,z]))}} for model 1,
    highest-occupancy altloc, standard residues only."""
    best = {}          # (chain,res) -> (occ, aa, xyz)
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            res3 = line[17:20].strip()
            aa = THREE2ONE.get(res3)
            if aa is None:
                continue
            if line[26] != " ":            # skip insertion codes
                continue
            chain = line[21]
            try:
                rs = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]),
                                float(line[46:54])])
                occ = float(line[54:60]) if line[54:60].strip() else 1.0
            except ValueError:
                continue
            key = (chain, rs)
            if key not in best or occ > best[key][0]:
                best[key] = (occ, aa, xyz)
    chains = {}
    for (chain, rs), (_, aa, xyz) in best.items():
        chains.setdefault(chain, {})[rs] = (aa, xyz)
    return chains


def pick_t4l_chain(chains):
    """Standalone full-length T4L: 150-175 std residues, max resSeq<=200
    (excludes GPCR/BRIL fusions and engineered insert/permutant chains)."""
    best, bestscore = None, 1e9
    for ch, res in chains.items():
        n = len(res)
        if n < 150 or n > 175 or max(res) > 200:
            continue
        score = abs(n - 164)
        if score < bestscore:
            best, bestscore = res, score
    return best


def window_bend(res, start):
    pts = []
    for r in range(start, start + 5):
        if r not in res:
            return None
        pts.append(res[r][1])
    pts = np.array(pts)
    if not is_continuous(pts):
        return None
    return bending_angle(pts)


def main():
    try:
        ids = fetch_id_list()
    except Exception as e:
        print(f"search API failed ({e}); aborting"); sys.exit(1)
    ids = ids[:MAX_DL]

    ok = []
    with cf.ThreadPoolExecutor(max_workers=12) as ex:
        for pid, good in ex.map(fetch_pdb, ids):
            if good:
                ok.append(pid)
    print(f"downloaded/cached {len(ok)}/{len(ids)} PDB files")

    # parse + pick T4L chain
    structs = {}
    for pid in ok:
        try:
            ch = pick_t4l_chain(parse_ca(os.path.join(PDB_DIR, f"{pid}.pdb")))
        except Exception:
            ch = None
        if ch:
            structs[pid] = ch
    print(f"usable standalone T4L chains: {len(structs)}")
    if len(structs) < 20:
        print("too few usable chains -- cannot assess feasibility"); sys.exit(1)

    # consensus WT* reference (majority residue per position, seen in >=10 structs)
    from collections import Counter, defaultdict
    counts = defaultdict(Counter)
    for res in structs.values():
        for rs, (aa, _) in res.items():
            counts[rs][aa] += 1
    consensus = {rs: c.most_common(1)[0][0]
                 for rs, c in counts.items() if sum(c.values()) >= 10}

    # mutation set per structure (vs consensus, over modeled consensus positions)
    mutsets = {}
    for pid, res in structs.items():
        muts = tuple(sorted(
            (rs, consensus[rs], aa) for rs, (aa, _) in res.items()
            if rs in consensus and aa != consensus[rs]))
        mutsets[pid] = muts

    wt_ids = [p for p, m in mutsets.items() if len(m) == 0]
    single = [(p, m[0]) for p, m in mutsets.items() if len(m) == 1]
    print(f"WT* (0-mutation) crystals: {len(wt_ids)}   single-mutant chains: {len(single)}")
    if len(wt_ids) < 3:
        print("not enough redundant WT crystals for a noise floor"); sys.exit(1)

    # per-window WT reference + robust sigma (1.4826*MAD)
    starts = range(min(consensus), max(consensus) - 4 + 1)
    wt_med, wt_sig = {}, {}
    for s in starts:
        vals = [window_bend(structs[p], s) for p in wt_ids]
        vals = np.array([v for v in vals if v is not None])
        if len(vals) >= 3:
            wt_med[s] = float(np.median(vals))
            wt_sig[s] = float(1.4826 * np.median(np.abs(vals - np.median(vals))))
    pooled_sig = float(np.median([v for v in wt_sig.values() if v > 0]))
    print(f"windows with WT floor: {len(wt_med)}   pooled per-window sigma = "
          f"{pooled_sig:.2f} deg")

    # single-mutation deltas at the window centered on the mutated residue
    rows = []   # (pid, resnum, delta, floor)
    for pid, (r, wtaa, mutaa) in single:
        s = r - 2                       # window r-2..r+2, mutation central
        if s not in wt_med:
            continue
        b = window_bend(structs[pid], s)
        if b is None:
            continue
        floor = wt_sig.get(s, pooled_sig) or pooled_sig
        rows.append((pid, r, b - wt_med[s], floor))
    if not rows:
        print("no scorable single-mutation windows"); sys.exit(1)

    D = np.array([x[2] for x in rows]); F = np.array([x[3] for x in rows])
    absd = np.abs(D)
    above2 = absd > 2 * F
    above3 = absd > 3 * F
    N = len(D); M2 = int(above2.sum()); M3 = int(above3.sum()); p = M2 / N

    # independent noise estimate: redundant crystals of the SAME variant
    var_groups = defaultdict(list)
    for pid, m in mutsets.items():
        var_groups[m].append(pid)
    rep_sigmas = []
    for m, pids in var_groups.items():
        if len(pids) < 2 or len(m) != 1:
            continue
        s = m[0][0] - 2
        vals = [window_bend(structs[p], s) for p in pids]
        vals = [v for v in vals if v is not None]
        if len(vals) >= 2:
            rep_sigmas.append(np.std(vals, ddof=1))
    rep_sig = float(np.median(rep_sigmas)) if rep_sigmas else float("nan")

    pz_full = float(np.sqrt(np.mean(D ** 2)))
    pz_above = float(np.sqrt(np.mean(D[above2] ** 2))) if M2 else 0.0
    diff_noise = np.sqrt(2) * pooled_sig

    print("\n" + "=" * 64 + "\n  FEASIBILITY REPORT (T4 lysozyme)\n" + "=" * 64)
    print(f"scorable single-mutation windows N      : {N}")
    print(f"median |delta|                          : {np.median(absd):.2f} deg")
    print(f"90th-pct |delta|                         : {np.percentile(absd,90):.2f} deg")
    print(f"per-window noise floor sigma (pooled)    : {pooled_sig:.2f} deg")
    print(f"redundant-variant sigma (independent)    : {rep_sig:.2f} deg")
    print(f"differencing noise ~ sqrt2*sigma         : {diff_noise:.2f} deg")
    print(f"above-floor events  >2sigma (M2)         : {M2}  ({100*p:.1f}% of N)")
    print(f"above-floor events  >3sigma (M3)         : {M3}  ({100*M3/N:.1f}% of N)")
    print(f"predict-zero RMSE  full set              : {pz_full:.2f} deg")
    print(f"predict-zero RMSE  above-floor subset    : {pz_above:.2f} deg "
          f"(this is what a model must beat)")
    print(f"median |delta| among above-floor         : "
          f"{np.median(absd[above2]) if M2 else float('nan'):.2f} deg")

    # confound check: are above-floor events spread along the sequence, or just
    # clustered at the T4L hinge (domain motion leaking through as "bending")?
    pos_above = [rows[i][1] for i in range(N) if above2[i]]
    hinge = set(range(9, 17)) | set(range(78, 84))
    frac_hinge = (sum(r in hinge for r in pos_above) / len(pos_above)
                  if pos_above else 0.0)
    print(f"\nabove-floor events: {len(set(pos_above))} distinct residue positions"
          f"  ({len(pos_above)} events); fraction in hinge band = {frac_hinge:.0%}")
    print(f"redundant single-mutant variant groups used for indep sigma: {len(rep_sigmas)}")

    # GO / NO-GO
    if M2 >= 50 and p >= 0.15:
        verdict = "GO -- build the pair-miner"
    elif M2 >= 10 or p >= 0.05:
        verdict = "PIVOT -- broaden target (general local geometry / loop sets) before building"
    else:
        verdict = "NO-GO as framed -- bending-from-point-mutation is near-null; rethink target"
    print("\nVERDICT:", verdict)
    print("Caveat: T4L is a rigid, mutation-tolerant fold -> this is a CONSERVATIVE\n"
          "        (low) estimate of the above-floor fraction vs the whole PDB.")

    # plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].hist(D, bins=40, color="#4C78A8", edgecolor="white")
    for k in (-2, 2):
        ax[0].axvline(k * pooled_sig, color="#E45756", ls="--", lw=1)
    ax[0].set_title("single-mutation local bending delta")
    ax[0].set_xlabel("delta bending (deg)  [dashed = +/-2 sigma floor]")
    ax[0].set_ylabel("count")
    ax[1].scatter(F, absd, s=12, alpha=0.5, color="#4C78A8")
    xx = np.linspace(0, F.max() * 1.05, 50)
    ax[1].plot(xx, 2 * xx, color="#E45756", ls="--", lw=1, label="2 sigma")
    ax[1].set_title(f"|delta| vs floor   (above 2sigma: {M2}/{N} = {100*p:.0f}%)")
    ax[1].set_xlabel("per-window floor sigma (deg)")
    ax[1].set_ylabel("|delta bending| (deg)")
    ax[1].legend()
    fig.tight_layout(); fig.savefig("feasibility_t4l.png", dpi=130)
    print("\nplot -> feasibility_t4l.png")

    import csv
    with open("feasibility_t4l.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["pdb", "resnum", "delta_deg", "floor_sigma_deg"])
        w.writerows(rows)
    print("rows -> feasibility_t4l.csv")


if __name__ == "__main__":
    main()
