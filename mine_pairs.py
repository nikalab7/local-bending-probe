"""
proteinX -- SYSTEMATIC WT/mutant pair miner (SPEC 2.1, primary source).

The v2 labels came from five hand-picked mutational series (T4L + 4). With
four usable proteins, the site signal could not be separated from WT
noise structure. This script finds every protein in the PDB that can supply
clean pairs, using only RCSB metadata before downloading any coordinates:

  1. CANDIDATES  UniProt accessions with >= MIN_SINGLE X-ray (<= 2.5 A)
                 polymer entities carrying exactly one substitution against
                 UniProt (entity_poly.rcsb_mutation_count, from SIFTS).
  2. METADATA    per entity (GraphQL): canonical sequence, chain ids, space
                 group, cell, resolution, number of protein entities.
                 Only single-protein-entity entries are kept: a partner
                 protein next to the window is a confound like a ligand.
  3. PRE-SCREEN  per accession, pick the reference construct (the canonical
                 sequence with the best yield), then keep crystal forms
                 (pairs.assign_forms rules) with >= pairs.MIN_WT identical-
                 sequence WT crystals AND >= 1 Hamming-1 single mutant.
                 Only those crystals are downloaded.
  4. FAMILIES    RCSB sequence search (mmseqs2, >= 30% identity) on each
                 reference sequence; accessions that hit each other are
                 merged (union-find). delta_model holds out whole families.

Everything is pinned to manifests/mined.json on the first run (delete it to
refresh). Coordinates go to mined_pdb/ as .pdb.gz (gitignored).

Proteins already covered by the gate scripts (T4L + the Gate 4 set) are kept
in the family clustering but not re-mined, so no mutation is counted twice.

Usage: python mine_pairs.py [--min-single 8] [--workers 8]
"""
from __future__ import annotations
import argparse
import concurrent.futures as cf
import json
import os
import time
import urllib.request
from collections import Counter, defaultdict

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
MINED_DIR = "mined_pdb"
MANIFEST = os.path.join("manifests", "mined.json")
MIN_SINGLE = 8          # single-substitution entities per accession (facet cut)
MAX_LEN = 500           # residues; longer entities rarely have PDB-format files
MAX_WT_PER_FORM = 40    # plenty for a noise floor; caps e.g. CA II
FAMILY_IDENTITY = 0.3
ACC_ATTR = ("rcsb_polymer_entity_container_identifiers."
            "reference_sequence_identifiers.database_accession")
# already mined by feasibility_t4l.py / powered_loop_gate.py
COVERED = {"P00720": "T4L", "P00644": "SNase", "P00648": "barnase",
           "P61626": "human_lysozyme", "P61823": "RNaseA"}

_XRAY = [{"type": "terminal", "service": "text", "parameters": {
             "attribute": "exptl.method", "operator": "exact_match",
             "value": "X-RAY DIFFRACTION"}},
         {"type": "terminal", "service": "text", "parameters": {
             "attribute": "rcsb_entry_info.resolution_combined",
             "operator": "less_or_equal", "value": 2.5}}]

GQL = """query($ids:[String!]!){ polymer_entities(entity_ids:$ids){
  rcsb_id
  rcsb_polymer_entity{ pdbx_description }
  entity_poly{ rcsb_mutation_count pdbx_seq_one_letter_code_can }
  rcsb_polymer_entity_container_identifiers{ auth_asym_ids
    reference_sequence_identifiers{ database_name database_accession } }
  entry{ rcsb_entry_info{ resolution_combined polymer_entity_count_protein }
    symmetry{ space_group_name_H_M }
    cell{ length_a length_b length_c angle_alpha angle_beta angle_gamma } }
}}"""


# ------------------------------- HTTP ---------------------------------------
def _post(url, payload, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            return json.loads(body) if body else {}
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 ** (k + 1))


def candidate_accessions(min_single):
    """[(accession, n_single_entities)] from a search facet."""
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": _XRAY + [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "entity_poly.rcsb_mutation_count",
                "operator": "equals", "value": 1}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers."
                             "reference_sequence_identifiers.database_name",
                "operator": "exact_match", "value": "UniProt"}}]},
         "return_type": "polymer_entity",
         "request_options": {"paginate": {"start": 0, "rows": 0},
                             "results_content_type": ["experimental"],
                             "facets": [{"name": "acc", "aggregation_type": "terms",
                                         "attribute": ACC_ATTR,
                                         "max_num_intervals": 2000,
                                         "min_interval_population": min_single}]}}
    res = _post(SEARCH, q)
    return [(b["label"], b["population"]) for b in res["facets"][0]["buckets"]]


def entity_ids(acc):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": _XRAY + [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": ACC_ATTR, "operator": "exact_match", "value": acc}}]},
         "return_type": "polymer_entity",
         "request_options": {"return_all_hits": True,
                             "results_content_type": ["experimental"]}}
    res = _post(SEARCH, q)
    return [r["identifier"] for r in res.get("result_set", [])]


def entity_meta(ids, batch=200):
    out = []
    for i in range(0, len(ids), batch):
        res = _post(GRAPHQL, {"query": GQL, "variables": {"ids": ids[i:i + batch]}})
        out += [e for e in res["data"]["polymer_entities"] if e]
    return out


def family_hits(seq):
    """UniProt accessions of entities >= FAMILY_IDENTITY identical to seq."""
    q = {"query": {"type": "terminal", "service": "sequence", "parameters": {
            "evalue_cutoff": 1e-3, "identity_cutoff": FAMILY_IDENTITY,
            "sequence_type": "protein", "value": seq}},
         "return_type": "polymer_entity",
         "request_options": {"return_all_hits": True,
                             "results_content_type": ["experimental"]}}
    ids = [r["identifier"] for r in _post(SEARCH, q).get("result_set", [])]
    return ids


# ----------------------------- pre-screen -----------------------------------
def _flatten(e):
    """GraphQL entity -> flat record, or None if unusable."""
    try:
        info = e["entry"]["rcsb_entry_info"]
        cell = e["entry"]["cell"]; sym = e["entry"]["symmetry"]
        seq = e["entity_poly"]["pdbx_seq_one_letter_code_can"]
        chains = e["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]
        res = min(info["resolution_combined"] or [99.0])
        return dict(entity=e["rcsb_id"], pdb=e["rcsb_id"].split("_")[0],
                    chain=sorted(chains)[0], seq=seq, resolution=res,
                    n_protein_entities=info["polymer_entity_count_protein"],
                    spacegroup=sym["space_group_name_H_M"],
                    cell=tuple(cell[k] for k in ("length_a", "length_b", "length_c",
                                                 "angle_alpha", "angle_beta",
                                                 "angle_gamma")),
                    description=(e.get("rcsb_polymer_entity") or {}).get(
                        "pdbx_description") or "")
    except (KeyError, TypeError, ValueError):
        return None


def hamming1(a, b):
    """Position (0-based) if a and b differ at exactly one position, else None."""
    if len(a) != len(b):
        return None
    diff = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    return diff[0] if len(diff) == 1 else None


def prescreen(records, min_wt=None):
    """Choose the reference construct and the usable crystal forms.

    records: flat entity records of ONE accession (already resolution-cut).
    Returns dict(reference, forms={form: dict(wt=[...], single=[...])},
    n_single) or None. A form is usable with >= min_wt WT crystals of the
    reference sequence and >= 1 single mutant of it.
    """
    from pairs import assign_forms, MIN_WT
    min_wt = MIN_WT if min_wt is None else min_wt
    recs = [r for r in records if r["n_protein_entities"] == 1
            and r["spacegroup"] and all(c for c in r["cell"])
            and 30 <= len(r["seq"]) <= MAX_LEN]
    one_per_pdb = {}
    for r in recs:                           # one entity per entry
        one_per_pdb.setdefault(r["pdb"], r)
    recs = list(one_per_pdb.values())
    if not recs:
        return None
    form = assign_forms({r["pdb"]: dict(spacegroup=r["spacegroup"], cell=r["cell"])
                         for r in recs})
    by_seq = Counter(r["seq"] for r in recs)
    best = None
    for ref, cnt in by_seq.most_common(10):  # candidate reference constructs
        if cnt < min_wt:
            break
        forms = defaultdict(lambda: dict(wt=[], single=[]))
        for r in recs:
            f = form.get(r["pdb"])
            if f is None:
                continue
            if r["seq"] == ref:
                forms[f]["wt"].append(r)
            elif hamming1(r["seq"], ref) is not None:
                forms[f]["single"].append(r)
        forms = {f: v for f, v in forms.items() if len(v["wt"]) >= min_wt and v["single"]}
        n_single = len({(r["seq"]) for v in forms.values() for r in v["single"]})
        if n_single and (best is None or n_single > best["n_single"]):
            best = dict(reference=ref, forms=forms, n_single=n_single)
    if best is None:
        return None
    for v in best["forms"].values():         # best-resolution WT first, capped
        v["wt"] = sorted(v["wt"], key=lambda r: r["resolution"])[:MAX_WT_PER_FORM]
    return best


# ------------------------------ families ------------------------------------
def families(ref_seqs, workers):
    """{accession: family_id} by union-find over >= 30% identity hits."""
    parent = {a: a for a in ref_seqs}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    ent_to_acc = {}

    def hits(a):
        return a, family_hits(ref_seqs[a])

    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        hit = dict(ex.map(hits, ref_seqs))
    # map hit entities to accessions. Chimeric entities (e.g. T4L-GPCR
    # fusions) map to several UniProt accessions and would chain unrelated
    # families together, so only single-accession entities link.
    all_ids = sorted({i for v in hit.values() for i in v})
    acc_of = {}
    for m in entity_meta(all_ids):
        accs = {r["database_accession"] for r in
                (m["rcsb_polymer_entity_container_identifiers"]
                 ["reference_sequence_identifiers"] or [])
                if r["database_name"] == "UniProt"}
        if len(accs) == 1 and next(iter(accs)) in parent:
            acc_of[m["rcsb_id"]] = next(iter(accs))
    for a, ids in hit.items():
        for i in ids:
            b = acc_of.get(i)
            if b is not None and find(a) != find(b):
                parent[find(b)] = find(a)
    roots = {}
    return {a: roots.setdefault(find(a), f"F{len(roots):03d}") for a in sorted(parent)}


# ------------------------------ download ------------------------------------
def fetch_gz(pdb):
    path = os.path.join(MINED_DIR, f"{pdb}.pdb.gz")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return pdb
    for k in range(3):
        try:
            with urllib.request.urlopen(
                    f"https://files.rcsb.org/download/{pdb}.pdb.gz", timeout=60) as r:
                data = r.read()
            with open(path, "wb") as fh:
                fh.write(data)
            return pdb
        except urllib.error.HTTPError as e:
            if e.code == 404:                 # no PDB-format file (large entry)
                return None
        except Exception:
            pass
        time.sleep(2 ** (k + 1))
    return None


# -------------------------------- main --------------------------------------
def build_manifest(min_single, workers):
    cands = candidate_accessions(min_single)
    print(f"candidate accessions (>= {min_single} single-substitution entities): {len(cands)}")

    def screen(acc):
        try:
            recs = [r for r in map(_flatten, entity_meta(entity_ids(acc))) if r]
        except Exception as e:                # network failure: skip, report
            return acc, None, f"error: {e}"
        return acc, prescreen([r for r in recs if r["resolution"] <= 2.5]), ""

    proteins, why = {}, Counter()
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        for acc, ps, err in ex.map(screen, [a for a, _ in cands]):
            if err:
                why["query_error"] += 1
            elif ps is None:
                why["no_usable_form"] += 1
            else:
                proteins[acc] = ps
    print(f"accessions with a usable form: {len(proteins)}  (dropped: {dict(why)})")

    ref = {a: p["reference"] for a, p in proteins.items()}
    for acc in COVERED:                       # cluster the gate proteins too
        if acc not in ref:
            recs = [r for r in map(_flatten, entity_meta(entity_ids(acc))) if r]
            if recs:
                ref[acc] = Counter(r["seq"] for r in recs).most_common(1)[0][0]
    fam = families(ref, workers)

    out = dict(min_single=min_single, family_identity=FAMILY_IDENTITY,
               covered={a: dict(name=n, family=fam.get(a, a)) for a, n in COVERED.items()},
               proteins={})
    for acc, p in sorted(proteins.items()):
        if acc in COVERED:
            continue
        entries = {}
        for f, v in p["forms"].items():
            for kind in ("wt", "single"):
                for r in v[kind]:
                    entries[r["pdb"]] = dict(chain=r["chain"], form=f, kind=kind,
                                             resolution=r["resolution"])
        desc = Counter(r["description"] for v in p["forms"].values() for r in v["wt"])
        out["proteins"][acc] = dict(
            family=fam[acc], description=desc.most_common(1)[0][0],
            reference=p["reference"],
            reference_length=len(p["reference"]), n_single=p["n_single"],
            forms=sorted(p["forms"]), entries=entries)
    return out


def load_manifest():
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as fh:
            return json.load(fh)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-single", type=int, default=MIN_SINGLE)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    man = load_manifest()
    if man is None:
        man = build_manifest(args.min_single, args.workers)
        os.makedirs(os.path.dirname(MANIFEST), exist_ok=True)
        with open(MANIFEST, "w") as fh:
            json.dump(man, fh, indent=1, sort_keys=True)
        print(f"[manifest] pinned {len(man['proteins'])} proteins -> {MANIFEST}")
    else:
        print(f"[manifest] {len(man['proteins'])} pinned proteins from {MANIFEST}")
    os.makedirs(MINED_DIR, exist_ok=True)
    pdbs = sorted({p for v in man["proteins"].values() for p in v["entries"]})
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        got = [p for p in ex.map(fetch_gz, pdbs) if p]
    n_fam = len({v["family"] for v in man["proteins"].values()})
    print(f"downloaded/cached {len(got)}/{len(pdbs)} entries for "
          f"{len(man['proteins'])} proteins in {n_fam} families -> {MINED_DIR}/")


if __name__ == "__main__":
    main()
