"""
proteinX -- shared statistics + reproducibility helpers.

Why this exists
---------------
1. CLUSTER BOOTSTRAP. Validation pairs are not independent: several mutations
   hit the same residue (72 T4L movers sit on only 38 residues), and pooled
   loop sets are dominated by one protein. Resampling individual pairs treats
   correlated rows as independent and makes CIs too narrow. We resample whole
   clusters (residue, or protein+residue) instead.
2. PAIRED DELTA-AUC. "Model B beats model A" is a claim about the DIFFERENCE,
   not about B's own CI. Both models are scored on the SAME bootstrap
   resample, and the CI is reported on AUC_B - AUC_A.
3. PINNED ID MANIFESTS. RCSB search results change as the PDB grows, so live
   queries make runs drift even with fixed seeds. The first run writes the
   ID list to manifests/<name>.json; later runs read it instead of querying.
"""
from __future__ import annotations
import json
import os
import numpy as np

N_BOOT = 3000
MANIFEST_DIR = "manifests"


# ------------------------------ bootstrap -----------------------------------
def _cluster_index(clusters):
    """Map cluster labels -> list of row-index arrays."""
    clusters = np.asarray([str(c) for c in clusters])
    labels = np.unique(clusters)
    return [np.flatnonzero(clusters == lab) for lab in labels]


def _resamples(clusters, n_boot, seed):
    """Yield row-index arrays, each a resample of whole clusters."""
    groups = _cluster_index(clusters)
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        yield np.concatenate([groups[g] for g in pick])


def cluster_auc_ci(y, score, clusters, n_boot=N_BOOT, seed=0, level=0.90):
    """ROC-AUC + cluster-bootstrap CI. Returns (auc, lo, hi, n_clusters)."""
    from sklearn.metrics import roc_auc_score
    y = np.asarray(y, bool); score = np.asarray(score, float)
    if not (0 < y.sum() < len(y)):
        return float("nan"), float("nan"), float("nan"), 0
    auc = roc_auc_score(y, score)
    b = []
    for bi in _resamples(clusters, n_boot, seed):
        if 0 < y[bi].sum() < len(bi):
            b.append(roc_auc_score(y[bi], score[bi]))
    a = (1 - level) / 2 * 100
    lo, hi = np.percentile(b, [a, 100 - a])
    return float(auc), float(lo), float(hi), len(_cluster_index(clusters))


def paired_delta_auc(y, score_a, score_b, clusters, n_boot=N_BOOT, seed=0,
                     level=0.90):
    """AUC_b - AUC_a with a paired cluster-bootstrap CI.

    Returns (delta, lo, hi, p_one_sided) where p is the bootstrap fraction of
    resamples with delta <= 0 (i.e. evidence that b is NOT better than a).
    """
    from sklearn.metrics import roc_auc_score
    y = np.asarray(y, bool)
    sa = np.asarray(score_a, float); sb = np.asarray(score_b, float)
    if not (0 < y.sum() < len(y)):
        return float("nan"), float("nan"), float("nan"), float("nan")
    delta = roc_auc_score(y, sb) - roc_auc_score(y, sa)
    d = []
    for bi in _resamples(clusters, n_boot, seed):
        if 0 < y[bi].sum() < len(bi):
            d.append(roc_auc_score(y[bi], sb[bi]) - roc_auc_score(y[bi], sa[bi]))
    d = np.array(d)
    a = (1 - level) / 2 * 100
    lo, hi = np.percentile(d, [a, 100 - a])
    return float(delta), float(lo), float(hi), float((d <= 0).mean())


# ------------------------------ manifests -----------------------------------
def pinned_ids(name, fetch_fn):
    """Return a pinned PDB-ID list; query RCSB (fetch_fn) only on first run.

    Delete manifests/<name>.json to deliberately refresh against today's PDB.
    """
    os.makedirs(MANIFEST_DIR, exist_ok=True)
    path = os.path.join(MANIFEST_DIR, f"{name}.json")
    if os.path.exists(path):
        with open(path) as fh:
            ids = json.load(fh)["ids"]
        print(f"[manifest] {name}: {len(ids)} pinned ids from {path}")
        return ids
    ids = list(fetch_fn())
    if ids:                                   # never pin an empty/failed query
        with open(path, "w") as fh:
            json.dump({"name": name, "ids": ids}, fh, indent=0)
        print(f"[manifest] {name}: pinned {len(ids)} ids -> {path}")
    return ids
