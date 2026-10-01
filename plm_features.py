"""
proteinX -- protein-language-model features (ESM-2) for the delta model.

The substitution features in delta_model are physicochemical and either
context-free or tied to hand-picked mechanisms. A protein language model
scores a substitution IN ITS SEQUENCE CONTEXT, from evolution: the masked
marginal log p(mut) - log p(wt) at the site is a strong zero-shot predictor
of mutational effects on fitness and stability. If anything about "what"
the substitution is predicts backbone movement, this should see it.

Per mutation site (WT scaffold sequence, site masked) we store the 20
amino-acid log-probabilities; delta_model derives
  esm_llr      log p(mut) - log p(wt)          "what given where"
  esm_entropy  entropy of the site distribution "where" (conservation)
  esm_wt_logp  log p(wt)                        "where" (how expected WT is)

Note this is NOT local information: the model reads the whole sequence and
encodes evolutionary constraints learned from UniRef. It answers a
different question from the local-sequence probe: does sequence-level
evolutionary context carry the information?

Requires torch + fair-esm (not in requirements.txt; optional):
  pip install torch --index-url https://download.pytorch.org/whl/cpu
  pip install fair-esm
Usage: python plm_features.py   -> results/esm_site_logp.csv (committed)
"""
from __future__ import annotations
import csv
import os
import numpy as np

MODEL = "esm2_t30_150M_UR50D"
OUT = os.path.join("results", "esm_site_logp.csv")
AA = "ACDEFGHIKLMNPQRSTVWY"
MAX_GAP = 30           # missing residues filled with X up to this many per gap


def scaffold_sequence(res):
    """(sequence, {resSeq: index}) of a scaffold chain; gaps filled with X."""
    nums = sorted(res)
    seq, idx = [], {}
    for i, n in enumerate(nums):
        if i:
            gap = n - nums[i - 1] - 1
            seq.extend("X" * min(max(gap, 0), MAX_GAP))
        idx[n] = len(seq)
        seq.append(res[n][0])
    return "".join(seq), idx


def site_logp(rows, batch=8):
    """{(protein, r): {aa: log p}} from one masked ESM-2 pass per site."""
    import torch
    import esm
    model, alphabet = getattr(esm.pretrained, MODEL)()
    model.eval()
    conv = alphabet.get_batch_converter()
    aa_idx = [alphabet.get_idx(a) for a in AA]
    jobs = {}
    for x in rows:
        key = (x["protein"], x["r"])
        if key not in jobs:
            seq, idx = scaffold_sequence(x["scaffold_res"])
            if x["r"] in idx:
                jobs[key] = (seq, idx[x["r"]])
    out = {}
    keys = sorted(jobs)
    with torch.no_grad():
        for i in range(0, len(keys), batch):
            chunk = keys[i:i + batch]
            data = [(f"{k[0]}:{k[1]}", jobs[k][0]) for k in chunk]
            _, _, toks = conv(data)
            for j, k in enumerate(chunk):
                toks[j, 1 + jobs[k][1]] = alphabet.mask_idx      # +1: BOS token
            logits = model(toks)["logits"]
            for j, k in enumerate(chunk):
                lp = torch.log_softmax(logits[j, 1 + jobs[k][1], aa_idx], -1).numpy()
                out[k] = dict(zip(AA, map(float, lp)))
    return out


def load(path=OUT):
    """{(protein, r): {aa: log p}} from the committed CSV, or {} if absent."""
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            out[(row["protein"], int(row["r"]))] = {a: float(row[a]) for a in AA}
    return out


def features(lp, wt, mut):
    if lp is None or wt not in lp or mut not in lp:
        return dict(esm_llr=np.nan, esm_entropy=np.nan, esm_wt_logp=np.nan)
    p = np.exp(np.array([lp[a] for a in AA]))
    return dict(esm_llr=lp[mut] - lp[wt],
                esm_entropy=float(-(p * np.log(p + 1e-12)).sum()),
                esm_wt_logp=lp[wt])


def main():
    from pairs import load_proteins
    rows = [x for rs, _ in load_proteins().values() for x in rs]
    lp = site_logp(rows)
    os.makedirs("results", exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["protein", "r"] + list(AA))
        for (prot, r), v in sorted(lp.items()):
            w.writerow([prot, r] + [f"{v[a]:.4f}" for a in AA])
    print(f"wrote {OUT}: {len(lp)} sites ({MODEL})")


if __name__ == "__main__":
    main()
