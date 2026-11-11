#!/usr/bin/env python3
"""Revision-2 control analysis (deposited in release v1.1.0).

Inputs: the derived tables in ``data/`` (shipped with this repository) and the
four public CGGA release files (2020-05-06) that must be downloaded from
http://www.cgga.org.cn, exactly as for ``replication_cgga693.py``.

Environment variables (all optional):

    MRG4_CGGA      directory holding the four CGGA .txt files
                   (default: <repo>/raw/cgga)
    MRG4_DERIVED   directory holding the deposited derived CSV tables
                   (default: <repo>/data)
    MRG4_OUT       directory for the JSON/CSV this script writes
                   (default: <repo>/replication)

The machine-readable result of the run already performed for the article is in
``replication/`` and can be checked without re-downloading any cohort.
"""
import os

"""R47 main problem 3: pool the two validation cohorts (592 patients, same
consortium, same platform), split at random 200 times, and compute the C-index
difference between halves under the frozen coefficients. This gives the
"same population, same model, different half-sample" distribution that the 0.095
fall should be read against.

If the split-half fall is materially smaller than 0.095, the difference between
cohorts is not explained by sampling within one population.
"""
import csv, json, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
from lifelines.utils import concordance_index

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = ROOT
EX = os.environ.get("MRG4_CGGA", os.path.join(ROOT, "raw", "cgga"))
S = os.environ.get("MRG4_DERIVED", os.path.join(ROOT, "data"))
OUT = os.environ.get("MRG4_OUT", os.path.join(ROOT, "replication"))
G4 = ["BMP1", "KIF15", "TRAF3", "PDGFA"]
B = 200
rng = np.random.default_rng(20261111)

m25 = {r["Gene"]: float(r["Multivariate_coefficient_ridge_L2"])
       for r in csv.DictReader(open(S + "/Table_S12_ModelA_25基因系数.csv", encoding="utf-8-sig"))}
G25 = list(m25)


def read_genes(path, genes):
    want = set(genes)
    out = {}
    with open(path, encoding="utf-8", errors="ignore") as fh:
        ids = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] in want:
                out[p[0]] = np.array([float(x) if x not in ("", "NA") else np.nan for x in p[1:]])
    return ids, out


def clin(path):
    return {r["CGGA_ID"]: r for r in csv.DictReader(
        open(path, encoding="utf-8-sig", errors="ignore"), delimiter="\t")}


def cohort(rsem, clinical, genes):
    ids, expr = read_genes(rsem, genes)
    cl = clin(clinical)
    keep = [i for i, pid in enumerate(ids)
            if pid in cl and cl[pid]["Grade"] in ("WHO II", "WHO III")
            and cl[pid]["OS"] not in ("", "NA")
            and cl[pid]["Censor (alive=0; dead=1)"] in ("0", "1")]
    t = np.array([float(cl[ids[i]]["OS"]) / 30.4375 for i in keep])
    e = np.array([float(cl[ids[i]]["Censor (alive=0; dead=1)"]) for i in keep], int)
    Z = {}
    for g in genes:
        v = np.log(np.asarray(expr[g], float) + 1.0).copy()
        v = (v - np.nanmean(v)) / np.nanstd(v, ddof=1)
        Z[g] = v[keep]
    return t, e, Z


c325 = cohort(EX + "/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",
              EX + "/CGGA.mRNAseq_325_clinical.20200506.txt", G4 + G25)
c693 = cohort(EX + "/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",
              EX + "/CGGA.mRNAseq_693_clinical.20200506.txt", G4 + G25)

# pool, keeping each cohort's own z-reference (the convention the frozen scores use)
t = np.concatenate([c325[0], c693[0]])
e = np.concatenate([c325[1], c693[1]])
cohort_id = np.array([0] * len(c325[0]) + [1] * len(c693[0]))
Z = {g: np.concatenate([c325[2][g], c693[2][g]]) for g in G4 + G25}
print(f"pooled n = {len(t)}, events = {int(e.sum())}")

# frozen coefficients: Model A from the deposited ridge table; MRG-4 as published
b25 = np.array([m25[g] for g in G25])
mr4 = {}
for r in csv.DictReader(open(S + "/Table_S17_MRG4部署常量.csv", encoding="utf-8-sig")):
    if r["item"] == "系数":
        mr4[r["detail"].split(" ")[0]] = float(r["value"])
b4 = np.array([mr4[g] for g in G4])
print("MRG-4 coefficients:", {g: round(v, 4) for g, v in zip(G4, b4)})


def cd(genes, beta, idx):
    X = np.column_stack([Z[g][idx] for g in genes])
    lp = X @ beta
    ok = np.isfinite(lp)
    return concordance_index(t[idx][ok], -lp[ok], e[idx][ok])


n = len(t)
obs4, obs25 = cd(G4, b4, np.arange(n)), cd(G25, b25, np.arange(n))
print(f"pooled C (frozen): MRG-4 {obs4:.4f}, 25-gene {obs25:.4f}")

d4, d25, bw = [], [], []
for b in range(B):
    perm = rng.permutation(n)
    h1, h2 = perm[: n // 2], perm[n // 2:]
    a4, a25 = cd(G4, b4, h1), cd(G25, b25, h1)
    c4, c25 = cd(G4, b4, h2), cd(G25, b25, h2)
    d4.append(a4 - c4)
    d25.append(a25 - c25)
    bw.append((a4 - a25) - (c4 - c25))

res = dict(
    n=int(n), events=int(e.sum()), reps=B, seed=20261111,
    pooled_c4=float(obs4), pooled_c25=float(obs25),
    fall_c4=dict(q=[float(np.percentile(np.abs(d4), 2.5)), float(np.percentile(np.abs(d4), 97.5))],
                 p95=float(np.percentile(np.abs(d4), 95)),
                 share_ge_0095=float((np.abs(d4) >= 0.095).mean())),
    fall_c25=dict(q=[float(np.percentile(np.abs(d25), 2.5)), float(np.percentile(np.abs(d25), 97.5))],
                  p95=float(np.percentile(np.abs(d25), 95)),
                  share_ge_0016=float((np.abs(d25) >= 0.016).mean())),
    between_model=dict(q=[float(np.percentile(bw, 2.5)), float(np.percentile(bw, 97.5))],
                       p95=float(np.percentile(np.abs(bw), 95))))
json.dump(res, open(OUT + "/analysis_splithalf.json", "w"), indent=2)
print(f"|fall| MRG-4: 2.5-97.5 pct {res['fall_c4']['q'][0]:.3f}-{res['fall_c4']['q'][1]:.3f}, "
      f"95th pct {res['fall_c4']['p95']:.3f}, share >= 0.095: {res['fall_c4']['share_ge_0095']:.1%}")
print(f"|fall| 25-gene: 2.5-97.5 pct {res['fall_c25']['q'][0]:.3f}-{res['fall_c25']['q'][1]:.3f}, "
      f"95th pct {res['fall_c25']['p95']:.3f}, share >= 0.016: {res['fall_c25']['share_ge_0016']:.1%}")
print(f"between-model diff across halves: {res['between_model']['q'][0]:+.3f} to "
      f"{res['between_model']['q'][1]:+.3f}")
print("written analysis_splithalf.json")
