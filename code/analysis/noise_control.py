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

"""R66 main 1 and main 9 -- two controls that decide whether v13's argument holds.

(1) Is the refit gain of the 25-gene index optimism of the parameter count, or
    information in the extra genes? Control: a 25-parameter model made of the
    four MRG-4 genes plus 21 permuted noise columns. Same number of parameters,
    no extra information about the outcome. If the noise model gains as much as
    the real 25-gene index, the gain is the parameter count.

(2) The pooled split-half reference of v13 mixes the two cohorts, so part of its
    spread is the between-cohort difference itself. Repeat the split inside each
    cohort: that is within-cohort sampling only.
"""
import csv, json, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = ROOT
EX = os.environ.get("MRG4_CGGA", os.path.join(ROOT, "raw", "cgga"))
S = os.environ.get("MRG4_DERIVED", os.path.join(ROOT, "data"))
OUT = os.environ.get("MRG4_OUT", os.path.join(ROOT, "replication"))
G4 = ["BMP1", "KIF15", "TRAF3", "PDGFA"]
B = 300
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


def fit(X, t, e, idx):
    d = pd.DataFrame(X[idx], columns=[f"g{i}" for i in range(X.shape[1])])
    d["T"] = t[idx]
    d["E"] = e[idx]
    d = d[d.notna().all(axis=1)].reset_index(drop=True)
    if len(d) < 60 or d["E"].sum() < 18:
        return None
    try:
        return CoxPHFitter(penalizer=0.0).fit(d, "T", "E").params_.values
    except Exception:
        return None


def score(X, beta, t, e, idx):
    lp = X[idx] @ beta
    ok = np.isfinite(lp)
    if ok.sum() < 40 or e[idx][ok].sum() < 12:
        return np.nan
    return concordance_index(t[idx][ok], -lp[ok], e[idx][ok])


c325 = cohort(EX + "/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",
              EX + "/CGGA.mRNAseq_325_clinical.20200506.txt", G4 + G25)
c693 = cohort(EX + "/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",
              EX + "/CGGA.mRNAseq_693_clinical.20200506.txt", G4 + G25)

out = {}
print("[1] parameter-matched noise control: 4 real genes + 21 permuted columns")
for lab, c in (("CGGA-325", c325), ("CGGA-693", c693)):
    t, e, Z = c
    n = len(t)
    X4 = np.column_stack([Z[g] for g in G4])
    X25 = np.column_stack([Z[g] for g in G25])
    frozen4 = concordance_index(t, -X4 @ np.array([0.2582, 0.4636, -0.2379, 0.4795]), e)
    g4, g25, gnoise, g25p = [], [], [], []
    allidx = np.arange(n)
    for b in range(B):
        ii = rng.integers(0, n, n)
        noise = np.column_stack([rng.permutation(X4[:, k % 4]) for k in range(21)])
        Xn = np.column_stack([X4, noise])
        b4 = fit(X4, t, e, ii)
        b25 = fit(X25, t, e, ii)
        bn = fit(Xn, t, e, ii)
        if b4 is None or b25 is None or bn is None:
            continue
        c4 = score(X4, b4, t, e, allidx)
        c25 = score(X25, b25, t, e, allidx)
        cn = score(Xn, bn, t, e, allidx)
        if any(np.isnan(v) for v in (c4, c25, cn)):
            continue
        g4.append(c4)
        g25.append(c25)
        gnoise.append(cn)
    # frozen baselines for the same models, on the full cohort
    out[lab] = dict(
        n=n, events=int(e.sum()), reps=len(g4),
        frozen_c4=float(frozen4),
        refit_c4_mean=float(np.mean(g4)), refit_c25_mean=float(np.mean(g25)),
        refit_cnoise_mean=float(np.mean(gnoise)),
        gain_noise=float(np.mean(gnoise) - np.mean(g4)),
        gain_25=float(np.mean(g25) - np.mean(g4)))
    o = out[lab]
    print(f"  {lab}: refit C4 {o['refit_c4_mean']:.3f}  C25 {o['refit_c25_mean']:.3f}  "
          f"C(4+21 noise) {o['refit_cnoise_mean']:.3f}")
    print(f"        gain over refit C4: real 25-gene +{o['gain_25']:.3f}, "
          f"noise 25-param +{o['gain_noise']:.3f}")

print("\n[2] split-half inside each cohort (within-cohort sampling only)")
mr4 = {r["detail"].split(" ")[0]: float(r["value"])
       for r in csv.DictReader(open(S + "/Table_S17_MRG4部署常量.csv", encoding="utf-8-sig"))
       if r["item"] == "系数"}
b4 = np.array([mr4[g] for g in G4])
b25 = np.array([m25[g] for g in G25])
for lab, c in (("CGGA-325", c325), ("CGGA-693", c693)):
    t, e, Z = c
    X4 = np.column_stack([Z[g] for g in G4])
    X25 = np.column_stack([Z[g] for g in G25])
    n = len(t)
    d4, d25, bw = [], [], []
    for _ in range(200):
        perm = rng.permutation(n)
        h1, h2 = perm[: n // 2], perm[n // 2:]
        a4 = concordance_index(t[h1], -(X4[h1] @ b4), e[h1])
        c4 = concordance_index(t[h2], -(X4[h2] @ b4), e[h2])
        a25 = concordance_index(t[h1], -(X25[h1] @ b25), e[h1])
        c25 = concordance_index(t[h2], -(X25[h2] @ b25), e[h2])
        d4.append(abs(a4 - c4))
        d25.append(abs(a25 - c25))
        bw.append((a4 - a25) - (c4 - c25))
    o = dict(reps=200,
             fall_c4_p95=float(np.percentile(d4, 95)),
             fall_c4_max=float(np.max(d4)),
             fall_c25_p95=float(np.percentile(d25, 95)),
             share_ge_0095=float((np.array(d4) >= 0.095).mean()),
             share_ge_0016=float((np.array(d25) >= 0.016).mean()),
             between_model=[float(np.percentile(bw, 2.5)), float(np.percentile(bw, 97.5))])
    out[lab]["splithalf_within"] = o
    print(f"  {lab}: |fall| MRG-4 p95 {o['fall_c4_p95']:.3f} (max {o['fall_c4_max']:.3f}), "
          f"25-gene p95 {o['fall_c25_p95']:.3f}; between-model {o['between_model'][0]:+.3f} to {o['between_model'][1]:+.3f}")
    print(f"        share >= 0.095: {o['share_ge_0095']:.1%};  >= 0.016: {o['share_ge_0016']:.1%}")

json.dump(out, open(OUT + "/analysis_noise_control.json", "w"), indent=2)
print("\nwritten analysis_noise_control.json")
