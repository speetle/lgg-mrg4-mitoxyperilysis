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

"""Bootstrap that carries coefficient uncertainty: refit both indices inside every
resample, then evaluate the fitted coefficients on the whole cohort."""
import csv, json, os, warnings
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
B = 500
rng = np.random.default_rng(20261004)

m25 = {r["Gene"]: float(r["Multivariate_coefficient_ridge_L2"])
       for r in csv.DictReader(open(S + "/Table_S12_ModelA_25基因系数.csv", encoding="utf-8-sig"))}


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


def fit_eval(Z, genes, t, e, idx_fit, idx_eval):
    """Fit Cox on the resample, evaluate C on idx_eval."""
    X = np.column_stack([Z[g][idx_fit] for g in genes])
    d = pd.DataFrame(X, columns=[f"g{i}" for i in range(len(genes))])
    d["T"] = t[idx_fit]
    d["E"] = e[idx_fit]
    ok = d.notna().all(axis=1)
    d = d[ok].reset_index(drop=True)
    if len(d) < 50 or d["E"].sum() < 20:
        return np.nan
    try:
        f = CoxPHFitter(penalizer=0.0).fit(d, "T", "E")
    except Exception:
        return np.nan
    beta = f.params_.values
    Xe = np.column_stack([Z[g][idx_eval] for g in genes])
    lp = Xe @ beta
    return concordance_index(t[idx_eval], -lp, e[idx_eval])


c325 = cohort(EX + "/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",
              EX + "/CGGA.mRNAseq_325_clinical.20200506.txt", G4 + list(m25))
c693 = cohort(EX + "/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",
              EX + "/CGGA.mRNAseq_693_clinical.20200506.txt", G4 + list(m25))
print("325", len(c325[0]), int(c325[1].sum()), " 693", len(c693[0]), int(c693[1].sum()))

res = {}
for lab, c in (("CGGA-325", c325), ("CGGA-693", c693)):
    t, e, Z = c
    n = len(t)
    all_idx = np.arange(n)
    d4, d25, both = [], [], []
    for b in range(B):
        ii = rng.integers(0, n, n)
        c4 = fit_eval(Z, G4, t, e, ii, all_idx)
        c25 = fit_eval(Z, list(m25), t, e, ii, all_idx)
        if np.isnan(c4) or np.isnan(c25):
            continue
        d4.append(c4)
        d25.append(c25)
        both.append(c4 - c25)
    res[lab] = dict(n=n, events=int(e.sum()), reps=len(both),
                    c4=[float(np.percentile(d4, 2.5)), float(np.percentile(d4, 97.5))],
                    c25=[float(np.percentile(d25, 2.5)), float(np.percentile(d25, 97.5))],
                    diff=[float(np.percentile(both, 2.5)), float(np.percentile(both, 97.5))],
                    diff_point=float(np.mean(both)))
    print(f"{lab}: refit-evaluated C4 {res[lab]['c4'][0]:.3f}-{res[lab]['c4'][1]:.3f}, "
          f"C25 {res[lab]['c25'][0]:.3f}-{res[lab]['c25'][1]:.3f}, "
          f"diff {res[lab]['diff'][0]:+.3f} to {res[lab]['diff'][1]:+.3f}")

# interaction with coefficient uncertainty in both cohorts
inter = []
for b in range(B):
    i = rng.integers(0, len(c325[0]), len(c325[0]))
    j = rng.integers(0, len(c693[0]), len(c693[0]))
    a4 = fit_eval(c325[2], G4, c325[0], c325[1], i, np.arange(len(c325[0])))
    a25 = fit_eval(c325[2], list(m25), c325[0], c325[1], i, np.arange(len(c325[0])))
    b4 = fit_eval(c693[2], G4, c693[0], c693[1], j, np.arange(len(c693[0])))
    b25 = fit_eval(c693[2], list(m25), c693[0], c693[1], j, np.arange(len(c693[0])))
    if any(np.isnan(v) for v in (a4, a25, b4, b25)):
        continue
    inter.append((a4 - a25) - (b4 - b25))
res["interaction"] = dict(reps=len(inter),
                          ci=[float(np.percentile(inter, 2.5)), float(np.percentile(inter, 97.5))],
                          p=float(2 * min((np.array(inter) <= 0).mean(),
                                          (np.array(inter) >= 0).mean())))
print(f"interaction with refitting: {res['interaction']['ci'][0]:+.3f} to "
      f"{res['interaction']['ci'][1]:+.3f}, P = {res['interaction']['p']:.3f}")
json.dump(res, open(OUT + "/analysis_refit_bootstrap.json", "w"), indent=2)
print("written analysis_refit_bootstrap.json")
