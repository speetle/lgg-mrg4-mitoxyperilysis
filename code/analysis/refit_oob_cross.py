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

"""Two additions to the refit bootstrap, both required to interpret it.

(a) Out-of-bag evaluation: fit on the resample, score only the patients left out.
    This removes the optimism that a 25-parameter fit enjoys over a 4-parameter fit
    when both are scored on the patients they were fitted on.
(b) Cross-cohort refit: fit in one validation cohort, score the other. This is a
    transport estimate in which nothing but the coefficients travels -- the
    z-reference is cohort-internal in both directions, so the two cohorts are on
    the same footing.
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
B = 500
rng = np.random.default_rng(20261004)

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


def fit(Z, genes, t, e, idx):
    X = np.column_stack([Z[g][idx] for g in genes])
    d = pd.DataFrame(X, columns=[f"g{i}" for i in range(len(genes))])
    d["T"] = t[idx]
    d["E"] = e[idx]
    d = d[d.notna().all(axis=1)].reset_index(drop=True)
    if len(d) < 60 or d["E"].sum() < 18:
        return None
    try:
        f = CoxPHFitter(penalizer=0.0).fit(d, "T", "E")
    except Exception:
        return None
    return f.params_.values


def score(beta, Z, genes, t, e, idx):
    Xe = np.column_stack([Z[g][idx] for g in genes])
    lp = Xe @ beta
    ok = np.isfinite(lp)
    if ok.sum() < 40 or e[idx][ok].sum() < 12:
        return np.nan
    return concordance_index(t[idx][ok], -lp[ok], e[idx][ok])


c325 = cohort(EX + "/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",
              EX + "/CGGA.mRNAseq_325_clinical.20200506.txt", G4 + G25)
c693 = cohort(EX + "/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",
              EX + "/CGGA.mRNAseq_693_clinical.20200506.txt", G4 + G25)
print("325", len(c325[0]), int(c325[1].sum()), " 693", len(c693[0]), int(c693[1].sum()))

out = {}

# ---------- (a) out-of-bag -------------------------------------------------
print("\n[a] out-of-bag: fit on resample, score the patients left out")
for lab, c in (("CGGA-325", c325), ("CGGA-693", c693)):
    t, e, Z = c
    n = len(t)
    d4, d25, both = [], [], []
    for b in range(B):
        ii = rng.integers(0, n, n)
        oob = np.setdiff1d(np.arange(n), np.unique(ii))
        if len(oob) < 45:
            continue
        b4 = fit(Z, G4, t, e, ii)
        b25 = fit(Z, G25, t, e, ii)
        if b4 is None or b25 is None:
            continue
        c4 = score(b4, Z, G4, t, e, oob)
        c25 = score(b25, Z, G25, t, e, oob)
        if np.isnan(c4) or np.isnan(c25):
            continue
        d4.append(c4)
        d25.append(c25)
        both.append(c4 - c25)
    out.setdefault(lab, {})["oob"] = dict(
        reps=len(both),
        c4=[float(np.percentile(d4, 2.5)), float(np.percentile(d4, 97.5))],
        c25=[float(np.percentile(d25, 2.5)), float(np.percentile(d25, 97.5))],
        diff=[float(np.percentile(both, 2.5)), float(np.percentile(both, 97.5))],
        diff_point=float(np.mean(both)),
        p=float(2 * min((np.array(both) <= 0).mean(), (np.array(both) >= 0).mean())))
    o = out[lab]["oob"]
    print(f"  {lab}: C4 {o['c4'][0]:.3f}-{o['c4'][1]:.3f}  C25 {o['c25'][0]:.3f}-{o['c25'][1]:.3f}"
          f"  diff {o['diff'][0]:+.3f} to {o['diff'][1]:+.3f}  P={o['p']:.3f}")

# ---------- (b) cross-cohort refit ----------------------------------------
print("\n[b] cross-cohort: fit in one validation cohort, score the other")
pairs = (("CGGA-325", c693, c325), ("CGGA-693", c325, c693))
for lab, ct, cv in pairs:
    tf, ef, Zf = ct
    tv, ev, Zv = cv
    nf = len(tf)
    d4, d25, both = [], [], []
    for b in range(B):
        ii = rng.integers(0, nf, nf)
        b4 = fit(Zf, G4, tf, ef, ii)
        b25 = fit(Zf, G25, tf, ef, ii)
        if b4 is None or b25 is None:
            continue
        c4 = score(b4, Zv, G4, tv, ev, np.arange(len(tv)))
        c25 = score(b25, Zv, G25, tv, ev, np.arange(len(tv)))
        if np.isnan(c4) or np.isnan(c25):
            continue
        d4.append(c4)
        d25.append(c25)
        both.append(c4 - c25)
    lab2 = "fit693_eval325" if lab == "CGGA-325" else "fit325_eval693"
    out.setdefault(lab, {})["cross"] = dict(
        reps=len(both), direction=lab2,
        c4=[float(np.percentile(d4, 2.5)), float(np.percentile(d4, 97.5))],
        c25=[float(np.percentile(d25, 2.5)), float(np.percentile(d25, 97.5))],
        diff=[float(np.percentile(both, 2.5)), float(np.percentile(both, 97.5))],
        diff_point=float(np.mean(both)),
        p=float(2 * min((np.array(both) <= 0).mean(), (np.array(both) >= 0).mean())))
    o = out[lab]["cross"]
    print(f"  {lab2}: C4 {o['c4'][0]:.3f}-{o['c4'][1]:.3f}  C25 {o['c25'][0]:.3f}-{o['c25'][1]:.3f}"
          f"  diff {o['diff'][0]:+.3f} to {o['diff'][1]:+.3f}  P={o['p']:.3f}")

json.dump(out, open(OUT + "/analysis_refit_oob_cross.json", "w"), indent=2)
print("\nwritten analysis_refit_oob_cross.json")
