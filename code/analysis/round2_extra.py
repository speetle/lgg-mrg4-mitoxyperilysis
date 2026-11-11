#!/usr/bin/env python3
"""RMST, grade-specific hazard ratios and IPCW time-dependent AUC in CGGA-693.

Inputs and environment variables are the same as for the other v1.1.0 control
scripts; see the README. The recorded output is replication/analysis_round2.json.
"""
import os

"""Round-2 (wave 2) analyses: grade-specific hazard ratios, RMST, time-dependent
AUC and the interval for the 0.016 movement in CGGA-693."""
import csv, json, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = ROOT
EX = os.environ.get("MRG4_CGGA", os.path.join(ROOT, "raw", "cgga"))
S = os.environ.get("MRG4_DERIVED", os.path.join(ROOT, "data"))
OUT = os.environ.get("MRG4_OUT", os.path.join(ROOT, "replication"))
rows = [r for r in csv.DictReader(open(OUT + "//Table_S34_CGGA693_per_patient.csv",
                                    encoding="utf-8-sig")) if r["Included"] == "yes"]
t = np.array([float(r["OS_months"]) for r in rows])
e = np.array([float(r["Death"]) for r in rows], int)
s4 = np.array([float(r["MRG4_score_raw_scale"]) for r in rows])
s25 = np.array([float(r["ModelA_score_raw_scale"]) for r in rows])
grade = np.array([r["Grade"].strip().endswith("III") for r in rows])
idh = np.array([r["IDH_status"] for r in rows])
out = {}

# ------------------------------------------------------- 1. HR within grade
print("=== per-SD hazard ratio within WHO grade ===")
out["hr_by_grade"] = {}
for g, lab in ((False, "II"), (True, "III")):
    m = grade == g
    for nm, sc in (("MRG-4", s4), ("ModelA", s25)):
        f = CoxPHFitter(penalizer=0.0).fit(
            pd.DataFrame({"T": t[m], "E": e[m], "s": sc[m]}).reset_index(drop=True), "T", "E")
        sd = float(np.std(sc[m], ddof=1))
        hr = float(np.exp(f.params_["s"] * sd))
        lo, hi = np.exp(f.confidence_intervals_.loc["s"] * sd).tolist()
        out["hr_by_grade"][f"{lab}|{nm}"] = dict(n=int(m.sum()), events=int(e[m].sum()),
                                                 hr_per_sd=hr, ci=[float(lo), float(hi)])
        print(f"  grade {lab} {nm:<7} n={int(m.sum()):3d} ev={int(e[m].sum()):3d} "
              f"HR/SD {hr:.2f} ({lo:.2f}-{hi:.2f})")

# ------------------------------------------------------------- 2. RMST
def rmst(time, event, tau):
    km = KaplanMeierFitter().fit(time, event)
    sf = km.survival_function_
    idx = sf.index[sf.index <= tau]
    vals = sf.loc[idx].iloc[:, 0].values
    xs = np.concatenate([[0.0], idx.values])
    ys = np.concatenate([[1.0], vals])
    return float(np.sum(np.diff(xs) * ys[:-1]))


print("\n=== restricted mean survival, median split, months ===")
rng = np.random.default_rng(20261004)
out["rmst"] = {}
for nm, sc in (("MRG-4", s4), ("ModelA", s25)):
    cut = np.median(sc)
    hi = sc >= cut
    for tau in (36, 60):
        a, b = rmst(t[hi], e[hi], tau), rmst(t[~hi], e[~hi], tau)
        boots = []
        n = len(t)
        for _ in range(1000):
            ii = rng.integers(0, n, n)
            h2 = sc[ii] >= cut
            if e[ii][h2].sum() < 5 or e[ii][~h2].sum() < 5:
                continue
            boots.append(rmst(t[ii][~h2], e[ii][~h2], tau) - rmst(t[ii][h2], e[ii][h2], tau))
        lo, up = np.percentile(boots, [2.5, 97.5])
        out["rmst"][f"{nm}|{tau}"] = dict(low=float(b), high=float(a), diff=float(b - a),
                                          ci=[float(lo), float(up)],
                                          n_low=int((~hi).sum()), n_high=int(hi.sum()),
                                          ev_low=int(e[~hi].sum()), ev_high=int(e[hi].sum()))
        print(f"  {nm:<7} {tau} mo: low-risk {b:5.1f}  high-risk {a:5.1f}  "
              f"difference {b - a:+5.1f} ({lo:+.1f} to {up:+.1f})")

# ------------------------------------------------- 3. time-dependent AUC
def td_auc(time, event, marker, horizons=(24, 36, 60)):
    """Cumulative/dynamic AUC with inverse-probability-of-censoring weights."""
    from sksurv.metrics import cumulative_dynamic_auc
    from sksurv.util import Surv
    y = Surv.from_arrays(event.astype(bool), time)
    auc, mean_auc = cumulative_dynamic_auc(y, y, marker, horizons)
    return auc


print("\n=== time-dependent AUC in CGGA-693 ===")
out["td_auc"] = {}
HZ = (24, 36, 60)
auc4 = td_auc(t, e, s4, HZ)
auc25 = td_auc(t, e, s25, HZ)
n = len(t)
bootstrap = {h: [] for h in HZ}
for _ in range(600):
    ii = rng.integers(0, n, n)
    if e[ii].sum() < 20:
        continue
    try:
        x = td_auc(t[ii], e[ii], s4[ii], HZ)
        y = td_auc(t[ii], e[ii], s25[ii], HZ)
    except Exception:
        continue
    for k, h in enumerate(HZ):
        bootstrap[h].append(x[k] - y[k])
for k, h in enumerate(HZ):
    lo, up = np.percentile(bootstrap[h], [2.5, 97.5])
    out["td_auc"][str(h)] = dict(mrg4=float(auc4[k]), modelA=float(auc25[k]),
                                 diff=float(auc4[k] - auc25[k]), ci=[float(lo), float(up)])
    print(f"  {h} mo: {auc4[k]:.3f} vs {auc25[k]:.3f}  diff {auc4[k] - auc25[k]:+.3f} "
          f"({lo:+.3f} to {up:+.3f})")

# ------------------------------------- 4. interval for the 0.016 movement
print("\n=== interval for the CGGA-325 -> CGGA-693 movement of the 25-gene index ===")
s5 = [r for r in csv.DictReader(open(S + "/Table_S5_CGGA325_MRG4逐患者评分.csv",
                                    encoding="utf-8-sig"))
      if r["Included_in_primary_validation"] == "True"]
s3 = {r["Patient_ID"]: r for r in csv.DictReader(
    open(S + "/Table_S3_CGGA325逐患者数据.csv", encoding="utf-8-sig"))}
t3 = np.array([float(r["OS_months"]) for r in s5])
e3 = np.array([float(r["Death"]) for r in s5], int)
a25 = np.array([float(s3[r["Patient_ID"]]["ModelA_score_25gene_cohort_z"]) for r in s5])
a4 = np.array([float(r["MRG4_lp_cohort_z"]) for r in s5])
d25 = concordance_index(t3, -a25, e3) - concordance_index(t, -s25, e)
d4 = concordance_index(t3, -a4, e3) - concordance_index(t, -s4, e)
boots = []
for _ in range(2000):
    i = rng.integers(0, len(t3), len(t3))
    j = rng.integers(0, len(t), len(t))
    if e3[i].sum() < 10 or e[j].sum() < 10:
        continue
    boots.append((concordance_index(t3[i], -a25[i], e3[i]) - concordance_index(t[j], -s25[j], e[j])))
lo, up = np.percentile(boots, [2.5, 97.5])
p = 2 * min((np.array(boots) <= 0).mean(), (np.array(boots) >= 0).mean())
out["movement"] = dict(modelA=float(d25), ci=[float(lo), float(up)], p=float(p),
                       mrg4=float(d4))
print(f"  25-gene index: {d25:+.4f} ({lo:+.4f} to {up:+.4f}), P = {p:.3f}")
print(f"  MRG-4        : {d4:+.4f}")

json.dump(out, open(OUT + "/analysis_round2.json", "w"), indent=2)
print("\nwritten analysis_round2.json")
