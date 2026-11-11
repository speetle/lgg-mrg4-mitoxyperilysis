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

"""Round 1 deliverable: Supplementary Table S37 (decision curves, both cohorts,
two calibrations) and the JSON of record."""
import csv, json, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter, KaplanMeierFitter

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = ROOT
EX = os.environ.get("MRG4_CGGA", os.path.join(ROOT, "raw", "cgga"))
S = os.environ.get("MRG4_DERIVED", os.path.join(ROOT, "data"))
OUT = os.environ.get("MRG4_OUT", os.path.join(ROOT, "replication"))
THRESH = np.round(np.arange(0.20, 0.61, 0.05), 2)
HORIZON = 36.0

rows = [r for r in csv.DictReader(open(OUT + "//Table_S34_CGGA693_per_patient.csv",
                                    encoding="utf-8-sig")) if r["Included"] == "yes"]
c693 = dict(t=np.array([float(r["OS_months"]) for r in rows]),
            e=np.array([float(r["Death"]) for r in rows], int),
            s4=np.array([float(r["MRG4_score_raw_scale"]) for r in rows]),
            s25=np.array([float(r["ModelA_score_raw_scale"]) for r in rows]),
            age=np.array([float(r["Age"]) if r["Age"] not in ("", "NA") else np.nan for r in rows]),
            male=np.array([1.0 if r["Sex"].startswith("Male") else 0.0 for r in rows]),
            idhwt=np.array([np.nan if r["IDH_status"] not in ("Mutant", "Wildtype")
                            else (1.0 if r["IDH_status"] == "Wildtype" else 0.0) for r in rows]),
            g3=np.array([1.0 if r["Grade"].strip().endswith("III") else 0.0 for r in rows]))

s5 = [r for r in csv.DictReader(open(S + "/Table_S5_CGGA325_MRG4逐患者评分.csv",
                                    encoding="utf-8-sig"))
      if r["Included_in_primary_validation"] == "True"]
s3 = {r["Patient_ID"]: r for r in csv.DictReader(
    open(S + "/Table_S3_CGGA325逐患者数据.csv", encoding="utf-8-sig"))}
def _g325(pid, col):
    v = s3[pid][col]
    return float(v) if v not in ("", "NA") else np.nan
c325 = dict(t=np.array([float(r["OS_months"]) for r in s5]),
            e=np.array([float(r["Death"]) for r in s5], int),
            s4=np.array([float(r["MRG4_lp_cohort_z"]) for r in s5]),
            s25=np.array([float(s3[r["Patient_ID"]]["ModelA_score_25gene_cohort_z"]) for r in s5]),
            age=np.array([_g325(r["Patient_ID"], "Age") for r in s5]),
            male=np.array([1.0 if s3[r["Patient_ID"]]["Sex"].startswith("Male") else 0.0 for r in s5]),
            idhwt=np.array([np.nan if s3[r["Patient_ID"]]["IDH_status"] not in ("Mutant", "Wildtype")
                            else (1.0 if s3[r["Patient_ID"]]["IDH_status"] == "Wildtype" else 0.0) for r in s5]),
            g3=np.array([1.0 if (s3[r["Patient_ID"]]["WHO_grade"].strip().endswith("III")) else 0.0 for r in s5]))

base = list(csv.DictReader(open(S + "/Table_S19_Breslow基线累积风险.csv",
                                encoding="utf-8-sig")))
bt = np.array([float(r["time_months"]) for r in base])
bs = np.array([float(r["S0"]) for r in base])


def risk_frozen(lp):
    return 1.0 - float(np.interp(HORIZON, bt, bs)) ** np.exp(lp)


def risk_local(c, lp):
    f = CoxPHFitter(penalizer=0.0).fit(pd.DataFrame({"T": c["t"], "E": c["e"], "lp": lp}),
                                       "T", "E")
    H = f.baseline_cumulative_hazard_
    idx = H.index[H.index <= HORIZON]
    return 1.0 - np.exp(-float(H.loc[idx[-1]].iloc[0]) * np.exp(lp))


def curves(c):
    t, e = c["t"], c["e"]
    known = (t > HORIZON) | ((t <= HORIZON) & (e == 1))
    y = np.where((t <= HORIZON) & (e == 1), 1.0, 0.0)
    n = len(t)
    km = KaplanMeierFitter().fit(t, e).survival_function_
    idx = km.index[km.index <= HORIZON]
    s_h = float(km.loc[idx[-1]].iloc[0])

    def nb(risk):
        out = []
        for pt in THRESH:
            m = known & (risk >= pt)
            tp = float(np.sum(y[m] == 1.0))
            fp = float(np.sum(y[m] == 0.0))
            out.append((tp - fp * pt / (1 - pt)) / n)
        return out

    r4f, r25f = risk_frozen(c["s4"]), risk_frozen(c["s25"])
    r4l, r25l = risk_local(c, c["s4"]), risk_local(c, c["s25"])
    # four-variable clinical model refitted inside the cohort (age, sex, IDH, grade)
    d = pd.DataFrame({"T": c["t"], "E": c["e"], "age": c["age"], "male": c["male"],
                      "idhwt": c["idhwt"], "g3": c["g3"]}).dropna().reset_index(drop=True)
    fc = CoxPHFitter(penalizer=0.0).fit(d, "T", "E")
    lp_all = fc.predict_log_partial_hazard(pd.DataFrame(
        {"age": c["age"], "male": c["male"], "idhwt": c["idhwt"], "g3": c["g3"]}).fillna(
        {"age": d["age"].mean(), "idhwt": 0.0})).values
    Hc = fc.baseline_cumulative_hazard_
    idxc = Hc.index[Hc.index <= HORIZON]
    rcl = 1.0 - np.exp(-float(Hc.loc[idxc[-1]].iloc[0]) * np.exp(lp_all))
    # both gene indices refitted inside the cohort, for a like-for-like comparison
    def refit(lp):
        f = CoxPHFitter(penalizer=0.0).fit(pd.DataFrame({"T": c["t"], "E": c["e"], "lp": lp}), "T", "E")
        H = f.baseline_cumulative_hazard_
        i2 = H.index[H.index <= HORIZON]
        return 1.0 - np.exp(-float(H.loc[i2[-1]].iloc[0]) * np.exp(f.params_["lp"] * lp))
    r4r, r25r = refit(c["s4"]), refit(c["s25"])
    res = {
        "frozen": {"MRG-4": nb(r4f), "ModelA": nb(r25f)},
        "recalibrated": {"MRG-4": nb(r4l), "ModelA": nb(r25l)},
        "refitted": {"MRG-4": nb(r4r), "ModelA": nb(r25r), "Clinical": nb(rcl)},
        "treat_all": [(1 - s_h) - s_h * pt / (1 - pt) for pt in THRESH],
        "prevalence36": float(1 - s_h),
        "mean_risk_frozen": [float(np.mean(r4f)), float(np.mean(r25f))],
        "mean_risk_local": [float(np.mean(r4l)), float(np.mean(r25l)), float(np.mean(rcl))],
    }
    return res


out = {"thresholds": [float(x) for x in THRESH], "cohorts": {}}
lines = [["Cohort", "Calibration", "Strategy"] + [f"NB_at_{t:.2f}" for t in THRESH]]
for lab, c in (("CGGA-325", c325), ("CGGA-693", c693)):
    r = curves(c)
    out["cohorts"][lab] = r
    for scale in ("frozen", "recalibrated", "refitted"):
        for strat in r[scale]:
            lines.append([lab, scale, strat] + [round(v, 4) for v in r[scale][strat]])
        lines.append([lab, scale, "treat_all"] + [round(v, 4) for v in r["treat_all"]])

with open(OUT + "/Table_S37_decision_curve.csv", "w", newline="", encoding="utf-8-sig") as fh:
    csv.writer(fh).writerows(lines)
json.dump(out, open(OUT + "/analysis_round1_dca.json", "w"), indent=2)

for lab in out["cohorts"]:
    r = out["cohorts"][lab]
    for scale in ("frozen", "recalibrated", "refitted"):
        a, b = r[scale]["MRG-4"], r[scale]["ModelA"]
        wins = sum(1 for x, y in zip(a, b) if x > y)
        gap = max(abs(x - y) for x, y in zip(a, b))
        print(f"{lab} {scale:<13} MRG-4 higher at {wins}/9 thresholds, max |gap| {gap:.4f}")
        print("   MRG-4 :", " ".join(f"{x:+.4f}" for x in a))
        print("   ModelA:", " ".join(f"{x:+.4f}" for x in b))
    print(f"   treat-all:", " ".join(f"{x:+.4f}" for x in r["treat_all"]))
    print(f"   36-month event rate (KM): {r['prevalence36']:.4f}")
print("written Table_S37_decision_curve.csv + analysis_round1_dca.json")
