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

"""R58 main 1-3 and R60 main 3: calibration at four horizons with interval
estimates, a recalibration interaction test across the two cohorts, and the
NRI/IDI that the replication cohort has so far been missing.

Prediction uses the frozen Breslow baseline: S(t | lp) = S0(t) ^ exp(lp), with
lp the raw score on the scale the baseline was estimated on. The convention is
verified against the deposited calibration table before anything is computed.
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
SUP = OUT
B = 2000
rng = np.random.default_rng(20261111)
HOR = [12, 24, 36, 60]


def rd(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))


# ---- frozen baseline ------------------------------------------------------
base = rd(S + "/Table_S19_Breslow基线累积风险.csv")
bt = np.array([float(r["time_months"]) for r in base])
bs = np.array([float(r["S0"]) for r in base])
order = np.argsort(bt)
bt, bs = bt[order], bs[order]


def S0(t):
    k = np.searchsorted(bt, t, side="right") - 1
    k = max(k, 0)
    return float(bs[k])


# ---- cohorts --------------------------------------------------------------
c693 = rd(SUP + "/Table_S34_CGGA693_per_patient.csv")
c693 = [r for r in c693 if r["Included"] == "yes"]
d693 = dict(t=np.array([float(r["OS_months"]) for r in c693]),
            e=np.array([int(float(r["Death"])) for r in c693]),
            lp=np.array([float(r["MRG4_score_raw_scale"]) for r in c693]),
            lp25=np.array([float(r["ModelA_score_raw_scale"]) for r in c693]),
            grade=[r["Grade"] for r in c693], idh=[r["IDH_status"] for r in c693],
            sex=[r["Sex"] for r in c693],
            age=np.array([float(r["Age"]) if r["Age"] not in ("", "NA") else np.nan for r in c693]))

c325 = rd(S + "/Table_S5_CGGA325_MRG4逐患者评分.csv")
c325 = [r for r in c325 if r["Included_in_primary_validation"] == "True"]
d325 = dict(t=np.array([float(r["OS_months"]) for r in c325]),
            e=np.array([int(float(r["Death"])) for r in c325]),
            lp=np.array([float(r["MRG4_lp_discoverySD_scaled"]) for r in c325]),
            grade=[r["Grade"] for r in c325])
print("693:", len(d693["t"]), int(d693["e"].sum()), " 325:", len(d325["t"]), int(d325["e"].sum()))

# ---- verify the prediction convention on the deposited table --------------
s36 = S0(36.0)
q = np.digitize(d693["lp"], np.percentile(d693["lp"], [33.333, 66.667]))
pred_q = [float(np.mean(s36 ** np.exp(d693["lp"][q == k]))) for k in range(3)]
dep = rd(SUP + "/Table_S35_calibration_CGGA693.csv")
dep_pred = [float(r["Predicted_S36"]) for r in dep]
print("convention check, predicted 3-year by tertile:")
print("  rebuilt:", [round(x, 4) for x in pred_q])
print("  deposited:", dep_pred)
assert all(abs(a - b) < 0.004 for a, b in zip(pred_q, dep_pred)), "convention mismatch"
print("  -> convention confirmed")


# ---- KM ------------------------------------------------------------------
def km(t, e, tt):
    t = np.asarray(t, float)
    e = np.asarray(e, int)
    if len(t) == 0:
        return np.nan, np.nan, np.nan
    o = np.argsort(t)
    t, e = t[o], e[o]
    uniq = np.unique(t[e == 1])
    s, var, n = 1.0, 0.0, len(t)
    at_risk = len(t)
    prev = 0.0
    for u in uniq:
        if u > tt:
            break
        d = int(((t == u) & (e == 1)).sum())
        if at_risk >= d and at_risk > 0:
            s *= 1 - d / at_risk
            var += d / (at_risk * (at_risk - d)) if at_risk > d else 0.0
        at_risk -= int((t < u).sum()) - int((t < prev).sum())
        at_risk = int((t >= u).sum())
        prev = u
    se = s * np.sqrt(var) if var > 0 else np.nan
    return s, s - 1.96 * se, s + 1.96 * se


def calib(d, horizons=HOR):
    out = []
    for h in horizons:
        p = S0(h) ** np.exp(d["lp"])
        obs, lo, hi = km(d["t"], d["e"], h)
        out.append(dict(horizon=h, n=len(p), events=int(d["e"].sum()),
                        pred=float(np.mean(p)), obs=float(obs),
                        lo=float(lo), hi=float(hi),
                        ratio=float(np.mean(p) / obs) if obs == obs else np.nan))
    return out


def slope(d):
    df = pd.DataFrame({"lp": d["lp"], "T": d["t"], "E": d["e"]}).dropna()
    f = CoxPHFitter(penalizer=0.0).fit(df, "T", "E")
    return float(f.params_["lp"]), tuple(float(x) for x in f.confidence_intervals_.loc["lp"])


def ici(d, k=10):
    p = S0(36.0) ** np.exp(d["lp"])
    edges = np.quantile(p, np.linspace(0, 1, k + 1))
    edges[-1] += 1e-9
    idx = np.clip(np.digitize(p, edges) - 1, 0, k - 1)
    tot = 0.0
    for g in range(k):
        m = idx == g
        if m.sum() == 0:
            continue
        pred = float(np.mean(p[m]))
        obs, _, _ = km(d["t"][m], d["e"][m], 36.0)
        if obs != obs:
            continue
        tot += m.sum() * abs(pred - obs)
    return tot / len(p)


def boot(fn, d, reps=400):
    n = len(d["t"])
    vals = []
    for _ in range(reps):
        ii = rng.integers(0, n, n)
        sub = {k: (v[ii] if isinstance(v, np.ndarray) else [v[j] for j in ii])
               for k, v in d.items()}
        try:
            v = fn(sub)
        except Exception:
            v = np.nan
        vals.append(v)
    vals = np.array([v for v in vals if v == v])
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


res = {}
for lab, d in (("CGGA-693", d693), ("CGGA-325", d325)):
    c = calib(d)
    sl, sl_ci = slope(d)
    ic = ici(d)
    ic_ci = boot(lambda x: ici(x), d)
    res[lab] = dict(calibration=c, slope=sl, slope_ci=list(sl_ci),
                    ici=ic, ici_ci=ic_ci)
    print(f"\n{lab}: slope {sl:.3f} ({sl_ci[0]:.3f}-{sl_ci[1]:.3f})  ICI {ic:.4f} ({ic_ci[0]:.4f}-{ic_ci[1]:.4f})")
    for r in c:
        print(f"  {r['horizon']:>3} mo: pred {r['pred']:.3f} obs {r['obs']:.3f} "
              f"({r['lo']:.3f}-{r['hi']:.3f}) ratio {r['ratio']:.3f}")

# ---- recalibration interaction across cohorts -----------------------------
print("\n[recalibration interaction] logit(obs) ~ lp + cohort + lp:cohort, 36 mo")
rows = []
for d, cid in ((d693, 0), (d325, 1)):
    p = np.clip(S0(36.0) ** np.exp(d["lp"]), 1e-6, 1 - 1e-6)
    # patient-level observed status at 36 months; use cloglog link via Cox on the
    # complementary log-log, which is the grouped-survival equivalent
    y = ((d["t"] <= 36.0) & (d["e"] == 1)).astype(int)
    atrisk = d["t"] >= np.minimum(d["t"], 36.0)
    for k in range(len(y)):
        rows.append(dict(lp=d["lp"][k], cohort=cid, y=int(y[k])))
inter = pd.DataFrame(rows)
# Cox model on the binary 36-month outcome is the complementary log-log fit
inter["T"] = np.where(inter["y"] == 1, 1.0, 2.0)
f0 = CoxPHFitter(penalizer=0.0).fit(inter[["lp", "T", "y"]].rename(columns={"y": "E"}), "T", "E")
inter["lpc"] = inter["lp"] * inter["cohort"]
f1 = CoxPHFitter(penalizer=0.0).fit(
    inter[["lp", "cohort", "lpc", "T", "y"]].rename(columns={"y": "E"}), "T", "E")
lr = 2 * (f1.log_likelihood_ - f0.log_likelihood_)
from scipy.stats import chi2
p_int = float(chi2.sf(lr, 2))
res["recalibration_interaction"] = dict(lr=float(lr), df=2, p=p_int,
                                        coef_cohort=float(f1.params_["cohort"]),
                                        coef_interaction=float(f1.params_["lpc"]))
print(f"  LR = {lr:.3f} on 2 df, P = {p_int:.4f}; cohort {f1.params_['cohort']:.3f}, "
      f"lp x cohort {f1.params_['lpc']:.3f}")

# ---- NRI / IDI at 36 months in CGGA-693 ----------------------------------
print("\n[NRI/IDI, CGGA-693, 36 months]")
d = d693
y = ((d["t"] <= 36.0) & (d["e"] == 1)).astype(int)
cens = d["t"] < 36.0
elig = ~cens | (y == 1)          # known status at 36 months
m = elig
p4 = np.clip(S0(36.0) ** np.exp(d["lp"]), 1e-6, 1 - 1e-6)[m]
p25 = np.clip(S0(36.0) ** np.exp(d["lp25"]), 1e-6, 1 - 1e-6)[m]
yy = y[m]
r4, r25 = 1 - p4, 1 - p25        # predicted risk
up4 = r4[yy == 1] >= np.median(r4)
up25 = r25[yy == 1] >= np.median(r25)
dn4 = r4[yy == 0] >= np.median(r4)
dn25 = r25[yy == 0] >= np.median(r25)
nri_e = up4.mean() - up25.mean()
nri_ne = dn25.mean() - dn4.mean()
nri = nri_e + nri_ne
idi = (r4[yy == 1].mean() - r4[yy == 0].mean()) - (r25[yy == 1].mean() - r25[yy == 0].mean())
res["nri_idi_693"] = dict(n=float(m.sum()), events=int(yy.sum()),
                          nri_event=float(nri_e), nri_nonevent=float(nri_ne),
                          nri=float(nri), idi=float(idi))
print(f"  eligible n = {int(m.sum())}, events = {int(yy.sum())}")
print(f"  NRI {nri:+.3f} (event {nri_e:+.3f}, non-event {nri_ne:+.3f}), IDI {idi:+.4f}")

json.dump(res, open(OUT + "/analysis_calibration_v14.json", "w"), indent=2)
print("\nwritten analysis_calibration_v14.json")
