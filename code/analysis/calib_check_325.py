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

"""Two checks before anything goes into the manuscript.

(a) CGGA-325 stores the score in five different scalings. Only one of them can be
    the one the frozen baseline was estimated against; identify it by requiring
    the mean predicted 3-year survival to reproduce the deposited CITL of 1.22.
(b) Continuous NRI has no cut-points; the median-split version is degenerate and
    returned exactly zero. Recompute with the continuous definition.
"""
import csv, json, warnings
import numpy as np
warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = ROOT
EX = os.environ.get("MRG4_CGGA", os.path.join(ROOT, "raw", "cgga"))
S = os.environ.get("MRG4_DERIVED", os.path.join(ROOT, "data"))
OUT = os.environ.get("MRG4_OUT", os.path.join(ROOT, "replication"))
SUP = OUT
rng = np.random.default_rng(20261111)


def rd(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))


base = rd(S + "/Table_S19_Breslow基线累积风险.csv")
bt = np.array([float(r["time_months"]) for r in base])
bs = np.array([float(r["S0"]) for r in base])
o = np.argsort(bt)
bt, bs = bt[o], bs[o]


def S0(t):
    return float(bs[max(np.searchsorted(bt, t, side="right") - 1, 0)])


def km(t, e, tt):
    t, e = np.asarray(t, float), np.asarray(e, int)
    uniq = np.unique(t[e == 1])
    s, var = 1.0, 0.0
    for u in uniq:
        if u > tt:
            break
        at = int((t >= u).sum())
        d = int(((t == u) & (e == 1)).sum())
        if at > d > 0:
            s *= 1 - d / at
            var += d / (at * (at - d))
    se = s * np.sqrt(var) if var > 0 else np.nan
    return s, s - 1.96 * se, s + 1.96 * se


# ---- (a) which 325 scaling belongs to the frozen baseline? -----------------
rows = [r for r in rd(S + "/Table_S5_CGGA325_MRG4逐患者评分.csv")
        if r["Included_in_primary_validation"] == "True"]
t = np.array([float(r["OS_months"]) for r in rows])
e = np.array([int(float(r["Death"])) for r in rows])
obs36, lo36, hi36 = km(t, e, 36.0)
print(f"CGGA-325: n={len(t)} events={int(e.sum())} observed 3-yr {obs36:.3f} ({lo36:.3f}-{hi36:.3f})")
print("target: mean predicted 3-yr = 1.22 x observed = %.3f\n" % (1.22 * obs36))
cands = ["MRG4_lp_cohort_z", "MRG4_lp_discoverySD_scaled", "MRG4_lp_cohort_z_perSD"]
best = None
for c in cands:
    lp = np.array([float(r[c]) for r in rows])
    pred = S0(36.0) ** np.exp(lp)
    r_ = float(np.mean(pred) / obs36)
    print(f"  {c:<32} mean pred {np.mean(pred):.3f}  ratio {r_:.3f}")
    if best is None or abs(r_ - 1.22) < abs(best[1] - 1.22):
        best = (c, r_)
col, ratio = best
print(f"\n  -> {col} (ratio {ratio:.3f}) is the scaling the baseline was built on\n")

lp = np.array([float(r[col]) for r in rows])
out = {}
cal = []
for h in (12, 24, 36, 60):
    p = float(np.mean(S0(h) ** np.exp(lp)))
    ob, l, u = km(t, e, h)
    cal.append(dict(horizon=h, pred=p, obs=float(ob), lo=float(l), hi=float(u),
                    ratio=p / ob))
    print(f"  {h:>3} mo: pred {p:.3f} obs {ob:.3f} ({l:.3f}-{u:.3f}) ratio {p/ob:.3f}")
out["CGGA-325_calibration"] = dict(column=col, rows=cal)

# ---- (b) continuous NRI / IDI in CGGA-693 ---------------------------------
r693 = [r for r in rd(SUP + "/Table_S34_CGGA693_per_patient.csv") if r["Included"] == "yes"]
t6 = np.array([float(r["OS_months"]) for r in r693])
e6 = np.array([int(float(r["Death"])) for r in r693])
lp4 = np.array([float(r["MRG4_score_raw_scale"]) for r in r693])
lp25 = np.array([float(r["ModelA_score_raw_scale"]) for r in r693])

y = ((t6 <= 36.0) & (e6 == 1)).astype(int)
known = (t6 < 36.0) | (y == 1)      # status at 36 months is known
p4 = (1 - S0(36.0) ** np.exp(lp4))[known]
p25 = (1 - S0(36.0) ** np.exp(lp25))[known]
yy = y[known]


def cnri(pa, pb, yy):
    """Continuous NRI: pa is the new model, pb the reference."""
    up = pa > pb
    dn = pa < pb
    return (up[yy == 1].mean() - dn[yy == 1].mean()) + (dn[yy == 0].mean() - up[yy == 0].mean())


def idi(pa, pb, yy):
    return (pa[yy == 1].mean() - pa[yy == 0].mean()) - (pb[yy == 1].mean() - pb[yy == 0].mean())


nri = cnri(p4, p25, yy)
idv = idi(p4, p25, yy)
n = len(yy)
bs_n, bs_i = [], []
for _ in range(2000):
    ii = rng.integers(0, n, n)
    if yy[ii].sum() < 10 or (yy[ii] == 0).sum() < 10:
        continue
    bs_n.append(cnri(p4[ii], p25[ii], yy[ii]))
    bs_i.append(idi(p4[ii], p25[ii], yy[ii]))
print(f"\nCGGA-693 continuous NRI at 36 months: {nri:+.3f} "
      f"({np.percentile(bs_n,2.5):+.3f} to {np.percentile(bs_n,97.5):+.3f})")
print(f"  IDI {idv:+.4f} ({np.percentile(bs_i,2.5):+.4f} to {np.percentile(bs_i,97.5):+.4f}); "
      f"n={n}, events={int(yy.sum())}")
p_nri = 2 * min((np.array(bs_n) <= 0).mean(), (np.array(bs_n) >= 0).mean())
p_idi = 2 * min((np.array(bs_i) <= 0).mean(), (np.array(bs_i) >= 0).mean())
print(f"  P(NRI) = {p_nri:.3f}, P(IDI) = {p_idi:.3f}")
out["CGGA693_nri_idi"] = dict(n=int(n), events=int(yy.sum()), nri=float(nri),
                              nri_ci=[float(np.percentile(bs_n, 2.5)), float(np.percentile(bs_n, 97.5))],
                              nri_p=float(p_nri), idi=float(idv),
                              idi_ci=[float(np.percentile(bs_i, 2.5)), float(np.percentile(bs_i, 97.5))],
                              idi_p=float(p_idi))
json.dump(out, open(OUT + "/analysis_v14b.json", "w"), indent=2)
print("\nwritten analysis_v14b.json")
