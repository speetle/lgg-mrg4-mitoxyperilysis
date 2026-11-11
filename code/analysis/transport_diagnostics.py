#!/usr/bin/env python3
"""Transport diagnostics for the frozen MRG-4 score in two CGGA cohorts.

Reads only deposited per-patient tables, so it can be run without the primary
expression matrices:

    data/Table_S34_CGGA693_per_patient.csv   replication cohort (from replication_cgga693.py)
    data/Table_S5_CGGA325_MRG4逐患者评分.csv    discovery-scale validation cohort, MRG-4
    data/Table_S3_CGGA325逐患者数据.csv         CGGA-325, 25-gene Model A score
    data/Table_S19_Breslow基线累积风险.csv      frozen Breslow baseline (TCGA-LGG)

It reports the five quantities the article uses to separate discrimination from
calibration, for both cohorts:

  1. calibration slope     Cox coefficient of the frozen linear predictor refitted
                           alone; 1.0 means the frozen scale is reproduced.
  2. CITL                  ratio of mean predicted to observed 36-month survival
                           (>1 over-predicts survival); patient bootstrap CI.
  3. ICI                   mean absolute predicted-minus-observed difference over
                           risk tertiles, 36 months.
  4. proportional hazards  Schoenfeld residual test on the score.
  5. discrimination        C-index within WHO grade II and grade III.

It also tests whether the *difference* between the two indices is the same in the
two cohorts (difference of differences, patient bootstrap).

Usage:

    python transport_diagnostics.py --data-dir ../data --out transport.json

Nothing here writes into ``data/``; the JSON is a convenience record only.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import warnings

import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import concordance_index
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")

BOOT = 2000
SEED = 20261004


def read(name: str, data_dir: str) -> list[dict]:
    with open(os.path.join(data_dir, name), encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def survival_at(km, t: float) -> float:
    idx = km.index[km.index <= t]
    return float(km.loc[idx[-1]].iloc[0]) if len(idx) else 1.0


def load_cohorts(data_dir: str) -> dict[str, dict]:
    s34 = [r for r in read("Table_S34_CGGA693_per_patient.csv", data_dir) if r["Included"] == "yes"]
    s5 = [r for r in read("Table_S5_CGGA325_MRG4逐患者评分.csv", data_dir)
          if r["Included_in_primary_validation"] == "True"]
    s3 = {r["Patient_ID"]: r for r in read("Table_S3_CGGA325逐患者数据.csv", data_dir)}
    c693 = {
        "label": "CGGA-693",
        "time": np.array([float(r["OS_months"]) for r in s34]),
        "event": np.array([float(r["Death"]) for r in s34], dtype=int),
        "mrg4": np.array([float(r["MRG4_score_raw_scale"]) for r in s34]),
        "modelA": np.array([float(r["ModelA_score_raw_scale"]) for r in s34]),
        "grade": np.array([r["Grade"].strip() for r in s34]),
    }
    c325 = {
        "label": "CGGA-325",
        "time": np.array([float(r["OS_months"]) for r in s5]),
        "event": np.array([float(r["Death"]) for r in s5], dtype=int),
        "mrg4": np.array([float(r["MRG4_lp_cohort_z"]) for r in s5]),
        "modelA": np.array([float(s3[r["Patient_ID"]]["ModelA_score_25gene_cohort_z"]) for r in s5]),
        "grade": np.array([r["Grade"].strip() for r in s5]),
    }
    return {"CGGA-325": c325, "CGGA-693": c693}


def slope(time, event, lp) -> tuple[float, float, float]:
    """Cox coefficient of the frozen linear predictor, with its coefficient-scale CI."""
    import pandas as pd
    fit = CoxPHFitter(penalizer=0.0).fit(
        pd.DataFrame({"T": time, "E": event, "lp": lp}), "T", "E")
    lo, hi = fit.confidence_intervals_.loc["lp"].tolist()
    return float(fit.params_["lp"]), float(lo), float(hi)


def calibration(time, event, lp, baseline_t, baseline_s, horizon: float = 36.0):
    km = KaplanMeierFitter().fit(time, event).survival_function_
    s0 = float(np.interp(horizon, baseline_t, baseline_s))
    pred = float(np.mean(s0 ** np.exp(lp)))
    obs = survival_at(km, horizon)
    return pred, obs, pred / obs


def ici(time, event, lp, baseline_t, baseline_s, horizon: float = 36.0) -> float:
    cuts = np.quantile(lp, [1 / 3, 2 / 3])
    grp = np.digitize(lp, cuts)
    s0 = float(np.interp(horizon, baseline_t, baseline_s))
    diffs = []
    for k in range(3):
        m = grp == k
        if m.sum() == 0:
            continue
        km = KaplanMeierFitter().fit(time[m], event[m]).survival_function_
        diffs.append(abs(float(np.mean(s0 ** np.exp(lp[m]))) - survival_at(km, horizon)))
    return float(np.mean(diffs))


def main(argv: list[str] | None = None) -> int:
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data-dir", default=os.path.join(here, "data"))
    ap.add_argument("--out", default=None, help="optional path for a JSON record")
    args = ap.parse_args(argv)

    baseline = read("Table_S19_Breslow基线累积风险.csv", args.data_dir)
    bt = np.array([float(r["time_months"]) for r in baseline])
    bs = np.array([float(r["S0"]) for r in baseline])

    cohorts = load_cohorts(args.data_dir)
    rng = np.random.default_rng(SEED)
    out: dict[str, object] = {}

    for label, c in cohorts.items():
        t, e, lp = c["time"], c["event"], c["mrg4"]
        sl, exp_lo, exp_hi = slope(t, e, lp)
        pred, obs, citl = calibration(t, e, lp, bt, bs)
        rec = {
            "n": int(len(t)),
            "events": int(e.sum()),
            "c_mrg4": float(concordance_index(t, -lp, e)),
            "c_modelA": float(concordance_index(t, -c["modelA"], e)),
            "calibration_slope": sl,
            "slope_ci": [exp_lo, exp_hi],
            "pred_s36": pred,
            "obs_s36": obs,
            "citl": citl,
            "ici36": ici(t, e, lp, bt, bs),
        }
        import pandas as pd
        ph = proportional_hazard_test(CoxPHFitter().fit(pd.DataFrame({"T": t, "E": e, "lp": lp}), "T", "E"),
                                      pd.DataFrame({"T": t, "E": e, "lp": lp}),
                                      time_transform="rank")
        rec["schoenfeld_rank_p"] = float(ph.summary["p"].iloc[0])

        # bootstrap CI for CITL
        vals = []
        n = len(t)
        for _ in range(BOOT):
            idx = rng.integers(0, n, n)
            if e[idx].sum() < 10:
                continue
            try:
                vals.append(calibration(t[idx], e[idx], lp[idx], bt, bs)[2])
            except Exception:
                continue
        rec["citl_ci"] = [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]

        strata = {}
        for g in ("WHO II", "WHO III"):
            m = c["grade"] == g
            if m.sum() < 10:
                continue
            strata[g] = {
                "n": int(m.sum()),
                "events": int(e[m].sum()),
                "c_mrg4": float(concordance_index(t[m], -lp[m], e[m])),
                "c_modelA": float(concordance_index(t[m], -c["modelA"][m], e[m])),
                "difference": float(concordance_index(t[m], -lp[m], e[m])
                                    - concordance_index(t[m], -c["modelA"][m], e[m])),
            }
        rec["by_grade"] = strata
        out[label] = rec

        print(f"{label}: n={rec['n']}, {rec['events']} deaths")
        print(f"  C-index        MRG-4 {rec['c_mrg4']:.4f}  Model A {rec['c_modelA']:.4f}"
              f"  difference {rec['c_mrg4'] - rec['c_modelA']:+.4f}")
        print(f"  slope          {sl:.3f} ({exp_lo:.3f}-{exp_hi:.3f})")
        print(f"  CITL           {citl:.3f} ({rec['citl_ci'][0]:.3f}-{rec['citl_ci'][1]:.3f}); "
              f"predicted {pred:.3f} vs observed {obs:.3f} at 36 months")
        print(f"  ICI(36)        {rec['ici36']:.3f}")
        print(f"  Schoenfeld     P = {rec['schoenfeld_rank_p']:.3f}")
        for g, s in strata.items():
            print(f"  {g:<9}      C {s['c_mrg4']:.3f} vs {s['c_modelA']:.3f}"
                  f"  (difference {s['difference']:+.3f}; n={s['n']}, {s['events']} deaths)")

    # ------------------------------------------------ difference of differences
    d325 = out["CGGA-325"]["c_mrg4"] - out["CGGA-325"]["c_modelA"]
    d693 = out["CGGA-693"]["c_mrg4"] - out["CGGA-693"]["c_modelA"]
    c325, c693 = cohorts["CGGA-325"], cohorts["CGGA-693"]
    rng = np.random.default_rng(SEED)
    diffs = []
    for _ in range(BOOT):
        i = rng.integers(0, len(c325["time"]), len(c325["time"]))
        j = rng.integers(0, len(c693["time"]), len(c693["time"]))
        if c325["event"][i].sum() < 10 or c693["event"][j].sum() < 10:
            continue
        a = (concordance_index(c325["time"][i], -c325["mrg4"][i], c325["event"][i])
             - concordance_index(c325["time"][i], -c325["modelA"][i], c325["event"][i]))
        b = (concordance_index(c693["time"][j], -c693["mrg4"][j], c693["event"][j])
             - concordance_index(c693["time"][j], -c693["modelA"][j], c693["event"][j]))
        diffs.append(a - b)
    diff_of_diff = float(d325 - d693)
    lo, hi = float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))
    p = float(2 * min((np.array(diffs) <= 0).mean(), (np.array(diffs) >= 0).mean()))
    out["difference_of_differences"] = {
        "difference_325": float(d325), "difference_693": float(d693),
        "diff_of_diff": diff_of_diff, "ci": [lo, hi], "p": p,
    }
    print(f"difference of differences: {diff_of_diff:+.4f} ({lo:+.3f} to {hi:+.3f}), P = {p:.3f}")

    if args.out:
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=2, ensure_ascii=False)
        print("record written to " + args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
