#!/usr/bin/env python3
"""Replication in a second, non-overlapping CGGA cohort (CGGA-693, WHO grade II-III).

Rebuilds Supplementary Tables S34-S36 from the public CGGA release plus the frozen
objects deposited in ``data/``. Nothing is refitted: the coefficients of MRG-4
(``Table_S17``), the 25-gene ridge index Model A (``Table_S12``) and the Breslow
baseline cumulative hazard (``Table_S19``) are used exactly as frozen at derivation.

Inputs (download from http://www.cgga.org.cn, release 2020-05-06):

    CGGA.mRNAseq_693.RSEM-genes.20200506.txt
    CGGA.mRNAseq_693_clinical.20200506.txt
    CGGA.mRNAseq_325.RSEM-genes.20200506.txt
    CGGA.mRNAseq_325_clinical.20200506.txt

Usage:

    python replication_cgga693.py --cgga-dir /path/to/cgga --out-dir ./replication

Outputs, written to ``--out-dir``:

    Table_S34_CGGA693_per_patient.csv        every profiled patient, with inclusion flag
    Table_S35_calibration_CGGA693.csv        frozen-model calibration by risk tertile
    Table_S36_scale_sensitivity.csv          raw RSEM versus log2(RSEM+1), both cohorts
    replication_summary.json                 headline numbers quoted in the article

Conventions (identical to the discovery and CGGA-325 analyses):

* expression is used on the **raw RSEM** scale in the primary analysis, and on
  log2(RSEM+1) only in the S36 sensitivity analysis;
* each gene is *z*-scored within the cohort being scored (all profiled samples of
  that cohort are the reference set), so the reported C-indices are cohort-*
  z* estimates rather than frozen-scale estimates;
* the risk score is the frozen weighted sum of those *z*-scores;
* survival times are converted from days to months with 30.4375 days per month.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index

GENES4 = ["BMP1", "KIF15", "TRAF3", "PDGFA"]
DAY = 30.4375
BOOT = 2000
SEED = 20261111


# --------------------------------------------------------------------------- io
def read_genes(path: str, genes: list[str]) -> tuple[list[str], dict[str, np.ndarray]]:
    """Read selected rows of a CGGA RSEM matrix (genes x samples, tab separated)."""
    want = set(genes)
    out: dict[str, np.ndarray] = {}
    with open(path, encoding="utf-8", errors="ignore") as fh:
        ids = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in want:
                out[parts[0]] = np.array(
                    [float(x) if x not in ("", "NA") else np.nan for x in parts[1:]]
                )
    return ids, out


def read_clinical(path: str) -> dict[str, dict]:
    with open(path, encoding="utf-8-sig", errors="ignore") as fh:
        return {r["CGGA_ID"]: r for r in csv.DictReader(fh, delimiter="\t")}


def read_coefficients(data_dir: str) -> tuple[dict[str, float], dict[str, float]]:
    coef4 = {
        r["detail"].split()[0]: float(r["value"])
        for r in csv.DictReader(
            open(os.path.join(data_dir, "Table_S17_MRG4部署常量.csv"), encoding="utf-8-sig")
        )
        if r["item"] == "系数"
    }
    coef25 = {
        r["Gene"]: float(r["Multivariate_coefficient_ridge_L2"])
        for r in csv.DictReader(
            open(os.path.join(data_dir, "Table_S12_ModelA_25基因系数.csv"), encoding="utf-8-sig")
        )
    }
    missing = [g for g in GENES4 if g not in coef4]
    if missing:
        raise SystemExit(f"MRG-4 coefficients missing from Table_S17: {missing}")
    return coef4, coef25


# ---------------------------------------------------------------------- scoring
def score(expr: dict[str, np.ndarray], coefs: dict[str, float], log2: bool = False) -> np.ndarray:
    """Weighted sum of cohort-*z*-scored genes over every profiled sample."""
    total = np.zeros(len(next(iter(expr.values()))))
    for gene, weight in coefs.items():
        v = expr[gene].astype(float)
        if log2:
            v = np.log2(v + 1.0)
        total += weight * (v - np.nanmean(v)) / np.nanstd(v, ddof=1)
    return total


def exclusion_reason(clin: dict[str, dict], pid: str) -> str:
    row = clin.get(pid)
    if row is None:
        return "no clinical record"
    if row["Grade"] not in ("WHO II", "WHO III"):
        return "grade IV or ungraded (not in the LGG analysis set)"
    if row["OS"] in ("", "NA"):
        return "overall survival missing"
    if row["Censor (alive=0; dead=1)"] not in ("0", "1"):
        return "censoring status missing"
    return "included"


def tertiles(score_values: np.ndarray) -> np.ndarray:
    cuts = np.quantile(score_values, [1 / 3, 2 / 3])
    return np.digitize(score_values, cuts)


# ------------------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cgga-dir", required=True, help="directory holding the four CGGA 2020-05-06 files")
    ap.add_argument("--out-dir", default=os.path.join(here, "replication"))
    ap.add_argument("--data-dir", default=os.path.join(here, "data"),
                    help="deposited derived tables (default: <repo>/data)")
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    coef4, coef25 = read_coefficients(args.data_dir)
    rsem693 = os.path.join(args.cgga_dir, "CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
    clin693_path = os.path.join(args.cgga_dir, "CGGA.mRNAseq_693_clinical.20200506.txt")
    for path in (rsem693, clin693_path):
        if not os.path.exists(path):
            raise SystemExit(f"missing input: {path}")

    ids, expr = read_genes(rsem693, GENES4 + list(coef25))
    clin = read_clinical(clin693_path)
    absent = [g for g in GENES4 + list(coef25) if g not in expr]
    if absent:
        raise SystemExit(f"genes absent from the 693 matrix: {absent}")

    keep = [i for i, pid in enumerate(ids) if exclusion_reason(clin, pid) == "included"]
    reasons = [exclusion_reason(clin, pid) for pid in ids]

    s4_all = score(expr, {g: coef4[g] for g in GENES4})
    s25_all = score(expr, coef25)
    s4, s25 = s4_all[keep], s25_all[keep]

    time = np.array([float(clin[ids[i]]["OS"]) / DAY for i in keep])
    event = np.array([float(clin[ids[i]]["Censor (alive=0; dead=1)"]) for i in keep], dtype=int)
    grp = tertiles(s4)

    c4 = concordance_index(time, -s4, event)
    c25 = concordance_index(time, -s25, event)

    rng = np.random.default_rng(SEED)
    n = len(keep)
    b4, b25, bd = [], [], []
    for _ in range(BOOT):
        idx = rng.integers(0, n, n)
        if event[idx].sum() < 10:
            continue
        x = concordance_index(time[idx], -s4[idx], event[idx])
        y = concordance_index(time[idx], -s25[idx], event[idx])
        b4.append(x)
        b25.append(y)
        bd.append(x - y)
    ci = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]

    # ---------------------------------------------------------------- S34
    def z_of(gene: str, i: int) -> float:
        v = expr[gene]
        return float((v[i] - np.nanmean(v)) / np.nanstd(v, ddof=1))

    keep_pos = {i: k for k, i in enumerate(keep)}
    s34 = os.path.join(args.out_dir, "Table_S34_CGGA693_per_patient.csv")
    with open(s34, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Patient_ID", "Included", "Exclusion_reason", "BMP1_z", "KIF15_z", "TRAF3_z",
                    "PDGFA_z", "MRG4_score_raw_scale", "ModelA_score_raw_scale", "OS_months",
                    "Death", "Grade", "IDH_status", "Sex", "Age", "Risk_tertile"])
        for i, pid in enumerate(ids):
            row = clin.get(pid, {})
            os_months = (
                "" if not row or row["OS"] in ("", "NA")
                else round(float(row["OS"]) / DAY, 4)
            )
            w.writerow([
                pid, "yes" if i in keep_pos else "no", reasons[i],
                round(z_of("BMP1", i), 4), round(z_of("KIF15", i), 4),
                round(z_of("TRAF3", i), 4), round(z_of("PDGFA", i), 4),
                round(float(s4_all[i]), 6), round(float(s25_all[i]), 6),
                os_months,
                "" if not row else row["Censor (alive=0; dead=1)"],
                "" if not row else row["Grade"],
                "" if not row else row["IDH_mutation_status"],
                "" if not row else row["Gender"].strip(),
                "" if not row else row["Age"],
                "" if i not in keep_pos else int(grp[keep_pos[i]] + 1),
            ])

    # ---------------------------------------------------------------- S35
    baseline = list(csv.DictReader(
        open(os.path.join(args.data_dir, "Table_S19_Breslow基线累积风险.csv"), encoding="utf-8-sig")))
    bt = np.array([float(r["time_months"]) for r in baseline])
    bs = np.array([float(r["S0"]) for r in baseline])

    def s0(t: float) -> float:
        return float(np.interp(t, bt, bs))

    def km_at(curve, t: float) -> float:
        idx = curve.index[curve.index <= t]
        return float(curve.loc[idx[-1]].iloc[0]) if len(idx) else 1.0

    def ci_at(ci_frame, t: float) -> tuple[float, float]:
        idx = ci_frame.index[ci_frame.index <= t]
        if not len(idx):
            return (1.0, 1.0)
        return (float(ci_frame.loc[idx[-1]].iloc[0]), float(ci_frame.loc[idx[-1]].iloc[1]))

    s35 = os.path.join(args.out_dir, "Table_S35_calibration_CGGA693.csv")
    with open(s35, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Tertile", "n", "events", "Predicted_S36", "Observed_S36",
                    "KM_lower95_S36", "KM_upper95_S36", "Predicted_S60", "Observed_S60",
                    "KM_lower95_S60", "KM_upper95_S60"])
        for k in range(3):
            m = grp == k
            kmf = KaplanMeierFitter().fit(time[m], event[m])
            pred36 = float(np.mean(s0(36) ** np.exp(s4[m])))
            pred60 = float(np.mean(s0(60) ** np.exp(s4[m])))
            lo36, hi36 = ci_at(kmf.confidence_interval_, 36)
            lo60, hi60 = ci_at(kmf.confidence_interval_, 60)
            w.writerow([f"Q{k + 1}", int(m.sum()), int(event[m].sum()),
                        round(pred36, 4), round(km_at(kmf.survival_function_, 36), 4),
                        round(lo36, 4), round(hi36, 4),
                        round(pred60, 4), round(km_at(kmf.survival_function_, 60), 4),
                        round(lo60, 4), round(hi60, 4)])
        kmf = KaplanMeierFitter().fit(time, event)
        w.writerow(["all", len(keep), int(event.sum()),
                    round(float(np.mean(s0(36) ** np.exp(s4))), 4),
                    round(km_at(kmf.survival_function_, 36), 4), "", "",
                    round(float(np.mean(s0(60) ** np.exp(s4))), 4),
                    round(km_at(kmf.survival_function_, 60), 4), "", ""])

    # ---------------------------------------------------------------- S36
    s36 = os.path.join(args.out_dir, "Table_S36_scale_sensitivity.csv")
    rows = []
    cohorts = [
        ("CGGA-325", "CGGA.mRNAseq_325.RSEM-genes.20200506.txt", "Table_S5_CGGA325_MRG4逐患者评分.csv"),
        ("CGGA-693", "CGGA.mRNAseq_693.RSEM-genes.20200506.txt", None),
    ]
    for label, matrix, per_patient in cohorts:
        i2, e2 = read_genes(os.path.join(args.cgga_dir, matrix), GENES4 + list(coef25))
        if per_patient is None:
            kk, tt, ee = keep, time, event
        else:
            s5 = {r["Patient_ID"]: r for r in csv.DictReader(
                open(os.path.join(args.data_dir, per_patient), encoding="utf-8-sig"))}
            kk = [k for k, pid in enumerate(i2)
                  if pid in s5 and s5[pid]["Included_in_primary_validation"] == "True"]
            tt = np.array([float(s5[i2[k]]["OS_months"]) for k in kk])
            ee = np.array([float(s5[i2[k]]["Death"]) for k in kk], dtype=int)
        for log2 in (False, True):
            a = score(e2, {g: coef4[g] for g in GENES4}, log2=log2)
            b = score(e2, coef25, log2=log2)
            ca = concordance_index(tt, -a[kk], ee)
            cb = concordance_index(tt, -b[kk], ee)
            rows.append([label, "log2(RSEM+1)" if log2 else "raw RSEM", len(kk), int(ee.sum()),
                         round(ca, 4), round(cb, 4), round(ca - cb, 4)])
    with open(s36, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["Cohort", "Expression scale", "n", "events", "C_MRG4", "C_ModelA", "Difference"])
        w.writerows(rows)

    summary = {
        "cohort": "CGGA-693, WHO grade II-III, independent of CGGA-325",
        "n_profiled": len(ids),
        "n": n,
        "events": int(event.sum()),
        "median_os_months": float(np.median(time)),
        "n_grade2": int(sum(clin[ids[i]]["Grade"] == "WHO II" for i in keep)),
        "n_grade3": int(sum(clin[ids[i]]["Grade"] == "WHO III" for i in keep)),
        "n_idh_mutant": int(sum(clin[ids[i]]["IDH_mutation_status"] == "Mutant" for i in keep)),
        "n_idh_missing": int(sum(clin[ids[i]]["IDH_mutation_status"] not in ("Mutant", "Wildtype")
                                for i in keep)),
        "c_mrg4": float(c4),
        "c_modelA": float(c25),
        "ci_mrg4": ci(b4),
        "ci_modelA": ci(b25),
        "difference": float(c4 - c25),
        "ci_difference": ci(bd),
        "p_difference": float(2 * min((np.array(bd) <= 0).mean(), (np.array(bd) >= 0).mean())),
        "sd_score": float(np.std(s4, ddof=1)),
        "exclusions": {k: reasons.count(k) for k in sorted(set(reasons))},
    }
    with open(os.path.join(args.out_dir, "replication_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    # ------------------------------------------------- reproduction assertions
    # Published values (article, Table 2 and Figure 4A). The bootstrap intervals
    # move in the third decimal with the seed, so only the point estimates and the
    # analysis-set sizes are asserted.
    expected = {"n": 420, "events": 197, "c_mrg4": 0.6736, "c_modelA": 0.6475}
    bad = [f"{k}: got {summary[k]!r}, expected {v!r}"
           for k, v in expected.items() if abs(float(summary[k]) - float(v)) > 5e-4]
    if bad:
        raise SystemExit("REPRODUCTION CHECK FAILED:\n  " + "\n  ".join(bad))

    print(f"analysis set      : {n} patients, {int(event.sum())} deaths, "
          f"median OS {summary['median_os_months']:.1f} months")
    print(f"MRG-4             : C = {c4:.4f} ({summary['ci_mrg4'][0]:.3f}-{summary['ci_mrg4'][1]:.3f})")
    print(f"Model A           : C = {c25:.4f} ({summary['ci_modelA'][0]:.3f}-{summary['ci_modelA'][1]:.3f})")
    print(f"difference        : {c4 - c25:+.4f} "
          f"({summary['ci_difference'][0]:+.3f} to {summary['ci_difference'][1]:+.3f})")
    print(f"scale sensitivity : "
          + "; ".join(f"{r[0]} {r[1]} {r[6]:+.4f}" for r in rows))
    print("wrote S34, S35, S36 and replication_summary.json to " + args.out_dir)
    print("REPRODUCTION CHECK: all green (4 assertions)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
