# Changelog

## v1.0.1 — 2026-11-11

Adds the replication analysis that the article now reports, and repairs three
stale cross-references left in a deposited table.

### Added

- `code/analysis/replication_cgga693.py` — builds the second, non-overlapping
  CGGA cohort (CGGA-693, WHO grade II–III), scores both frozen indices, and writes
  Supplementary Tables S34–S36 plus `replication_summary.json`. The script asserts
  the published C-indices to 5 × 10⁻⁴ and stops if they are not reproduced.
- `code/analysis/transport_diagnostics.py` — calibration slope, CITL (with patient
  bootstrap interval), ICI at 36 months, Schoenfeld residual test, grade-stratified
  C-indices in both cohorts, and the difference-of-differences test.
- `code/figures/fig4_second_cohort.py` — replica of Figure 4 with six layout
  assertions (no vertical overlap between panels, panel letter on its own line
  above the title, titles top-aligned, no text clipped by the canvas).
- `data/Table_S34_CGGA693_per_patient.csv` — every profiled patient in CGGA-693
  with the inclusion flag and the exclusion reason, plus gene *z*-scores, both
  linear predictors, survival and covariates. 694 rows.
- `data/Table_S35_calibration_CGGA693.csv` — frozen-model calibration by risk
  tertile, 36 and 60 months, with Kaplan–Meier 95% intervals.
- `data/Table_S36_scale_sensitivity.csv` — raw RSEM versus log2(RSEM+1) in both
  CGGA cohorts.
- `replication/` — the JSON records from which every replication number in the
  article is taken.
- `requirements.txt` — the exact package versions used for this release.

### Fixed

- `data/Table_S30_不确定性.csv` carried three cross-references to the numbering of
  an earlier manuscript draft. `Table 3` (hazard ratios) is now `Table 2`;
  `Table 4` (the uncertainty table itself) is now `Supplementary Table S32`;
  `Figure 3B` (the 26-gene apparent C-index) is now `Figure 2B`. No value in the
  table changed.

### Changed

- `README.md` — v1.0.1 release line, the CGGA-693 replication section, the
  reproduction recipe for the three new scripts, and a scope note that two CGGA
  cohorts are not two independent resources.
- `CITATION.cff` — version 1.0.1, release date, keyword `external validation`.
- `SHA256SUMS.txt`, `data/SHA256SUMS.txt`, `data/_补充表清单.md` — regenerated over
  the new file set.

## v1.0.0 — 2026-10-05

First release: derivation of Model A (25-gene ridge Cox) and MRG-4 (four-gene
stability-selected score), external validation in CGGA-325, single-cell
localisation, mutation and immune analyses, PROBAST self-assessment, and the
deposited derived tables S1–S30.
