# MRG-4 and Model A — analysis code and derived data

Code and derived data for the manuscript on a mitoxyperilysis-related gene (MRG)
prognostic score in lower-grade glioma, submitted to *Gene*.

- **Author:** Bin Lian (ORCID [0000-0002-1477-9137](https://orcid.org/0000-0002-1477-9137))
- **Affiliation:** School of Health, Guangzhou Vocational and Technical University of Science and Technology, No. 1038 Guangcong 9th Road, Zhongluotan, Baiyun District, Guangzhou 510550, Guangdong, China
- **Correspondence:** drmilo@gkd.edu.cn
- **Release:** v1.0.0 — 2026-10-05

## What is here

| Path | Contents |
|---|---|
| `code/analysis/` | Model derivation, stability selection, external validation, statistics |
| `code/figures/` | Figure panel rendering, composite assembly and overlap/overflow checks |
| `code/deposits/` | Builders for the deposited derived tables |
| `data/` | The deposited derived tables, as supplied as Supplementary Tables, with their own `data/SHA256SUMS.txt` |
| `SHA256SUMS.txt` | SHA-256 of every file in this release |

## Primary data are not redistributed here

All primary data are public and were retrieved by the authors:

- **TCGA-LGG** — cBioPortal study `lgg_tcga_pan_can_atlas_2018` (RNA-Seq V2 and clinical data).
- **CGGA-325** — Chinese Glioma Genome Atlas, <http://www.cgga.org.cn> (RSEM and clinical annotation).
- **Single-cell** — CZ CELLxGENE Discover census, dataset
  `21722308-5091-4b63-9c07-4f116ab9a7b3` (SL040), glioblastoma, grade IV.

`data/` therefore holds **derived** quantities only — per-patient scores, fitted
coefficients, analysis-set logs, calibration and utility summaries. It contains no
raw expression matrix and no controlled-access data.

## Analysis sets and the two models

- **Model A** — 25-gene ridge (L2, lambda = 1.0) Cox index, fitted on the full
  89-gene MRG panel in TCGA-LGG. Apparent C = 0.8211; external CGGA-325 C = 0.6637;
  SD of the linear predictor in TCGA-LGG = 1.062 (n = 502).
- **Model B (MRG-4)** — four genes frozen by LASSO-Cox plus stability selection
  (pi = 0.75): `+0.258164 BMP1 + 0.463556 KIF15 - 0.237887 TRAF3 + 0.479533 PDGFA`.
  TCGA-LGG apparent 0.794; SD = 0.903031; external CGGA-325 C = 0.769
  (n = 172, 90 deaths).

Every per-SD hazard ratio is per one SD of **that model's own score** in that
model's own analysis set.

## Environment

Python 3.13 with scikit-survival 0.28.0 and lifelines 0.30.3. Exact versions,
random seeds and the software environment are recorded in `data/Table_S28_*.csv`.

## Reproducing

1. Download the three primary datasets listed above.
2. Edit the input paths at the top of the scripts in `code/analysis/` — they point
   to the authors' local working directory and are not portable as shipped.
3. Run `code/analysis/` in the order given by the script names; the deposit
   builders in `code/deposits/` regenerate `data/` from the analysis outputs.
4. Verify integrity with `shasum -a 256 -c SHA256SUMS.txt`.

Figure scripts assume the composite panels are assembled with a vector editor;
`code/figures/check_*.py` perform the automated overlap, overflow and occlusion
checks reported for each composite figure.

## Licence

- Code (`code/`): MIT — see `LICENSE`.
- Derived data (`data/`): CC BY 4.0 — see `LICENSE-DATA.md`.

## Citation

If you use this code or these derived tables, please cite the article (see
`CITATION.cff`) and this release.

## Scope note

This is a computational re-analysis of public cohorts. The score is a prognostic
association, not a validated clinical test; no wet-lab validation was performed,
and the "mitoxyperilysis-related" label is inherited from the screening panel
rather than from any measurement of lysis. Please read the limitations section of
the article before reusing the score.
