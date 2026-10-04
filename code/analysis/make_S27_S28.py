# -*- coding: utf-8 -*-
"""R3-M2 / R3-M5 关闭：Table_S27（PROBAST 自评）+ Table_S28（随机种子与软件环境）

两表均被 manuscript_v46.md 正文引用，必须在送第四轮盲审前落盘。
"""
import sys, platform, importlib.metadata as md
import pandas as pd

OUT = 'Table_S27_PROBAST_self_assessment.csv'
OUT2 = 'Table_S28_random_seeds_and_environment.csv'

# ---------------------------------------------------------------- S27 PROBAST
# 口径：PROBAST 域内任一 signalling question 判为 "Probably no"/"No" → 该域 High；
#       全部 "Yes"/"Probably yes" → Low；出现 "No information" → Unclear。
DOM = [
    # (域, 问题编号, 问题, 判定, 理由)
    ('1 Participants', '1.1', 'Appropriate data sources (cohort, RCT, nested case-control)',
     'Yes',
     'TCGA-LGG (population-based cohort, n=514 with mRNA) and CGGA-325 (independent cohort, n=325) are both public, consecutive series. Neither is a case-control or a convenience subset.'),
    ('1 Participants', '1.2', 'All inclusions and exclusions appropriate',
     'Yes',
     'Four of 514 TCGA patients were excluded for absent expression or absent survival; the external set was restricted a priori to WHO grade II-III to match the discovery population (all low-grade). Every exclusion is itemised per patient in Supplementary Tables S21 and S22.'),
    ('1 Participants', 'JUDGEMENT', 'Domain 1 - Participants', 'LOW',
     'All signalling questions answered Yes.'),

    ('2 Predictors', '2.1', 'Predictors defined and assessed in a similar way for all participants',
     'Probably no',
     'The four genes are the predictors, but their values are standardised against two different reference distributions: the cBioPortal study profile (..._Zscores) in the discovery cohort, and each cohort\'s own distribution in external validation. The same nominal predictor therefore enters on a different scale in the two cohorts. This is why the manuscript carries no absolute cross-cohort risk claim and reports the score only as a relative, within-cohort quantity (Sections 2.12 and 3.9).'),
    ('2 Predictors', '2.2', 'Predictor assessment made without knowledge of outcome data',
     'Yes',
     'Expression was measured by the depositing consortia independently of follow-up; no outcome information entered the assay, the standardisation or the gene values.'),
    ('2 Predictors', '2.3', 'All predictors available at the time the model is intended to be used',
     'Yes',
     'All four are bulk RNA-seq values available at diagnosis from the same specimen used for routine grading.'),
    ('2 Predictors', 'JUDGEMENT', 'Domain 2 - Predictors', 'HIGH',
     'Q2.1 is Probably no: the covariates are not on a common standardisation across the development and validation cohorts.'),

    ('3 Outcome', '3.1', 'Outcome determined appropriately',
     'Yes',
     'Overall survival, taken from the clinical files of the two depositing studies.'),
    ('3 Outcome', '3.2', 'Prespecified or standard outcome definition used',
     'Yes',
     'Death from any cause, measured from diagnosis; the same definition in both cohorts, requiring no adjudication by the investigators.'),
    ('3 Outcome', '3.3', 'Predictors excluded from the outcome definition',
     'Yes',
     'Gene expression plays no part in the survival endpoint.'),
    ('3 Outcome', '3.4', 'Outcome defined and determined in a similar way for all participants',
     'Yes',
     'Identical definition, identical time origin, in both cohorts.'),
    ('3 Outcome', '3.5', 'Time to outcome or time of prediction appropriate',
     'Yes',
     'Prediction at diagnosis; follow-up is adequate in both cohorts (CGGA-325 median 130.2 months; TCGA-LGG 125 deaths among 510).'),
    ('3 Outcome', '3.6', 'If outcome taken from existing records, were assessors blinded to predictors',
     'Probably yes',
     'Deaths were recorded by the registry, not collected or adjudicated by the investigators, who had no contact with participants.'),
    ('3 Outcome', 'JUDGEMENT', 'Domain 3 - Outcome', 'LOW',
     'No signalling question answered No or Probably no.'),

    ('4 Analysis', '4.1', 'Reasonable number of participants with the outcome',
     'No for Model A / Yes for MRG-4',
     'Model A: 125 deaths for 25 predictors, 5.0 events per variable, below the conventional minimum of 10. MRG-4: 125 deaths for 4 predictors, 31.2 events per variable. The two models are therefore not equally reliable, and the paper reports the Model A comparator with that caveat attached.'),
    ('4 Analysis', '4.2', 'Continuous and categorical predictors handled appropriately',
     'Yes',
     'Expression standardised as described; the continuous score is additionally reported per SD of the discovery cohort, so that no hazard ratio depends on an arbitrary dichotomisation.'),
    ('4 Analysis', '4.3', 'All enrolled participants included in the analysis',
     'Probably yes',
     '510 of 514 TCGA patients and 172 of 325 CGGA patients; the two dropped TCGA patients and the 12 without survival time are identified individually in Supplementary Tables S21 and S22.'),
    ('4 Analysis', '4.4', 'Participants with missing data handled appropriately',
     'Yes',
     'Complete-case analysis, with the derivation of every analysis set deposited per patient (Supplementary Tables S21 and S22) rather than described in prose only.'),
    ('4 Analysis', '4.5', 'Selection of predictors based on univariable analysis avoided',
     'No',
     'Model A\'s 25 predictors were taken from a univariable screen at P < 0.05, and ten of the 25 then carry a coefficient opposite in sign to that screen (Table 2). MRG-4 is derived by LASSO with stability selection rather than by a univariable filter, but its candidate set still entered through the same screen. This is the study\'s clearest source of analysis bias and is reported as a result, not only as a caveat.'),
    ('4 Analysis', '4.6', 'Complexities in the data accounted for appropriately',
     'Probably yes',
     'Censoring is handled by the Cox model throughout and by the concordance and RMST estimators in validation. Competing risks (death without glioma progression) are not modelled; absolute-risk statements are accordingly framed as life-table quantities (RMST), not as cause-specific probabilities.'),
    ('4 Analysis', '4.7', 'Relevant model performance measures evaluated appropriately',
     'No',
     'Discrimination was evaluated thoroughly (apparent, cross-validated, fully nested and external C-indices; time-dependent AUC with a patient-level paired bootstrap; RMST differences). Calibration and clinical utility (net benefit, decision curves) were not evaluated, and no comparison against an externally published model was possible because no comparable published model of these MRGs existed at the time.'),
    ('4 Analysis', '4.8', 'Overfitting, underfitting and optimism accounted for',
     'Yes',
     'The apparent C-index, the cross-validated value and the fully nested value are reported separately and explicitly distinguished, and the internal comparison between the two models is declared uninformative rather than quoted as evidence.'),
    ('4 Analysis', '4.9', 'Predictors and weights correspond to the reported multivariable analysis',
     'Yes',
     'The four coefficients are frozen from the multivariable fit and deposited (Table S12). For Model A the discrepancy between the archived score and a score rebuilt from the deposited coefficients is quantified rather than hidden (r = 0.9148; Supplementary Table S20).'),
    ('4 Analysis', 'JUDGEMENT', 'Domain 4 - Analysis', 'HIGH',
     'Q4.1 No for Model A, Q4.5 No, Q4.7 No.'),

    ('Overall', '-', 'Overall risk of bias judgement', 'HIGH',
     'Two of the four domains are at high risk of bias (Predictors and Analysis). The overall assessment is therefore High. This is stated in the manuscript rather than left to the reader, and the consequences are quantified in the corresponding sections: the external comparison is the one to carry, and the absolute score is not transferable between platforms.'),
    ('Overall', '-', 'Applicability concern', 'LOW',
     'Participants (adult diffuse low-grade glioma), predictors (bulk RNA-seq at diagnosis) and outcome (overall survival) all match the intended clinical question. Applicability was therefore not judged to be a concern.'),
]
pd.DataFrame(DOM, columns=['Domain', 'Question', 'Signalling_question', 'Answer', 'Justification']).to_csv(
    OUT, index=False, encoding='utf-8-sig')
print('-> %s  %d 行' % (OUT, len(DOM)))

# ------------------------------------------------- S28 随机种子与软件环境
def ver(p):
    try:
        return md.version(p)
    except Exception:
        return 'not installed'

ENV = [
    ('Environment', 'Operating system', 'platform', platform.platform(), 'Apple silicon (arm64)'),
    ('Environment', 'Python', 'interpreter', sys.version.split()[0], 'CPython'),
    ('Environment', 'numpy', 'installed version', ver('numpy'), 'array and linear algebra'),
    ('Environment', 'scipy', 'installed version', ver('scipy'), 'optimisation, interpolation'),
    ('Environment', 'pandas', 'installed version', ver('pandas'), 'tables'),
    ('Environment', 'scikit-learn', 'installed version', ver('scikit-learn'), 'KFold, KMeans, StandardScaler'),
    ('Environment', 'scikit-survival', 'installed version', ver('scikit-survival'), 'CoxnetSurvivalAnalysis, concordance_index_censored, cumulative_dynamic_auc'),
    ('Environment', 'lifelines', 'installed version', ver('lifelines'), 'CoxPHFitter, Breslow baseline cumulative hazard'),
    ('Environment', 'matplotlib', 'installed version', ver('matplotlib'), 'figure rendering'),
    ('Environment', 'anndata', 'installed version', ver('anndata'), 'single-cell h5ad reading'),
    ('Environment', 'h5py', 'installed version', ver('h5py'), 'HDF5 access'),
    ('Environment', 'openpyxl', 'installed version', ver('openpyxl'), 'spreadsheet export'),
]

SEED = [
    ('Seed', 'k-means subtypes of the 25 Model A genes (k = 2, Figure 8A)',
     'sklearn.cluster.KMeans(n_clusters=2, n_init=50, random_state=42)', '42',
     'n_init = 50 restarts; the partition is stable across restarts (silhouette 0.276).'),
    ('Seed', '10-fold cross-validation along the LASSO path (Figure 4A, Table S7)',
     'sklearn.model_selection.KFold(n_splits=10, shuffle=True, random_state=7)', '7',
     'Shuffled folds; the curve is flat-topped, so the fold assignment is not decisive.'),
    ('Seed', 'Fully nested cross-validation (the 0.731 conservative estimate)',
     'outer KFold(random_state=7); inner alpha selection seed = 1', '7 (outer), 1 (inner)',
     'Gene selection, alpha selection and performance estimation are all inside the outer loop.'),
    ('Seed', 'Sensitivity of the cross-validated estimate to fold assignment',
     'KFold(random_state=42), KFold(random_state=11)', '42, 11',
     'Reported only as a robustness check; no claim rests on the difference.'),
    ('Seed', '200-bootstrap stability selection, original run (Table 7, Table S6)',
     'numpy.random.seed(20261004)', '20261004',
     '200 patient-level bootstrap resamples of the 88-gene candidate set at alpha = 0.02432298416536324.'),
    ('Seed', '200-bootstrap stability selection, replication run (Table 7, Table S6)',
     'numpy.random.default_rng(7)', '7',
     'Independent seed over the full 89-gene candidate set; the candidate sets differ, so this is a check that the same core is recovered, not a pure seed-dependence test.'),
    ('Seed', 'Randomised-panel benchmark, run 1 (Table 10)',
     'numpy.random.default_rng(20261004 + k + int(penalty * 100))', '20261008 (k = 4), 20261129 (k = 25)',
     '3,000 draws per panel size; the seed is a deterministic function of k and the penalty so that each panel is independently reproducible.'),
    ('Seed', 'Randomised-panel benchmark, run 2 (Table 10)',
     'as above + 977', '20261985 (k = 4), 20262106 (k = 25)',
     'Independent replication of the same benchmark; a conclusion is drawn only where it holds in both runs.'),
    ('Seed', 'NRI and IDI (Section 3.6)',
     'numpy.random.seed(42)', '42',
     'Patient-level resampling for the continuous-NRI and IDI confidence intervals.'),
    ('Seed', 'Time-dependent AUC and RMST, paired bootstrap (Figures 6C, 6D; Table S23, Table S24)',
     'numpy.random.default_rng(20261004) for the shared patient resample, numpy.random.default_rng(7) for the paired draw', '20261004, 7',
     '2,000 resamples; the same resampled patients are scored by both models, so the difference is paired rather than independent.'),
    ('Seed', 'Single-cell aggregation (Figures 9A, 9B; Tables S10, S11, S16)',
     'none - deterministic', '-',
     'Cell-type and state means are computed over the full deposited cell set; no subsampling or stochastic step is involved.'),
]

rows = [[t, a, b, c, d] for (t, a, b, c, d) in ENV] + [[t, a, b, c, d] for (t, a, b, c, d) in SEED]
pd.DataFrame(rows, columns=['Type', 'Item', 'Setting', 'Value', 'Note']).to_csv(
    OUT2, index=False, encoding='utf-8-sig')
print('-> %s  %d 环境行 + %d 种子行' % (OUT2, len(ENV), len(SEED)))
