# -*- coding: utf-8 -*-
"""R3-M1 / R1-m7 / R1-m2 关闭
A. MRG-4 部署包：常数表 + 逐患者四基因标准值 + Breslow 基线累积风险
B. Model A 线性预测子溯源（存档版 vs 由 S12 系数重导出）
C. 分析集推导表（510 / 502 / 410；325 / 313 / 182 / 172）
"""
import json, io, numpy as np, pandas as pd, urllib.request
from lifelines import CoxPHFitter
from sksurv.metrics import concordance_index_censored as cic

C = lambda t, e, r: cic(np.asarray(e).astype(bool), np.asarray(t), np.asarray(r))[0]
B = 'https://www.cbioportal.org/api'
def post(u, p):
    r = urllib.request.Request(u, data=json.dumps(p).encode(),
        headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
    return json.load(urllib.request.urlopen(r, timeout=300))

G4 = ['BMP1', 'KIF15', 'TRAF3', 'PDGFA']
COEF = {'BMP1': 0.2581639216398504, 'KIF15': 0.4635562109878421,
        'TRAF3': -0.23788676758461966, 'PDGFA': 0.4795327193081988}
PROF_Z = 'lgg_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_all_sample_Zscores'
PROF_RAW = 'lgg_tcga_pan_can_atlas_2018_rna_seq_v2_mrna'
SL = 'lgg_tcga_pan_can_atlas_2018_rna_seq_v2_mrna'

E = pd.read_csv('panel_expr_tcga.csv', index_col=0, encoding='utf-8-sig')
G = pd.read_csv('lgg_g1_dataset.csv', encoding='utf-8-sig', low_memory=False)
d = G.merge(E, left_on='pat', right_index=True, how='inner').rename(columns={'e': 'ev', 't': 'time'})
d['lp'] = sum(d[g] * COEF[g] for g in G4)

# ---------- A1. 原始表达 vs z 谱的关系（如实报告：非同一分布） ----------
g = post(B + '/genes/fetch?geneIdType=HUGO_GENE_SYMBOL', G4)
inv = {x['entrezGeneId']: x['hugoGeneSymbol'] for x in g}
def fetch(prof):
    rows = post(f'{B}/molecular-profiles/{prof}/molecular-data/fetch?projection=SUMMARY',
                {'entrezGeneIds': sorted(inv), 'sampleListId': SL})
    return pd.DataFrame([{'upk': x.get('uniquePatientKey'), 'gene': inv[x['entrezGeneId']],
                          'v': x['value']} for x in rows]).pivot_table(index='upk', columns='gene', values='v')
RAW, ZP = fetch(PROF_RAW), fetch(PROF_Z)
cc = RAW.index.intersection(ZP.index); RAW, ZP = RAW.loc[cc], ZP.loc[cc]
rel = {}
for gname in G4:
    x, z = RAW[gname].values, ZP[gname].values
    rel[gname] = dict(raw_median=float(np.median(x)), raw_mad=float(np.median(np.abs(x - np.median(x)))),
                      pearson_r=float(np.corrcoef(x, z)[0, 1]),
                      spearman_r=float(pd.Series(x).corr(pd.Series(z), method='spearman')))
    print('%-6s 原始与 z 谱：Pearson r=%.4f  Spearman %.4f' % (gname, rel[gname]['pearson_r'], rel[gname]['spearman_r']))

# ---------- A2. 常数表 ----------
sd_t = float(d['lp'].std(ddof=1))
cg = pd.read_csv('cgga_mrg4_perpatient.csv', encoding='utf-8-sig')
sd_c = float(cg['lp'].std(ddof=1))
rows = []
for gname in G4:
    rows.append(['系数', f'{gname} (标准值)', '%.10f' % COEF[gname], '对 cBioPortal z 谱值加权求和'])
rows += [
    ['表达谱（标准值）', 'z-score profile', PROF_Z, '模型直接定义在该公共已沉积对象上；作者未自估标准化常数'],
    ['表达谱（原始）', 'raw profile', PROF_RAW, '用于检验两者关系；与 z 谱非同一分布（见下）'],
    ['样本清单', 'sample list', SL, 'cBioPortal 该研究下全部有 mRNA 数据的样本'],
    ['评分均值（发现队列）', 'TCGA-LGG n=510', '%.6f' % d['lp'].mean(), '未标准化评分的中心'],
    ['评分 SD（发现队列）', 'TCGA-LGG n=510', '%.6f' % sd_t, 'per-SD 分母；与正文 0.903 一致'],
    ['评分 SD（验证队列）', 'CGGA-325 WHO II-III n=172', '%.6f' % sd_c, '与正文 0.591 一致'],
    ['基线累积风险', 'H0(t)', 'Table S19', 'Breslow 估计，取自 lp 中心化后的冻结模型（510 例）'],
]
for gname in G4:
    rows.append(['原始/z 关系', f'{gname} 的 Pearson r', '%.4f' % rel[gname]['pearson_r'],
                 '两者非同一分布，故跨平台的绝对风险不可迁移（正文已声明只用相对量）'])
pd.DataFrame(rows, columns=['item', 'detail', 'value', 'note']).to_csv(
    'Table_S17_MRG4_deployment_constants.csv', index=False, encoding='utf-8-sig')
print('-> Table_S17  %d 行' % len(rows))

# ---------- A3. 逐患者四基因标准值 ----------
tc = d[['pat'] + G4 + ['lp', 'time', 'ev']].copy()
tc = tc.rename(columns={'pat': 'Patient_ID', 'lp': 'MRG4_score', 'time': 'OS_months', 'ev': 'Death'})
tc.columns = ['Patient_ID'] + [f'{x}_z' for x in G4] + ['MRG4_score', 'OS_months', 'Death']
tc.insert(1, 'Cohort', 'TCGA-LGG')
tc.insert(2, 'Convention', 'cBioPortal z-score profile')

cexpr = pd.read_csv('cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt', sep='\t', index_col=0)
_raw = open('cgga/CGGA.mRNAseq_325_clinical.20200506.txt', encoding='utf-8-sig', errors='ignore').read().replace('\r', '\n')
cl = pd.read_csv(io.StringIO(_raw), sep='\t'); cl.columns = [c.strip() for c in cl.columns]
cl = cl.set_index('CGGA_ID')
cl['ev'] = pd.to_numeric(cl['Censor (alive=0; dead=1)'], errors='coerce')
cl['mos'] = pd.to_numeric(cl['OS'], errors='coerce') / 30.4375
M = cexpr.loc[G4].T; M = M[~M.index.duplicated()]
Zc = (M.values - M.values.mean(0)) / (M.values.std(0) + 1e-9)
vc = pd.DataFrame(Zc, index=M.index, columns=[f'{x}_z' for x in G4])
vc['MRG4_score'] = vc[[f'{x}_z' for x in G4]].values @ np.array([COEF[x] for x in G4])
vc = vc.join(cl[['mos', 'ev', 'Grade']]).dropna(subset=['MRG4_score', 'mos', 'ev'])
vc = vc[vc['Grade'].isin(['WHO II', 'WHO III'])].copy()
vc = vc.reset_index().rename(columns={'index': 'Patient_ID', 'mos': 'OS_months', 'ev': 'Death'})
vc.insert(1, 'Cohort', 'CGGA-325 WHO II-III')
vc.insert(2, 'Convention', 'cohort-internal z (n=325 reference)')
vc = vc[['Patient_ID', 'Cohort', 'Convention'] + [f'{x}_z' for x in G4] + ['MRG4_score', 'OS_months', 'Death']]
both = pd.concat([tc, vc], ignore_index=True)
both.to_csv('Table_S18_MRG4_per_patient_gene_values.csv', index=False, encoding='utf-8-sig')
print('-> Table_S18  TCGA %d 行 + CGGA %d 行 = %d' % (len(tc), len(vc), len(both)))
print('   冻结常数复现：TCGA C=%.4f（发表 0.7938）| CGGA C=%.4f（发表 0.7688）'
      % (C(tc['OS_months'], tc['Death'], tc['MRG4_score']),
         C(vc['OS_months'], vc['Death'], vc['MRG4_score'])))

# ---------- A4. Breslow 基线累积风险 ----------
mean_lp = float(d['lp'].mean())
d['lp0'] = d['lp'] - mean_lp
f = CoxPHFitter().fit(d[['time', 'ev', 'lp0']], 'time', 'ev')
b = f.baseline_cumulative_hazard_
bh = pd.DataFrame({'time_months': b.index.values, 'H0_cumhaz': b.iloc[:, 0].values})
bh['S0'] = np.exp(-bh['H0_cumhaz'])
bh['is_event_time'] = True
grid = [12, 24, 36, 48, 60, 72, 84, 96, 120]
H0 = lambda t: float(np.interp(t, bh['time_months'].values, bh['H0_cumhaz'].values))
extra = pd.DataFrame({'time_months': grid, 'H0_cumhaz': [H0(t) for t in grid],
                      'S0': [np.exp(-H0(t)) for t in grid], 'is_event_time': False})
bh = pd.concat([bh, extra], ignore_index=True).drop_duplicates('time_months').sort_values('time_months')
bh['S0'] = np.exp(-bh['H0_cumhaz'])
bh.to_csv('Table_S19_baseline_cumulative_hazard.csv', index=False, encoding='utf-8-sig')
print('-> Table_S19  %d 行（其中 %d 个事件时点）| H0(60)=%.4f' % (len(bh), int(bh['is_event_time'].sum()), H0(60)))

json.dump({'mean_lp_TCGA': mean_lp, 'sd_lp_TCGA': sd_t, 'sd_lp_CGGA172': sd_c,
           'H0_grid': {str(t): round(H0(t), 6) for t in grid},
           'raw_vs_z': rel, 'n_reference': int(len(RAW))},
          open('deploy_mrg4_context.json', 'w'), ensure_ascii=False, indent=1)

# ---------- B. Model A 线性预测子溯源 ----------
S12 = pd.read_csv('Table_S12_ModelA_coefficients.csv', encoding='utf-8-sig')
print('\nS12 列:', list(S12.columns))
gcol = [c for c in S12.columns if '基因' in c or 'gene' in c.lower()][0]
ccol = [c for c in S12.columns if '系数' in c or 'coef' in c.lower()][0]
co25 = dict(zip(S12[gcol].astype(str), S12[ccol].astype(float)))
have = [x for x in co25 if x in E.columns]
Za = (E[have].values - E[have].values.mean(0)) / (E[have].values.std(0) + 1e-9)
lpA_re = pd.Series(Za @ np.array([co25[x] for x in have]), index=E.index)
chk = d[['pat', 'risk']].copy(); chk['lpA_re'] = chk['pat'].map(lpA_re)
chk = chk.dropna()
r_arch = float(np.corrcoef(chk['risk'], chk['lpA_re'])[0, 1])
r_arch_s = float(pd.Series(chk['risk'].values).corr(pd.Series(chk['lpA_re'].values), method='spearman'))
z_direct = pd.Series(E[have].values @ np.array([co25[x] for x in have]), index=E.index)
chk['lpA_zdir'] = chk['pat'].map(z_direct)
r_dir = float(np.corrcoef(chk['risk'], chk['lpA_zdir'])[0, 1])
print('存档 Model A LP vs 由 S12 系数在原矩阵上重导出（先 z 后加权）：r=%.4f (Spearman %.4f)' % (r_arch, r_arch_s))
print('存档 Model A LP vs 直接加权已沉积 z 值                   ：r=%.4f' % r_dir)
print('存档 LP 的 C-index = %.4f' % C(d['time'], d['ev'], d['risk']))
AUD = pd.DataFrame([
    ['archived_ModelA_LP', 'n=%d' % len(chk), 'r_vs_recompute_zfirst=%.6f' % r_arch,
     'Spearman %.6f' % r_arch_s],
    ['archived_ModelA_LP', 'n=%d' % len(chk), 'r_vs_recompute_weighting_deposited_z=%.6f' % r_dir,
     '直接对已沉积 z 值加权'],
], columns=['quantity', 'n', 'check', 'note'])
AUD.to_csv('Table_S20_ModelA_LP_provenance.csv', index=False, encoding='utf-8-sig')
print('-> Table_S20')

# ---------- C. 分析集推导 ----------
S4 = pd.read_csv('Table_S4_TCGA_LGG_MRG4_per_patient.csv', encoding='utf-8-sig')
keycol = 'uniquePatientKey'
mrg = pd.read_csv('mrg4_tcga_perpatient.csv', encoding='utf-8-sig')
mrg = mrg.rename(columns={'pat': keycol})
der = S4[[keycol, 'OS_months', 'Death', 'Age', 'Sex', 'Molecular_subtype', 'WHO_grade']].copy()
der['has_expression'] = der[keycol].isin(E.index)
der['has_OS'] = der['OS_months'].notna() & der['Death'].notna()
sub = der['Molecular_subtype'].astype(str)
der['has_IDH_and_1p19q'] = sub.str.startswith('LGG_IDH')
der['has_grade'] = der['WHO_grade'].notna()
der['in_510'] = der['has_expression'] & der['has_OS']
der['in_502'] = der['in_510'] & der['has_IDH_and_1p19q'] & der['has_grade']
der['in_410'] = der['in_502'] & (sub.str.contains('IDHmut'))
excl = []
for _, r in der.iterrows():
    if not r['in_502']:
        why = []
        if not r['in_510']: why.append('no expression or no survival')
        if not r['has_IDH_and_1p19q']: why.append('missing IDH/1p19q subtype')
        if not r['has_grade']: why.append('missing WHO grade')
        excl.append('; '.join(why))
    else:
        excl.append('')
der['exclusion_reason_from_502'] = excl
der['analysis_set'] = np.where(der['in_410'], '502+410', np.where(der['in_502'], '502', np.where(der['in_510'], '510', 'excluded')))
der.to_csv('Table_S21_analysis_set_derivation_TCGA.csv', index=False, encoding='utf-8-sig')
print('\n-> Table_S21  TCGA 510=%d  502=%d  410=%d' % (der['in_510'].sum(), der['in_502'].sum(), der['in_410'].sum()))

raw_cl = cl.copy()
# 注意：cg_der 必须以 CGGA_ID 为索引，否则赋值时与 raw_cl 的索引对不齐，
# notna() 结果会被 pandas 对齐成整列 NaN（float），后续 ~ 取反即报
# TypeError: bad operand type for unary ~: 'float'。
cg_der = pd.DataFrame({'Patient_ID': raw_cl.index}, index=raw_cl.index)
cg_der['has_OS'] = raw_cl['OS'].notna()
cg_der['grade_II_III'] = raw_cl['Grade'].isin(['WHO II', 'WHO III'])
cg_der['in_313'] = cg_der['has_OS']
cg_der['in_172'] = cg_der['has_OS'] & cg_der['grade_II_III']
cg_der['exclusion_reason'] = np.where(~cg_der['has_OS'], 'missing overall-survival time',
                              np.where(cg_der['grade_II_III'], '', 'grade IV (outside primary validation set)'))
cg_der.to_csv('Table_S22_analysis_set_derivation_CGGA.csv', index=False, encoding='utf-8-sig')
print('-> Table_S22  CGGA 325 原始=%d  313=%d  182(WHO II-III)=%d  172=%d'
      % (len(cg_der), cg_der['in_313'].sum(), ((raw_cl['Grade'].isin(['WHO II', 'WHO III'])).sum()),
         cg_der['in_172'].sum()))
