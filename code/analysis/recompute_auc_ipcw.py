# -*- coding: utf-8 -*-
"""重算时间依赖 AUC：把 IPCW 的删失分布从"发现队列"换成"被评估队列自身"。

背景（第四轮盲审 + 独立复核双确认）
  round4_addendum.py 里写的是
      cumulative_dynamic_auc(tr, te, ...)   # tr = TCGA 发现队列, te = CGGA 172
  即用**发现队列**的删失分布给**外部队列**的 AUC 做逆概率删失加权。
  两队列随访深度差近 5 倍（反向 KM 中位随访 27.8 月 vs 130.2 月），
  权重对 CGGA 样本是严重错配的。

本脚本输出两套口径，最终以"被评估队列自身删失分布"为准：
  A 口径 = 发现队列删失分布（旧口径，保留供对照）
  B 口径 = 被评估队列自身删失分布（新口径，采纳）

产出
  Table_S23_timedep_AUC_paired_bootstrap.csv   （覆盖：列含两套均值口径）
  auc_ipcw_both.json                           （差值、CI、逐点排除零计数）
"""
import json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from sksurv.util import Surv
from sksurv.metrics import cumulative_dynamic_auc

RNG = np.random.default_rng(20261004)
TIMES = np.arange(12, 121, 6)               # 19 个时点，与旧表一致

# ---------------------------------------------------------------- 数据准备（同 round4_addendum.py）
E = pd.read_csv('panel_expr_tcga.csv', index_col=0, encoding='utf-8-sig')
G = pd.read_csv('lgg_g1_dataset.csv', encoding='utf-8-sig', low_memory=False)
S12 = pd.read_csv('Table_S12_ModelA_coefficients.csv', encoding='utf-8-sig')
co = dict(zip(S12['Gene'], S12['Multivariate_coefficient_ridge_L2']))
G25 = list(co)
C4 = {'BMP1': 0.2581639216398504, 'KIF15': 0.4635562109878421,
      'TRAF3': -0.23788676758461966, 'PDGFA': 0.4795327193081988}

cexpr = pd.read_csv('cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt', sep='\t', index_col=0)
_raw = open('cgga/CGGA.mRNAseq_325_clinical.20200506.txt', encoding='utf-8-sig',
            errors='ignore').read().replace('\r', '\n')
cl = pd.read_csv(__import__('io').StringIO(_raw), sep='\t')
cl.columns = [c.strip() for c in cl.columns]
cl = cl.set_index('CGGA_ID')
cl['ev'] = pd.to_numeric(cl['Censor (alive=0; dead=1)'], errors='coerce')
cl['mos'] = pd.to_numeric(cl['OS'], errors='coerce') / 30.4375
cl['age'] = pd.to_numeric(cl['Age'], errors='coerce')
cl['male'] = (cl['Gender'].astype(str).str[0] == 'M').astype(float)
cl['idhwt'] = cl['IDH_mutation_status'].astype(str).str.lower().str.startswith('wild').astype(float)
cl['grade_hi'] = cl['Grade'].isin(['WHO III', 'WHO IV']).astype(float)

ALL = G25 + [g for g in C4 if g not in G25]
M = cexpr.loc[ALL].T
Z = pd.DataFrame((M.values - M.values.mean(0)) / (M.values.std(0) + 1e-9),
                 index=M.index, columns=M.columns)
v = pd.DataFrame({'lp25': Z[G25].values @ np.array([co[g] for g in G25]),
                  'lp4': Z[list(C4)].values @ np.array([C4[g] for g in C4])}, index=M.index)
v = v.join(cl[['mos', 'ev', 'age', 'male', 'idhwt', 'grade_hi']])
v = v[v.index.isin([i for i in cl.index if cl.loc[i, 'Grade'] in ('WHO II', 'WHO III')])]
v = v.dropna(subset=['mos', 'ev', 'lp4', 'lp25'])
v['lp4z'] = (v['lp4'] - v['lp4'].mean()) / v['lp4'].std()
v['lp25z'] = (v['lp25'] - v['lp25'].mean()) / v['lp25'].std()
fcl = CoxPHFitter().fit(v[['mos', 'ev', 'age', 'male', 'idhwt', 'grade_hi']], 'mos', 'ev')
v['lpclinz'] = fcl.predict_log_partial_hazard(v[['age', 'male', 'idhwt', 'grade_hi']]).values
v['lpclinz'] = (v['lpclinz'] - v['lpclinz'].mean()) / v['lpclinz'].std()
print('验证队列 n=%d 死亡=%d' % (len(v), int(v['ev'].sum())))

TCGA_SURV = Surv.from_arrays(G['e'].astype(bool).values, G['t'].values)
CGGA_SURV = Surv.from_arrays(v['ev'].astype(bool).values, v['mos'].values)


def auc_set(df, survival_train):
    te = Surv.from_arrays(df['ev'].astype(bool).values, df['mos'].values)
    out = {}
    for k in ('lp4z', 'lp25z', 'lpclinz'):
        a, m = cumulative_dynamic_auc(survival_train, te, df[k].values, TIMES)
        out[k] = (float(m), np.asarray(a, dtype=float))
    return out


A = auc_set(v, TCGA_SURV)      # 旧口径：发现队列删失分布
Bs = auc_set(v, CGGA_SURV)     # 新口径：被评估队列自身删失分布


def unweighted(arr):
    return float(np.mean(arr))


print('\n口径A（发现队列删失分布） MRG-4 %.4f  Model A %.4f  Clinical %.4f' %
      (A['lp4z'][0], A['lp25z'][0], A['lpclinz'][0]))
print('口径B（评估队列删失分布） MRG-4 %.4f  Model A %.4f  Clinical %.4f' %
      (Bs['lp4z'][0], Bs['lp25z'][0], Bs['lpclinz'][0]))

# ---------------------------------------------------------------- 配对自举（两套口径各跑一次）
B = 2000
n = len(v)


def bootstrap(surv_train, seed_shift=0):
    rng = np.random.default_rng(20261004 + seed_shift)
    d25, dcl = [], []
    ok = 0
    for _ in range(B):
        idx = rng.integers(0, n, n)
        df = v.iloc[idx]
        if df['ev'].sum() < 8:
            continue
        try:
            s = auc_set(df, surv_train)
        except Exception:
            continue
        d25.append(s['lp4z'][0] - s['lp25z'][0])
        dcl.append(s['lp4z'][0] - s['lpclinz'][0])
        ok += 1
    return (np.array(d25), np.array(dcl), ok)


dA25, dAcl, okA = bootstrap(TCGA_SURV)
dB25, dBcl, okB = bootstrap(CGGA_SURV)


def ci(x):
    return float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))


# 逐时点配对自举（新口径）——统计 CI 下界 > 0 的时点数
rng = np.random.default_rng(20261004)
pt_diff = []
for _ in range(B):
    idx = rng.integers(0, n, n)
    df = v.iloc[idx]
    if df['ev'].sum() < 8:
        continue
    s = auc_set(df, CGGA_SURV)
    pt_diff.append(s['lp4z'][1] - s['lp25z'][1])
pt_diff = np.array(pt_diff)
lo = np.percentile(pt_diff, 2.5, axis=0)
hi = np.percentile(pt_diff, 97.5, axis=0)
n_excl = int((lo > 0).sum())

rng2 = np.random.default_rng(20261004)
ptA = []
for _ in range(B):
    idx = rng2.integers(0, n, n)
    df = v.iloc[idx]
    if df['ev'].sum() < 8:
        continue
    s = auc_set(df, TCGA_SURV)
    ptA.append(s['lp4z'][1] - s['lp25z'][1])
ptA = np.array(ptA)
loA = np.percentile(ptA, 2.5, axis=0)
n_exclA = int((loA > 0).sum())

summary = {
    'times_months': TIMES.tolist(),
    'A_ipcw_discovery': {
        'mean_MRG4': A['lp4z'][0], 'mean_ModelA': A['lp25z'][0], 'mean_Clinical': A['lpclinz'][0],
        'unw_MRG4': unweighted(A['lp4z'][1]), 'unw_ModelA': unweighted(A['lp25z'][1]),
        'auc_MRG4': A['lp4z'][1].tolist(), 'auc_ModelA': A['lp25z'][1].tolist(),
        'auc_Clinical': A['lpclinz'][1].tolist(),
        'd_vs_ModelA': float(np.mean(dA25)), 'd_vs_ModelA_ci': ci(dA25),
        'd_vs_Clinical': float(np.mean(dAcl)), 'd_vs_Clinical_ci': ci(dAcl),
        'n_valid_boot': okA, 'n_timepoints_ci_gt0': n_exclA},
    'B_ipcw_evaluation': {
        'mean_MRG4': Bs['lp4z'][0], 'mean_ModelA': Bs['lp25z'][0], 'mean_Clinical': Bs['lpclinz'][0],
        'unw_MRG4': unweighted(Bs['lp4z'][1]), 'unw_ModelA': unweighted(Bs['lp25z'][1]),
        'auc_MRG4': Bs['lp4z'][1].tolist(), 'auc_ModelA': Bs['lp25z'][1].tolist(),
        'auc_Clinical': Bs['lpclinz'][1].tolist(),
        'd_vs_ModelA': float(np.mean(dB25)), 'd_vs_ModelA_ci': ci(dB25),
        'd_vs_Clinical': float(np.mean(dBcl)), 'd_vs_Clinical_ci': ci(dBcl),
        'n_valid_boot': okB, 'n_timepoints_ci_gt0': n_excl},
    'pointwise_lo_gt0_A': int(n_exclA), 'pointwise_lo_gt0_B': int(n_excl),
    'n_validation_patients': n,
}
json.dump(summary, open('auc_ipcw_both.json', 'w'), indent=2)

# 覆盖 Table_S23：列保留原命名，另加两套均值口径列头说明
tab = pd.DataFrame({
    'time_months': TIMES,
    'AUC_MRG4': Bs['lp4z'][1],
    'AUC_ModelA': Bs['lp25z'][1],
    'AUC_ModelA_ipcw_discovery': A['lp25z'][1],
    'AUC_MRG4_ipcw_discovery': A['lp4z'][1],
    'Diff_MRG4_minus_ModelA': Bs['lp4z'][1] - Bs['lp25z'][1],
    'Boot_CI_low': lo,
    'Boot_CI_high': hi,
})
tab.to_csv('Table_S23_timedep_AUC_paired_bootstrap.csv', index=False, encoding='utf-8-sig')

print('\n--- 两套口径对照 ---')
for tag, d in (('A 发现队列删失分布', summary['A_ipcw_discovery']),
               ('B 评估队列删失分布', summary['B_ipcw_evaluation'])):
    print('%s: MRG-4 %.3f | Model A %.3f | ΔvsModelA %+.3f (%.3f-%.3f) | ΔvsClin %+.3f (%.3f-%.3f) | 逐点CI>0 %d/%d | 自举有效 %d'
          % (tag, d['mean_MRG4'], d['mean_ModelA'], d['d_vs_ModelA'], d['d_vs_ModelA_ci'][0],
             d['d_vs_ModelA_ci'][1], d['d_vs_Clinical'], d['d_vs_Clinical_ci'][0],
             d['d_vs_Clinical_ci'][1], d['n_timepoints_ci_gt0'], len(TIMES), d['n_valid_boot']))
print('\n已写 auc_ipcw_both.json 与 Table_S23_timedep_AUC_paired_bootstrap.csv')
