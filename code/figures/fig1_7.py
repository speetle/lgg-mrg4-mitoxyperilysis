# -*- coding: utf-8 -*-
"""
Fig 1–7 — regenerated for 学生1 LGG 交付 v3.1

原则（与"数据红线"一致）：
  · 凡稿件/补充表已沉积的统计量 —— 直接作图，不重算（Fig 1 取自 Table_S1；Fig 7 取自 Table 5）。
  · 凡仅依赖已沉积风险评分与临床数据的量 —— 本地复算（Fig 2 KM、Fig 3 时间依赖 ROC、Fig 4 三联图）。
  · Fig 5/6 用到表达矩阵（cBioPortal RNA-Seq V2 median all-sample Z-scores）。
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import concordance_index
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from scipy.stats import norm

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.linewidth': 0.9,
                     'xtick.direction': 'out', 'ytick.direction': 'out', 'figure.dpi': 300})
CB = '#1f5fa8'; CR = '#c0392b'; CG = '#888888'; CO = '#e08a1e'

COH = pd.read_csv('lgg_g1_dataset.csv', encoding='utf-8-sig')
D = COH.dropna(subset=['risk', 't', 'e']).copy()
D['hi'] = (D['组'] == '高风险').astype(int)
EX = pd.read_csv('panel_expr_tcga.csv', index_col=0)
IDX = [i for i in EX.index if i in set(COH['upk'])]
EX = EX.loc[IDX]
G25 = ['KIF15', 'CHAF1A', 'PDGFA', 'PPFIA1', 'WDR19', 'PPM1M', 'ASF1B', 'LOXL3', 'RALGPS1',
       'GALNT2', 'MXRA8', 'OSGIN2', 'BMP1', 'HSPA4L', 'TLR8', 'POLM', 'EPS8', 'PALLD',
       'MCM5', 'LIG1', 'RUSC2', 'KIF24', 'SLC7A1', 'CALCRL', 'TAF12']
NRISK = 8


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(f'{name}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('OK ', name)


# ============ Figure 1 — univariate Cox forest (from deposited Table_S1) ============
s1 = pd.read_csv('Table_S1_MRG_panel_LGG.csv')
s1 = s1.dropna(subset=['LGG_univariate_HR', 'LGG_univariate_p']).copy()
# 95% CI reconstructed from the reported HR and two-sided Wald P (no value invented)
z = norm.isf(s1['LGG_univariate_p'].clip(lower=1e-300) / 2)
se = np.log(s1['LGG_univariate_HR']).abs() / z
s1['lo'] = np.exp(np.log(s1['LGG_univariate_HR']) - 1.96 * se)
s1['hi'] = np.exp(np.log(s1['LGG_univariate_HR']) + 1.96 * se)
s1 = s1.sort_values('LGG_univariate_p', ascending=True).reset_index(drop=True)
n = len(s1)
fig, ax = plt.subplots(figsize=(5.4, 0.155 * n + 1.5))
yy = np.arange(n)[::-1]


def _gapmask(lo, hi, xr):
    mx = np.maximum(lo, 1 / hi); mn = np.minimum(hi, 1 / lo)
    m = (mn > xr[1]) | (mx < xr[0])
    return m


xr = (0.35, 3.2)
for i, r in s1.iterrows():
    y = yy[i]
    p = r['LGG_univariate_p']; hr = r['LGG_univariate_HR']
    col = CR if p < 0.05 else CG
    lo, hi = r['lo'], r['hi']
    if _gapmask(lo, hi, xr):   # CI 超出轴范围 → 用箭头截断绘制
        lox, hix = max(lo, xr[0]), min(hi, xr[1])
        ax.plot([lox, hix], [y, y], color=col, lw=0.9, solid_capstyle='butt', zorder=2)
        if lo < xr[0]:
            ax.annotate('', xy=(xr[0], y), xytext=(xr[0] * 1.09, y),
                        arrowprops=dict(arrowstyle='-|>', color=col, lw=0.9))
        if hi > xr[1]:
            ax.annotate('', xy=(xr[1], y), xytext=(xr[1] / 1.09, y),
                        arrowprops=dict(arrowstyle='-|>', color=col, lw=0.9))
    else:
        ax.plot([lo, hi], [y, y], color=col, lw=0.9, solid_capstyle='butt', zorder=2)
    ax.plot([np.clip(hr, *xr)], [y], 'o', ms=2.5, color=col, zorder=3)
ax.axvline(1.0, color='#333333', lw=0.8, ls='--', zorder=1)
ax.set_yticks(yy); ax.set_yticklabels(list(s1['Human_symbol']), fontsize=5.0)
ax.set_ylim(-1, n)
ax.set_xscale('log')
ax.set_xlim(*xr)
ax.set_xticks([0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0])
ax.set_xticklabels(['0.4', '0.6', '0.8', '1.0', '1.5', '2.0', '3.0'], fontsize=7.5)
ax.set_xlabel('Hazard ratio per SD (95% CI, log scale)', fontsize=8.5)
ax.set_title('Univariate Cox regression of %d mitoxyperilysis-related genes' % n, fontsize=9.4)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax.grid(axis='x', alpha=0.12, lw=0.5)
ax.text(0.99, 0.004, 'red: $P$ < 0.05    grey: $P$ $\\geq$ 0.05', transform=ax.transAxes,
        ha='right', va='bottom', fontsize=6.6, color=CG)
save(fig, 'Fig1_单因素Cox森林图')

# ============ Figure 2 — Kaplan–Meier, Model A ============
fig, (ax, axr) = plt.subplots(2, 1, figsize=(5.4, 4.4), height_ratios=[4, 1], sharex=True)
for g, col, lab in ((0, CB, 'Low risk'), (1, CR, 'High risk')):
    s = D[D['hi'] == g]
    k = KaplanMeierFitter().fit(s['t'], s['e'])
    ax.step(k.survival_function_.index, k.survival_function_.iloc[:, 0], where='post',
            color=col, lw=1.7, label=f'{lab} (n = {len(s)})')
    ci = k.confidence_interval_
    ax.fill_between(ci.index, ci.iloc[:, 0], ci.iloc[:, 1], color=col, alpha=0.13, lw=0)
tt = np.arange(0, 181, 12)
for g, col in ((0, CB), (1, CR)):
    k = KaplanMeierFitter().fit(D.loc[D['hi'] == g, 't'], D.loc[D['hi'] == g, 'e'])
    sf = k.survival_function_.reindex(k.survival_function_.index.union(tt)).ffill().reindex(tt)
    ax.plot(tt, sf.iloc[:, 0], '+', ms=5, color=col, mew=0.9)
r = multivariate_logrank_test(D['t'], D['hi'], D['e'])
c = CoxPHFitter().fit(D[['t', 'e', 'hi']], 't', 'e').summary.loc['hi']
ax.text(0.97, 0.93,
        'log-rank $P$ = %.2e\nHR %.2f (%.2f–%.2f)\nmedian 105.2 vs 51.9 mo' %
        (r.p_value, c['exp(coef)'], c['exp(coef) lower 95%'], c['exp(coef) upper 95%']),
        transform=ax.transAxes, ha='right', va='top', fontsize=7.8)
ax.set_ylabel('Overall survival', fontsize=9)
ax.set_ylim(0, 1.03); ax.set_xlim(0, 180)
ax.legend(frameon=False, fontsize=8.2, loc='lower left')
ax.set_title('Model A (25 genes), TCGA-LGG, median risk-score split', fontsize=9.4)
ax.grid(alpha=0.15, lw=0.5)
for row, (g, col) in enumerate(((0, CB), (1, CR))):
    s = D[D['hi'] == g]
    k = KaplanMeierFitter().fit(s['t'], s['e'])
    yv = k.survival_function_.reindex(k.survival_function_.index.union(tt)).ffill().reindex(tt).iloc[:, 0]
    yc = 1 - row
    axr.plot(tt, np.full(len(tt), yc), color=col, lw=0)
    for x, v in zip(tt, yv):
        axr.text(x, yc, '%.0f' % round(v * len(s)), ha='center', va='center',
                 fontsize=6.2, color=col)
axr.set_ylim(-0.6, 1.6); axr.set_yticks([]); axr.set_xlabel('Time (months)', fontsize=9)
axr.set_ylabel('No. at risk', fontsize=7.5, rotation=0, ha='right', va='center')
for s_ in ('top', 'right', 'left'):
    axr.spines[s_].set_visible(False)
axr.grid(alpha=0.10, lw=0.4)
save(fig, 'Fig2_KM曲线_ModelA')

# ============ Figure 3 — time-dependent ROC (1/3/5 years) ============
def km_cens(t, e):
    km = KaplanMeierFitter().fit(t, 1 - np.asarray(e).astype(int))
    xs = km.survival_function_.index.values
    ys = km.survival_function_.iloc[:, 0].values
    return lambda q: np.interp(q, xs, ys, left=1.0)


def td_roc(score, t, e, horizon):
    G = km_cens(t, e)
    cases = (e == 1) & (t <= horizon)
    ctrl = t > horizon
    if cases.sum() == 0 or ctrl.sum() == 0:
        return None, None, np.nan
    wc = 1.0 / G(t[cases]); wk = 1.0 / G(np.full(ctrl.sum(), horizon))
    sc, sk = score[cases], score[ctrl]
    thr = np.unique(np.r_[sc, sk])[::-1]
    tpr = np.array([wc[sc >= c].sum() for c in thr]) / wc.sum()
    fpr = np.array([wk[sk >= c].sum() for c in thr]) / wk.sum()
    auc = np.trapezoid(tpr, fpr)
    return fpr, tpr, auc


fig, ax = plt.subplots(figsize=(4.9, 4.5))
res = {}
for yr, col in ((1, CB), (3, CO), (5, CR)):
    fpr, tpr, auc = td_roc(D['risk'].values, D['t'].values, D['e'].values, yr * 12)
    res[yr] = auc
    ax.plot(fpr, tpr, color=col, lw=1.7, label='%d-year (AUC %.3f)' % (yr, auc))
ax.plot([0, 1], [0, 1], color=CG, lw=0.9, ls='--')
ax.set_xlabel('1 − specificity', fontsize=9); ax.set_ylabel('Sensitivity', fontsize=9)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
ax.set_title('Time-dependent ROC, Model A (TCGA-LGG)', fontsize=9.4)
ax.legend(frameon=False, fontsize=8.2, loc='lower right')
ax.grid(alpha=0.15, lw=0.5)
for s_ in ('top', 'right'):
    ax.spines[s_].set_visible(False)
save(fig, 'Fig3_时间依赖ROC')
print('   时间依赖 AUC:', {k: round(v, 4) for k, v in res.items()}, '(稿件 0.907 / 0.861 / 0.801)')

# ============ Figure 4 — risk-score triple plot ============
dz = D.sort_values('risk').reset_index(drop=True)
X4 = EX.loc[dz['upk']][G25].values.T
fig = plt.figure(figsize=(6.4, 5.4))
gs = GridSpec(3, 1, height_ratios=[1, 1, 2.9], hspace=0.10)
a1 = fig.add_subplot(gs[0])
a1.fill_between(np.arange(len(dz)), 0, dz['risk'], where=dz['hi'] == 0, color=CB, lw=0, alpha=0.9)
a1.fill_between(np.arange(len(dz)), 0, dz['risk'], where=dz['hi'] == 1, color=CR, lw=0, alpha=0.9)
a1.axhline(0, color='#333', lw=0.6)
a1.set_ylabel('Risk score', fontsize=8.4); a1.set_xlim(0, len(dz))
a1.set_xticks([]); a1.set_title('Model A risk-score triple plot (TCGA-LGG, n = %d)' % len(dz),
                                fontsize=9.4)
a2 = fig.add_subplot(gs[1], sharex=a1)
a2.scatter(np.arange(len(dz)), dz['t'], s=2.4,
           c=np.where(dz['e'] == 1, CR, CB), linewidths=0)
a2.set_ylabel('Survival (mo)', fontsize=8.4); a2.set_yscale('log')
a2.set_yticks([12, 36, 120, 240]); a2.set_yticklabels(['12', '36', '120', '240'], fontsize=7.5)
a2.set_ylim(5, 320)
a2.set_xticks([])
a2.scatter([], [], s=8, c=CR, label='Dead'); a2.scatter([], [], s=8, c=CB, label='Alive')
a2.legend(frameon=False, fontsize=7.4, loc='lower left', ncol=2)
a3 = fig.add_subplot(gs[2], sharex=a1)
im = a3.imshow(X4, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2, interpolation='nearest')
a3.set_yticks(np.arange(len(G25))); a3.set_yticklabels(G25, fontsize=6.2)
a3.set_xlabel('Patients ordered by risk score', fontsize=8.6)
for a in (a1, a2, a3):
    for s_ in ('top', 'right'):
        a.spines[s_].set_visible(False)
cb = fig.colorbar(im, ax=a3, fraction=0.030, pad=0.012)
cb.set_label('Expression (Z)', fontsize=7.6); cb.ax.tick_params(labelsize=6.8)
save(fig, 'Fig4_风险评分三联图')

# ============ Figure 5 — heat map of the 25 model genes ============
order = D.sort_values('risk')['upk'].tolist()
X5 = EX.loc[order][G25].values.T
fig, ax = plt.subplots(figsize=(6.6, 2.9))
im = ax.imshow(X5, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2, interpolation='nearest')
ax.set_yticks(np.arange(len(G25))); ax.set_yticklabels(G25, fontsize=6.4)
ax.set_xticks([0, 100, 200, 300, 400, 509])
ax.set_xticklabels(['0', '100', '200', '300', '400', '510'], fontsize=7.5)
ax.set_xlabel('Patients (ordered by Model A risk score)', fontsize=8.6)
ax.set_title('Z-score heat map of the 25 Model A genes', fontsize=9.4)
cb = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.012)
cb.set_label('Expression (Z)', fontsize=7.6); cb.ax.tick_params(labelsize=6.8)
save(fig, 'Fig5_25基因热图')

# ============ Figure 6 — PCA / k-means subtypes ============
X = EX[G25].values
p = PCA(n_components=2).fit(X); Z = p.transform(X)
km = KMeans(n_clusters=2, n_init=50, random_state=42).fit(X)
lab = km.labels_
if np.bincount(lab)[0] < np.bincount(lab)[1]:
    lab = 1 - lab
sil = silhouette_score(X, lab)
sz = np.bincount(lab)
dd = COH.set_index('upk').loc[EX.index].copy(); dd['cl'] = lab
lr = multivariate_logrank_test(dd['t'], dd['cl'], dd['e'])
fig, ax = plt.subplots(figsize=(5.0, 4.6))
for c, col, nm in ((0, CB, 'C1'), (1, CR, 'C2')):
    m = lab == c
    ax.scatter(Z[m, 0], Z[m, 1], s=6, color=col, alpha=0.75, linewidths=0,
               label='%s (n = %d, %d deaths)' % (nm, sz[c], int(dd.loc[dd['cl'] == c, 'e'].sum())))
ax.set_xlabel('PC1 (%.1f%%)' % (p.explained_variance_ratio_[0] * 100), fontsize=9)
ax.set_ylabel('PC2 (%.1f%%)' % (p.explained_variance_ratio_[1] * 100), fontsize=9)
ax.set_title('k-means subtypes of the 25 model genes (k = 2)', fontsize=9.4)
ax.text(0.03, 0.03, 'silhouette = %.3f\nbetween-subtype log-rank $P$ = %.2e' % (sil, lr.p_value),
        transform=ax.transAxes, fontsize=7.8, va='bottom')
ax.legend(frameon=False, fontsize=7.6, loc='upper right')
ax.grid(alpha=0.12, lw=0.5)
for s_ in ('top', 'right'):
    ax.spines[s_].set_visible(False)
save(fig, 'Fig6_分子分型PCA')
print('   k-means 复算: PC1 %.1f%% PC2 %.1f%% | sizes %s | silhouette %.3f | log-rank %.2e'
      % (p.explained_variance_ratio_[0] * 100, p.explained_variance_ratio_[1] * 100,
         sz.tolist(), sil, lr.p_value))

# ============ Figure 7 — immune features (from deposited Table 5) ============
# 数值取自稿件 Table 5（= LGG示范结果报告 表6），逐条照抄，未重算
IMM = [('CD8+ T cells', 0.311, 2.31e-12), ('CD4+ T cells', 0.278, 4.86e-12),
       ('Regulatory T cells', 0.184, 1.97e-5), ('NK cells', 0.259, 5.53e-11),
       ('B cells', 0.167, 1.76e-16), ('M1 macrophages', 0.166, 2.23e-8),
       ('M2 macrophages', 0.153, 0.0284), ('Dendritic cells', 0.303, 2.38e-15),
       ('Neutrophils', 0.238, 4.86e-10), ('Immune checkpoints', 0.248, 8.21e-15),
       ('Interferon signature', 0.241, 4.98e-10)]
IMM = sorted(IMM, key=lambda x: x[1])
fig, ax = plt.subplots(figsize=(5.6, 3.9))
yy = np.arange(len(IMM))
ax.barh(yy, [x[1] for x in IMM], height=0.62, color=CR, alpha=0.88)
for i, (nm, v, pv) in enumerate(IMM):
    st = '***' if pv < 1e-3 else '**' if pv < 0.01 else '*'
    ax.text(v + 0.006, i, '%.3f  %s' % (v, st), va='center', fontsize=7.6)
ax.set_yticks(yy); ax.set_yticklabels([x[0] for x in IMM], fontsize=8.2)
ax.set_xlim(0, max(x[1] for x in IMM) * 1.32)
ax.set_xlabel('log$_2$ fold change (high- vs low-risk)', fontsize=8.8)
ax.set_title('Immune feature scores, high- vs low-risk (n = 255 per group)', fontsize=9.4)
ax.grid(axis='x', alpha=0.15, lw=0.5); ax.set_axisbelow(True)
for s_ in ('top', 'right'):
    ax.spines[s_].set_visible(False)
fig.text(0.5, -0.015, 'all 11 features elevated; *** $P$ < 0.001, ** $P$ < 0.01, * $P$ < 0.05',
         ha='center', va='top', fontsize=6.6, color=CG)
save(fig, 'Fig7_免疫浸润')
print('\n全部完成。')
