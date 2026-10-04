# -*- coding: utf-8 -*-
"""v5.0 图形模块：6 张主图 + 3 张补充图。

与 make_composite_figs.py 的关系
-------------------------------
数据管线**复用**原脚本（import 时其模块级 data prep 会执行），
只重做「分组与版式」：9 张旧图拆成 21 张主图面板 + 4 张补充面板再重新组合。

版式规则（本次全部版面缺陷的根因都在这里被堵住）
-----------------------------------------------
R1  任何注释文字一律经 annot() 放入坐标轴四角之一，并带不透明白底——
    保证不压曲线、不被裁。四角都放不下的信息改写进标题或图注。
R2  面板字母不进标题，改由 fig.text 放在坐标轴外左上，避免长标题互撞。
R3  分组柱的数值标签按柱序交错上下偏移，避免相邻标签互撞。
R4  面板导出为**独立矢量 PDF**（save_panel），供 Illustrator 按坐标拼版；
    同时输出整图 PNG(600 dpi) + PDF。
R5  Fig8B 的免疫差异值改为读 Table_S26c（此前是写死的旧数，与 Table 5 不符）。
R6  Fig5A 的 tdAUC 改为读 Table_S23（此前读 nri_idi_auc.json 的旧数组，
    Model A 均值印成 0.735）。
R7  Fig6D 色标标签由 "mean log1p expression" 改为 "mean raw UMI count"
    （矩阵 sl040_expr.csv 实测为整数原始计数）。
"""
import os, json, io, csv
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import Bbox
from matplotlib.patches import Patch
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from scipy.stats import norm, spearmanr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

import make_composite_figs as M          # 复用其模块级数据准备

OUTD = '/tmp/mitoxy/figs_v5'
PDIR = os.path.join(OUTD, 'panels')
os.makedirs(PDIR, exist_ok=True)
CB, CR, CG, CO, CK = M.CB, M.CR, M.CG, M.CO, M.CK
D, COH, EX, G25, MRG4, GD, GGA = M.D, M.COH, M.EX, M.G25, M.MRG4, M.GD, M.GGA
NULL, obs4, obs25 = M.NULL, M.obs4, M.obs25
S26C = pd.read_csv('Table_S26c_immune_features_by_group.csv', encoding='utf-8-sig')
S23 = pd.read_csv('Table_S23_timedep_AUC_paired_bootstrap.csv', encoding='utf-8-sig')
LOG = []

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.8,
                     'xtick.direction': 'out', 'ytick.direction': 'out',
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'figure.dpi': 300})

# ------------------------------------------------------------------ 版式工具件
CORNERS = {'nw': (0.030, 0.965, 'left', 'top'), 'ne': (0.970, 0.965, 'right', 'top'),
           'sw': (0.030, 0.035, 'left', 'bottom'), 'se': (0.970, 0.035, 'right', 'bottom')}


def annot(ax, s, loc='nw', fs=6.8, col='#222', box=True):
    """R1：注释只进四角，且带不透明白底，不可能压住曲线或被裁。"""
    x, y, ha, va = CORNERS[loc]
    ax.text(x, y, s, transform=ax.transAxes, ha=ha, va=va, fontsize=fs, color=col,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.88, pad=1.8) if box else None,
            zorder=6)


def cleantop(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def vlabel(ax, x, y, s, fs=6.2, col='#222', va='bottom', box=False):
    """柱顶数值标签。

    v6.0 关键改动：**默认不再加不透明白底**（v5.0 默认加）。
    v5.0 给每个数值标签垫了一块白底，本意是防相邻标签互压，实际效果却是
    标签把底下的柱体、曲线、散点挖掉一块——这正是 sir 看到的「遮挡」。
    正确做法不是拿白底盖住数据，而是把每个标签放到**它自己那根柱子顶端之上**，
    于是白底不再需要。`va` 放开，便于同一面板内把相邻两个标签错到上下两侧。
    """
    kw = dict(ha='center', va=va, fontsize=fs, color=col, zorder=6)
    if box:
        kw['bbox'] = dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.1)
    ax.text(x, y, s, **kw)


def annot_below(ax, s, fs=6.6, y=-0.30, ha='left', x=0.0, col='#222'):
    """把注释放到坐标轴**外面**的下方，彻底不参与绘图区。

    v6.0：annot() 的四角方案仍会把白底盖在数据上（PCA 散点图尤其明显）。
    凡属"方法学说明"性质的文字，一律走这里，放到轴下方。
    """
    ax.text(x, y, s, transform=ax.transAxes, ha=ha, va='top',
            fontsize=fs, color=col, clip_on=False)


def savefig(fig, name, tight=True):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUTD, '%s.%s' % (name, ext)), dpi=600,
                    bbox_inches='tight' if tight else None, pad_inches=0.06)
    plt.close(fig)
    LOG.append('  [整图] %s' % name)


def save_panel(axs, name):
    """R4：把一张面板导出为**整幅画布尺寸**的独立矢量 PDF/PNG，供 Illustrator 拼版。

    v5.0 关键设计：导出时隐藏其余面板、但**不裁剪**画布——
    于是每一张面板 PDF 都是同一页尺寸（如 7.4 × 7.4 in），
    在 AI 里只需把每张都放在原点、100% 缩放，拼出来的版面与 matplotlib 的
    整图逐像素一致，无需换算坐标。双轴面板（twiny/twinx）必须一并传入，
    否则会被当成"其余面板"隐藏掉。
    """
    axs = [a for a in np.atleast_1d(axs).ravel() if a is not None]
    fig = axs[0].get_figure()
    keep = {id(a) for a in axs}
    fig.patch.set_alpha(0.0)
    others = [(a, a.get_visible()) for a in fig.axes if id(a) not in keep]
    texts = [(t, t.get_visible()) for t in fig.texts]
    for a, _ in others:
        a.set_visible(False)
    for t, _ in texts:
        t.set_visible(False)
    kept_patch = [(a, a.patch.get_alpha()) for a in axs]
    for a, _ in kept_patch:
        a.patch.set_alpha(0.0)
    fig.canvas.draw()
    for ext in ('pdf', 'png'):
        # transparent=True 是关键：每张面板都是整幅画布尺寸，若带不透明白底，
        # 在 Illustrator 里叠放时上层页面会把下层全部盖掉（实测只有最后放进的
        # 那张可见）。透明底才能让 9 张面板叠加成完整版面。
        fig.savefig(os.path.join(PDIR, '%s.%s' % (name, ext)), dpi=600,
                    bbox_inches=None, pad_inches=0, transparent=True)
    for a, v in others:
        a.set_visible(v)
    for t, v in texts:
        t.set_visible(v)
    for a, v in kept_patch:
        a.patch.set_alpha(v)
    LOG.append('  [面板] %s' % name)


def km(ax, df, tcol, ecol, gcol, ttl, colors=(CB, CR), labels=('Low risk', 'High risk'),
       leg=None, note=None, note_loc='sw'):
    for g, col, lab in zip(('低风险', '高风险'), colors, labels):
        s = df[df[gcol] == g]
        if len(s) == 0:
            s = df[df[gcol] == lab]
        k = KaplanMeierFitter().fit(s[tcol], s[ecol])
        ax.step(k.survival_function_.index, k.survival_function_.iloc[:, 0], where='post',
                color=col, lw=1.45, label='%s (n = %d)' % (lab, len(s)))
        ci = k.confidence_interval_
        ax.fill_between(ci.index, ci.iloc[:, 0], ci.iloc[:, 1], color=col, alpha=0.12, lw=0)
    r = multivariate_logrank_test(df[tcol], df[gcol], df[ecol])
    ax.set_ylim(0, 1.02); ax.set_xlabel('Time (months)', fontsize=7.8)
    ax.set_ylabel('Overall survival', fontsize=7.8)
    title(ax, ttl)
    if leg:
        ax.legend(frameon=False, fontsize=6.8, loc=leg)
    if note is None:
        note = 'log-rank $P$ = %.1e' % r.p_value
    annot(ax, note, note_loc, fs=6.8)
    ax.grid(alpha=0.15, lw=0.5)
    cleantop(ax)
    return r.p_value


def forest(axT, axP, axN, rows, xlim, xticks, xticklabels, xlabel, fs_lab=7.2, fs_val=6.6):
    """森林图三段式；y 轴全部由文本承担，故横向不溢出、纵向不互撞。"""
    ys, y = [], 0
    for r in rows:
        y -= 1.35 if r[0] == 'H' else 1.0; ys.append(y)
    ymin, ymax = min(ys) - 1.0, max(ys) + 1.6
    for a in (axT, axP, axN):
        a.set_ylim(ymin, ymax); a.set_yticks([])
    for a in (axT, axN):
        for s_ in ('top', 'right', 'left', 'bottom'):
            a.spines[s_].set_visible(False)
        a.set_xticks([]); a.tick_params(bottom=False, left=False)
    axP.spines[['top', 'right']].set_visible(False)
    axP.axvline(1, color='#999', ls='--', lw=0.8); axP.set_xscale('log'); axP.set_xlim(*xlim)
    tr = axT.get_yaxis_transform()
    for r, yy in zip(rows, ys):
        if r[0] == 'H':
            axT.text(1.0, yy, r[1], ha='right', va='center', fontsize=fs_lab + 0.2,
                     fontweight='bold', color='#222', transform=tr); continue
        _, lab, hr, lo, hi, p = r
        axT.text(0.985, yy, lab, ha='right', va='center', fontsize=fs_lab, transform=tr)
        col = CR if (lo > 1 or hi < 1) else CG
        axP.plot([lo, hi], [yy, yy], color=col, lw=1.3, solid_capstyle='round')
        axP.plot([hr], [yy], 's', color=col, ms=4.0)
        axN.text(0.01, yy, '%.2f (%.2f–%.2f)' % (hr, lo, hi), ha='left', va='center',
                 fontsize=fs_val, transform=axN.get_yaxis_transform())
        axN.text(0.995, yy, ('%.1e' % p) if p < 0.001 else ('%.3f' % p), ha='right', va='center',
                 fontsize=fs_val, transform=axN.get_yaxis_transform())
    axT.text(1.0, ymax - 0.10, 'Covariate', ha='right', va='center', fontsize=fs_lab,
             style='italic', color='#444', transform=tr)
    axN.text(0.01, ymax - 0.10, 'HR (95% CI)', ha='left', va='center', fontsize=fs_lab,
             style='italic', color='#444', transform=axN.get_yaxis_transform())
    axN.text(0.995, ymax - 0.10, '$P$', ha='right', va='center', fontsize=fs_lab,
             style='italic', color='#444', transform=axN.get_yaxis_transform())
    axP.set_xticks(xticks); axP.set_xticklabels(xticklabels, fontsize=6.8)
    axP.set_xlabel(xlabel, fontsize=7.6); axP.grid(axis='x', alpha=0.12, lw=0.5)


def maxticks(ax, n=5, x=True, y=True):
    """只给*未显式设定*刻度的那条轴加 MaxNLocator。

    v5.0 修正：旧写法对已 set_xticks/set_xticklabels 的轴（panel_2C、panel_3D）
    覆盖 locator，导致 FixedFormatter 的少量标签被循环贴到更多刻度上，
    于是出现 "IDH-mutant(n=410) IDH-wildtype(n=92)" 这类标签互撞。
    """
    from matplotlib.ticker import MaxNLocator
    if x:
        ax.xaxis.set_major_locator(MaxNLocator(n))
    if y:
        ax.yaxis.set_major_locator(MaxNLocator(n))


def title(ax, s, fs=8.2, pad=4):
    ax.set_title(s, fontsize=fs, pad=pad)


def panellabel(fig, ax, L, dx=None, dy=0.006):
    """R2：面板字母置于坐标轴外左上，与标题解耦。

    规则（v5.0 修正版）：
      · 偏移量按坐标轴宽度缩放——固定 -0.085 在窄面板（1×3 图）会跑到隔壁面板头顶；
      · 只有当标题左缘确实让出了位置时才把字母放在标题左侧；
        若标题很宽（左缘已越过坐标轴左缘），就把字母抬到标题上方。
      · 包围盒必须用 fig.transFigure 转成图幅比例；用 dpi_scale_trans 会得到英寸，
        从而把 fig.text 放到图外，bbox_inches='tight' 会把画布撑到几万像素高。
    """
    fig.canvas.draw()
    bb = ax.get_position()
    if dx is None:
        dx = -max(0.030, min(0.085, 0.14 * bb.width))
    try:
        r = fig.canvas.get_renderer()
        tb = ax.title.get_window_extent(r).transformed(fig.transFigure.inverted())
    except Exception:
        tb = None
    fs = plt.rcParams['font.size']
    default_x = max(bb.x0 + dx, 0.004)
    if tb is None or not ax.get_title() or tb.x0 - 0.012 > default_x:
        x, y = default_x, min(bb.y1 + dy, 0.995)
    else:                                   # 标题太宽 → 抬到标题上方
        x = max(bb.x0 - 0.012, 0.004)
        y = min(tb.y1 + 0.004, 0.995)
    fig.text(x, y, L, fontsize=fs + 2.4, fontweight='bold', ha='left', va='bottom')


# =====================================================================================
# 面板：发现与 Model A
# =====================================================================================
def panel_1A(fig, spec):
    # v5.0：改读沉积集的唯一 S1 源表。原先读的 Table_S1_MRG_panel_LGG.csv
    # 是同一份数据的 10 列旧副本（缺 BH q 列），既非沉积源、又在一次目录清理中
    # 被回收；凡图 1A 用到的列（Human_symbol / LGG_univariate_HR / _p）本表都有。
    s1 = pd.read_csv('Table_S1_v4.csv')
    s1 = s1.dropna(subset=['LGG_univariate_HR', 'LGG_univariate_p']).copy()
    z = norm.isf(s1['LGG_univariate_p'].clip(lower=1e-300) / 2)
    se = np.log(s1['LGG_univariate_HR']).abs() / z
    s1['lo'] = np.exp(np.log(s1['LGG_univariate_HR']) - 1.96 * se)
    s1['hi'] = np.exp(np.log(s1['LGG_univariate_HR']) + 1.96 * se)
    s1 = s1.sort_values('LGG_univariate_p').reset_index(drop=True)
    n = len(s1); yy = np.arange(n)[::-1]
    ax = fig.add_subplot(spec); xr = (0.35, 3.2)
    for i, r in s1.iterrows():
        y0 = yy[i]; p = r['LGG_univariate_p']; col = CR if p < 0.05 else CG
        lo, hi = r['lo'], r['hi']
        mx = max(lo, 1 / hi); mn = min(hi, 1 / lo)
        if (mn > xr[1]) or (mx < xr[0]):
            ax.plot([max(lo, xr[0]), min(hi, xr[1])], [y0, y0], color=col, lw=0.8, zorder=2)
            if lo < xr[0]:
                ax.annotate('', xy=(xr[0], y0), xytext=(xr[0] * 1.09, y0),
                            arrowprops=dict(arrowstyle='-|>', color=col, lw=0.8))
            if hi > xr[1]:
                ax.annotate('', xy=(xr[1], y0), xytext=(xr[1] / 1.09, y0),
                            arrowprops=dict(arrowstyle='-|>', color=col, lw=0.8))
        else:
            ax.plot([lo, hi], [y0, y0], color=col, lw=0.8, zorder=2)
        ax.plot([np.clip(r['LGG_univariate_HR'], *xr)], [y0], 'o', ms=2.1, color=col, zorder=3)
    ax.axvline(1.0, color='#333', lw=0.8, ls='--', zorder=1)
    ax.set_yticks(yy); ax.set_yticklabels(list(s1['Human_symbol']), fontsize=4.5)
    ax.set_ylim(-1, n); ax.set_xscale('log'); ax.set_xlim(*xr)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0])
    ax.set_xticklabels(['0.4', '0.6', '0.8', '1.0', '1.5', '2.0', '3.0'], fontsize=7)
    ax.set_xlabel('Hazard ratio (95% CI, log scale)', fontsize=7.6)
    cleantop(ax); ax.grid(axis='x', alpha=0.12, lw=0.5)
    title(ax, 'Univariate Cox regression, %d MRGs' % n)
    return [ax]


def panel_1B(fig, spec):
    dz = D.sort_values('risk').reset_index(drop=True)
    sub = spec.subgridspec(2, 1, height_ratios=[1, 1], hspace=0.10)
    b1 = fig.add_subplot(sub[0]); b2 = fig.add_subplot(sub[1], sharex=b1)
    b1.fill_between(np.arange(len(dz)), 0, dz['risk'], where=dz['hi'] == 0, color=CB, lw=0, alpha=0.9)
    b1.fill_between(np.arange(len(dz)), 0, dz['risk'], where=dz['hi'] == 1, color=CR, lw=0, alpha=0.9)
    b1.axhline(0, color='#333', lw=0.5); b1.set_ylabel('Risk score', fontsize=7.0)
    b1.set_xticks([]); b1.set_xlim(0, len(dz))
    b2.scatter(np.arange(len(dz)), dz['t'], s=1.5, c=np.where(dz['e'] == 1, CR, CB), linewidths=0)
    b2.set_ylabel('Survival (mo)', fontsize=7.0); b2.set_yscale('log')
    b2.set_yticks([12, 36, 120, 240]); b2.set_yticklabels(['12', '36', '120', '240'], fontsize=6.2)
    b2.set_ylim(5, 320); b2.set_xlabel('Patients ordered by Model A risk score', fontsize=7.0)
    # R1/R2 修正：图例移到坐标轴**外**的左上，不再压住散点；标题承担“存活/死亡”说明
    b2.scatter([], [], s=7, c=CR, label='Dead'); b2.scatter([], [], s=7, c=CB, label='Alive')
    b2.legend(frameon=False, fontsize=6.2, loc='lower left', bbox_to_anchor=(0.0, 1.02),
              ncol=2, borderaxespad=0.0, handletextpad=0.4, columnspacing=1.2)
    for a in (b1, b2):
        cleantop(a)
    title(b1, 'Risk score and survival status', pad=16)
    return [b1, b2]


def panel_1C(fig, spec):
    c = fig.add_subplot(spec)
    order = D.sort_values('risk')['upk'].tolist()
    X5 = EX.loc[order][G25].values.T
    im = c.imshow(X5, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2, interpolation='nearest')
    c.set_yticks(np.arange(len(G25)))
    c.set_yticklabels(G25, fontsize=4.6)          # 基因名按 25 行均分，逐行可读
    c.set_xlabel('Patients ordered by Model A risk score', fontsize=7.2)
    cb = fig.colorbar(im, ax=c, fraction=0.028, pad=0.014)
    cb.set_label('Expression ($z$)', fontsize=6.4); cb.ax.tick_params(labelsize=6.0)
    title(c, 'Expression heat map, 25 Model A genes')
    # v6.2 关键修正（2026-10-04）：色条轴**必须一并返回**，否则 save_panel 会把它
    # 当成「其余面板」隐藏掉——导出的面板 PDF 里热图右侧干干净净，既没有色条也
    # 没有 "Expression (z)" 标签。而 savefig(整图) 走的是另一条路径，画的是
    # fig.axes 全体，所以**正文里那张 PDF/PNG 有色条、Illustrator 拼版件没有**，
    # 两者对不上。panel_6B 早就返回了 [b, cb.ax]，这里当时漏了。
    # 是 compose_ai_fallback.py 的墨迹宽度比对（3709 vs 3909 px）把它抓出来的。
    return [c, cb.ax]


def panel_1D(fig, spec):
    ax = fig.add_subplot(spec)
    km(ax, D, 't', 'e', '组', 'Kaplan–Meier, Model A, TCGA-LGG', leg='upper right',
       note='HR 3.27 (2.21–4.83)\nmedian 105.2 vs 51.9 mo', note_loc='sw')
    return [ax]


# =====================================================================================
# 面板：Model A 表现与随机面板基准
# =====================================================================================
def _td_roc(score, t, e, horizon):
    kmc = KaplanMeierFitter().fit(t, 1 - np.asarray(e).astype(int))
    Gf = lambda q: np.interp(q, kmc.survival_function_.index.values,
                             kmc.survival_function_.iloc[:, 0].values, left=1.0)
    cases = (e == 1) & (t <= horizon); ctrl = t > horizon
    wc = 1.0 / Gf(t[cases]); wk = 1.0 / Gf(np.full(ctrl.sum(), horizon))
    sc, sk = score[cases], score[ctrl]
    thr = np.unique(np.r_[sc, sk])[::-1]
    tpr = np.array([wc[sc >= c].sum() for c in thr]) / wc.sum()
    fpr = np.array([wk[sk >= c].sum() for c in thr]) / wk.sum()
    return fpr, tpr, np.trapezoid(tpr, fpr)


def panel_2A(fig, spec):
    a = fig.add_subplot(spec); aucs = {}
    for yr, col in ((1, CB), (3, CO), (5, CR)):
        fpr, tpr, auc = _td_roc(D['risk'].values, D['t'].values, D['e'].values, yr * 12)
        aucs[yr] = auc
        a.plot(fpr, tpr, color=col, lw=1.5, label='%d-year (AUC %.3f)' % (yr, auc))
    a.plot([0, 1], [0, 1], color=CG, lw=0.8, ls='--')
    a.set_xlabel('1 − specificity', fontsize=7.6); a.set_ylabel('Sensitivity', fontsize=7.6)
    a.set_xlim(0, 1); a.set_ylim(0, 1.02)
    a.legend(frameon=False, fontsize=6.8, loc='lower right'); a.grid(alpha=0.15, lw=0.5)
    cleantop(a); title(a, 'Time-dependent ROC, Model A')
    LOG.append('  [2A] tdROC %.3f / %.3f / %.3f' % (aucs[1], aucs[3], aucs[5]))
    return [a]


def panel_2B(fig, spec):
    b = fig.add_subplot(spec)
    rho, pv = spearmanr(GD['risk_hat'], GD['prolif'])
    b.scatter(GD['prolif'], GD['risk_hat'], s=3.2, alpha=0.45, color=CB, lw=0)
    b.set_xlabel('Proliferation score (mean $z$)', fontsize=7.6)
    b.set_ylabel('Model A risk score', fontsize=7.6)
    # v6.1 修 2B 遮挡：v5.0 用 annot(..., 'nw') 带 88% 不透明白底，而这幅是
    # **铺满整个绘图区**的散点图，'nw' 角下边就是数据点——白底把点挖掉了一个方块
    # （字形级审计报出「散点 1 点」）。散点图没有"空白角"可用，注记只能移出绘图区：
    # 并进标题第二行。标题在轴外，永远不可能压住数据。
    b.grid(alpha=0.15, lw=0.5); cleantop(b); maxticks(b)
    title(b, 'Not a proliferation score\n'
             'Spearman $\\rho$ = %.2f; $P$ = %.1e' % (rho, pv))
    return [b]


def panel_2C(fig, spec):
    A_ = M.A_
    c = fig.add_subplot(spec)
    bp = c.boxplot([A_[A_['idhwt'] == 0]['risk'], A_[A_['idhwt'] == 1]['risk']],
                   patch_artist=True, widths=0.55, medianprops=dict(color='black'),
                   tick_labels=['IDH-mutant\n(n = %d)' % int((A_['idhwt'] == 0).sum()),
                                'IDH-wildtype\n(n = %d)' % int((A_['idhwt'] == 1).sum())])
    for p, col in zip(bp['boxes'], [CB, CR]):
        p.set_facecolor(col); p.set_alpha(0.55)
    c.set_ylabel('Model A risk score', fontsize=7.6)
    c.grid(alpha=0.15, lw=0.5); cleantop(c); maxticks(c, 4, x=False)
    c.tick_params(axis='x', labelsize=6.6)
    title(c, 'Risk score by IDH status')
    return [c]


def panel_2D(fig, spec):
    """随机面板零分布 vs 观测值。

    v5.0 修正：旧版把说明塞进超宽图例、把两条竖线标签放在 y = 20.6，
    三者挤在同一条带上互压成一团。现把纵轴上限抬到 32，
    曲线最高约 17，于是顶部 24–32 整条只归图例、21 归竖线标签，互不相犯。
    """
    a = fig.add_subplot(spec)
    xs = np.linspace(0.45, 0.90, 500)
    dens = lambda mu, sd: np.exp(-((xs - mu) ** 2) / (2 * sd ** 2)) / (sd * np.sqrt(2 * np.pi))
    k4 = NULL['k=4, matched (unpenalised)']; k25 = NULL['k=25, matched (pen 1.0)']
    k25u = NULL['k=25, unpenalised']
    h = [plt.Line2D([], [], color=CR, lw=1.3,
                    label='4-gene panels, unpenalised (n = %d); %.2f%% match/beat MRG-4'
                          % (k4['valid'], k4['frac_ge_obs'] * 100)),
         plt.Line2D([], [], color=CB, lw=1.3,
                    label='25-gene panels, $\\lambda$ = 1.0 (n = %d); %.1f%% match/beat Model A'
                          % (k25['valid'], k25['frac_ge_obs'] * 100)),
         plt.Line2D([], [], color=CB, lw=1.1, ls=':',
                    label='25-gene panels, unpenalised (n = %d); %.1f%% at/above Model A'
                          % (k25u['valid'], k25u['frac_ge_obs'] * 100))]
    a.plot(xs, dens(k4['mean'], k4['sd']), color=CR, lw=1.3)
    a.plot(xs, dens(k25['mean'], k25['sd']), color=CB, lw=1.3)
    a.plot(xs, dens(k25u['mean'], k25u['sd']), color=CB, lw=1.1, ls=':')
    # v5.0：竖虚线只画到 y = 19.5（曲线最高约 17），不再穿过顶部图例文字
    a.vlines(obs4, 0, 19.5, colors=CR, linewidths=1.8, linestyles='--')
    a.vlines(obs25, 0, 19.5, colors=CB, linewidths=1.8, linestyles='--')
    a.text(obs4 + 0.005, 21.0, 'MRG-4 %.3f' % obs4, color=CR, fontsize=6.8, va='center', ha='left',
           bbox=dict(facecolor='white', edgecolor='none', alpha=0.92, pad=1.5))
    a.text(obs25 - 0.005, 21.0, 'Model A %.3f' % obs25, color=CB, fontsize=6.8, va='center', ha='right',
           bbox=dict(facecolor='white', edgecolor='none', alpha=0.92, pad=1.5))
    a.set_xlabel('External C-index in CGGA-325', fontsize=7.6)
    a.set_ylabel('Null density (illustrative)', fontsize=7.6)
    a.set_ylim(0, 32.0); a.set_xlim(0.45, 0.90)
    a.legend(handles=h, frameon=True, facecolor='white', framealpha=0.95, edgecolor='none',
             fontsize=5.6, loc='upper left', bbox_to_anchor=(0.005, 1.005),
             borderaxespad=0.0, handlelength=1.9, labelspacing=0.30, borderpad=0.25)
    cleantop(a); a.grid(alpha=0.12, lw=0.5)
    title(a, 'Observed panels vs matched random panels')
    return [a]


# =====================================================================================
# 面板：MRG-4 的推导
# =====================================================================================
def panel_3A(fig, spec):
    P = M.PAR; cv = P['cv']
    alphas = np.array(cv['alphas']); ks = np.array(cv['k'])
    mean = np.array(cv['mean']); se = np.array(cv['se'])
    a = fig.add_subplot(spec)
    x = np.log10(alphas); m = ks <= 60
    a.errorbar(x[m], mean[m], yerr=se[m], fmt='o', ms=2.2, lw=0.8, color=CB, capsize=1.4, alpha=0.85)
    x1, xm = np.log10(cv['alpha_1se']), np.log10(cv['alpha_min'])
    a.axvline(x1, color=CR, ls='--', lw=1.0)
    a.axvline(xm, color=CK, ls=':', lw=1.0)
    # v6.0：两条竖线的说明文字从绘图区**移到轴下方**做成图例。
    # 起因（sir 2026-10-04 指出 3A 有遮挡）：v5.0 把 λ_1SE / λ_min 两行字放在
    # y≈0.85/0.805 处，白色底框把那一带的误差棒整段挖空，而且红色虚线正好
    # 从 "λ_1SE (26)" 中间穿过。竖线本身必须贯穿全高才能表意，所以文字只能让位。
    from matplotlib.lines import Line2D
    hh = [Line2D([], [], color=CR, ls='--', lw=1.0,
                 label='$\\lambda_{1\\mathrm{SE}}$ = %.5f (%d genes)' % (cv['alpha_1se'], cv['k_1se'])),
          Line2D([], [], color=CK, ls=':', lw=1.0, label='$\\lambda_{\\min}$ (43 genes)')]
    a.set_xlabel('log$_{10}(\\lambda)$ (=$\\log_{10}\\alpha$)', fontsize=7.4)
    a.set_ylabel('10-fold CV C-index', fontsize=7.6)
    a.set_ylim(0.48, 0.86); a.grid(alpha=0.15, lw=0.5); cleantop(a)
    a.legend(handles=hh, frameon=False, fontsize=6.3, loc='upper center',
             bbox_to_anchor=(0.5, -0.375), ncol=1, borderaxespad=0.0,
             handlelength=1.8, columnspacing=1.4, labelspacing=0.25)
    axt = a.twiny(); axt.set_xlim(a.get_xlim())
    axt.set_xticks([np.log10(alphas[np.argmin(np.abs(ks - k))]) for k in [5, 10, 20, 30]])
    axt.set_xticklabels(['5', '10', '20', '30'], fontsize=6.8)
    axt.set_xlabel('Number of genes', fontsize=7.2)
    title(a, 'LASSO-Cox cross-validation', pad=8)
    return [a, axt]          # axt 必须一并返回，否则导出面板时被隐藏


def panel_3B(fig, spec):
    b = fig.add_subplot(spec); P = M.PAR
    st = sorted(P['stability'].items(), key=lambda kv: -kv[1])[:12][::-1]
    cols = [CR if v >= 0.75 else CG for _, v in st]
    b.barh([g for g, _ in st], [v * 100 for _, v in st], color=cols, alpha=0.88, height=0.66)
    b.axvline(75, color='black', ls='--', lw=0.9)
    annot(b, '$\\pi$ = 0.75', 'se', fs=6.6, box=True)
    b.set_xlabel('Selection frequency (%)', fontsize=7.6); b.set_xlim(0, 100)
    b.tick_params(labelsize=6.8); b.grid(axis='x', alpha=0.15, lw=0.5); cleantop(b)
    title(b, 'Stability selection')
    return [b]


def panel_3C(fig, spec):
    """阈值敏感性 + Meinshausen–Bühlmann 界。

    v5.0 修正两处：
      · q 由 26 改为 **38**——正文 §2.4 与 Figure 4 图注都声明 q = 38
        （200 次重抽中最大选入基因数），旧图却按 q = 26 画，图与正文两个数；
      · 右轴标签改成短名，且右轴线不再带独立 ylabel 挤进右侧面板的 ylabel。
    """
    c = fig.add_subplot(spec)
    PS = M.json.load(open('rev_extra_results.json'))['E_pi_sensitivity']
    pv_ = [r['pi'] for r in PS]; kv_ = [r['k'] for r in PS]
    c.bar(range(3), kv_, width=0.52, color=[CG, CO, CR], alpha=0.9)
    for i, v in enumerate(kv_):
        vlabel(c, i, v + 0.25, '%d' % v, fs=7.2)
    c.set_xticks(range(3)); c.set_xticklabels(['0.70', '0.75', '0.80'], fontsize=7.4)
    c.set_xlabel('Stability threshold $\\pi$', fontsize=7.6)
    c.set_ylabel('Genes retained', fontsize=7.6)
    c.set_ylim(0, 16.5); c.grid(axis='y', alpha=0.15, lw=0.5)
    for s_ in ('top',):
        c.spines[s_].set_visible(False)
    c2 = c.twinx()
    q, p = 38, 89                                   # 与正文 §2.4 / 图注一致
    mb = [q ** 2 / ((2 * t - 1) * p) for t in pv_]
    c2.plot(range(3), mb, 'o--', color='#444444', ms=4, lw=1.1)
    for i, v in enumerate(mb):
        c2.text(i, v * 1.045, '%.1f' % v, fontsize=6.3, color='#444444', ha='center', va='bottom')
    c2.set_ylim(0, 52)
    c2.tick_params(axis='y', labelsize=6.4)
    c2.set_ylabel('E[V] permitted', fontsize=6.6)
    c2.spines['top'].set_visible(False)
    # v6.0：说明文字由 'se' 改到 'ne'，并且**挂在 twinx 轴 c2 上**。
    # ① 'se' 原位正好压在 π = 0.80 那根红柱上（白底把柱体挖掉一块）；
    #    红柱只有 4 高，而 ylim 到 16.5，右上方那一片是空的，标在那里不压任何东西。
    # ② c2 是在 c 之后加入的坐标轴，会整体后绘；挂在 c 上会被 c2 的虚线/数值压住。
    annot(c2, 'the four MRG-4\nat $\\pi$ = 0.80', 'ne', fs=6.5, col=CR)
    title(c, 'Threshold sensitivity')
    return [c, c2]           # c2 同上


def panel_3D(fig, spec):
    """留一基因：外部队列 C-index。

    v6.0 修三处：
      ① **数据错**：浅灰柱过去用 fz[i]（每个基因被剔除后的值）当柱高，
         但图注写的是"浅灰柱 = 全文模型在 discovery-SD 口径下的值 0.731"。
         现在浅灰柱统一取 G7['full_cgga_frozen_C'] = 0.7314，蓝方块才是剔基因值。
      ② 数值标签过去贴在"点的上方 0.005"，而红点在灰柱内部，于是标签连白底
         一起压在柱体上；PDGFA 组红蓝两标签还互压。改为一律抬高到柱顶之上。
      ③ "Model A 0.664" 原来放在 'se'，位置正好落在最后一组柱子里；改为写进图例。
    """
    G7 = json.load(open('rev_extra_results.json'))['G_leave_one_out']
    LOO = json.load(open('rev_round1_add.json'))['loo']['cohort_z']
    b = fig.add_subplot(spec)
    genes = [r['dropped'] for r in LOO['rows']]
    cz = [r['C'] for r in LOO['rows']]
    fz = [r['cgga_frozen_C'] for r in G7['rows']]
    full_z = LOO['full']                      # 0.7688，cohort-z 全文模型
    full_d = G7['full_cgga_frozen_C']         # 0.7314，discovery-SD 全文模型
    xp = np.arange(len(genes)); w = 0.44
    b.bar(xp - w / 2, [full_z] * len(genes), w, color='#dddddd',
          label='full MRG-4, cohort $z$')
    b.bar(xp + w / 2, [full_d] * len(genes), w, color='#c6c6c6', edgecolor='none',
          label='full MRG-4, disc. SD')
    b.plot(xp - w / 2, cz, 'o', color=CR, ms=5.5, label='dropped, cohort $z$')
    b.plot(xp + w / 2, fz, 's', color=CB, ms=4.5, label='dropped, disc. SD')
    # v6.1：两个数值标签统一抬到**全场最高柱顶之上**，并且柱宽 0.34 → 0.44。
    #   ① 旧版蓝标签放在 max(full_d, fz)+0.008，而 full_d (0.7314) 恒高于剔基因值，
    #      于是标签的实际落点就在浅灰柱面上——审计器报「柱/面积 9%」。
    #   ② 更深一层：6 pt 的 "0.680" 字形盒宽 0.414 个类别单位，而 w = 0.34 的柱
    #      只有 0.34 宽，标签**比它自己那根柱还宽**，左右各溢出 0.04 撞进邻柱。
    #      柱宽加到 0.44（半宽 0.22 > 字形半宽 0.21）后两处一起解决。
    #   ③ 但两个标签一左一右摆在同一高度上，中心只隔 0.44 个类别单位，
    #      而每个字形盒本身就有 0.41 宽——渲染出来是 "0.7470.680" 粘成一片。
    #      类别间距是 1.0，两枚标签就吃掉 0.82，横向无论如何挤不下，
    #      所以**改为纵向错开**：蓝标签再抬高 0.019 个单位（≈9 pt，
    #      大于字形盒高 4.5 pt），上下两排完全分离，横向便不必再让。
    ytop = max([full_z, full_d] + list(cz) + list(fz)) + 0.006
    for i in range(len(genes)):
        vlabel(b, i - w / 2, ytop, '%.3f' % cz[i], fs=5.8, col=CR)
        vlabel(b, i + w / 2, ytop + 0.019, '%.3f' % fz[i], fs=5.8, col=CB)
    b.axhline(obs25, color=CK, ls='--', lw=1.0, label='Model A (%.3f)' % obs25)
    b.set_xticks(xp); b.set_xticklabels(['−' + g for g in genes], fontsize=6.6)
    b.set_ylim(0.60, 0.86); b.set_ylabel('External C-index', fontsize=7.4)
    # R1/R2：图例移到坐标轴下方（两列），不再与柱顶数值标签同层。
    # v6.0：标签改短，避免图例比面板还宽横向溢出；位置上抬到 -0.17，免得第三行落到画布下缘外。
    b.legend(frameon=False, fontsize=5.9, loc='upper center', bbox_to_anchor=(0.5, -0.17),
             ncol=2, borderaxespad=0.0, handlelength=1.4, columnspacing=1.0, labelspacing=0.25)
    b.grid(axis='y', alpha=0.15, lw=0.5); cleantop(b)
    title(b, 'Leave-one-gene-out, CGGA-325')
    return [b]


# =====================================================================================
# 面板：外部验证
# =====================================================================================
def panel_4A(fig, spec):
    a = fig.add_subplot(spec)
    km(a, GGA, 'time', 'ev', 'grp4_c',
       'CGGA-325 WHO II–III\n(external, n = %d)' % len(GGA),
       colors=(CB, CR), labels=('Low MRG-4', 'High MRG-4'), leg='upper right')
    return [a]


def panel_4B(fig, spec):
    M4 = json.load(open('mrg4_results.json')); FD = json.load(open('fig13_data.json'))['bars']
    b = fig.add_subplot(spec)
    labels = ['25-gene\nModel A', 'LASSO-26', 'MRG-4\n4 genes']
    series = [('TCGA apparent', [FD['25-gene\noriginal'][0], FD['LASSO\n(26 genes)'][0], M4['tcga']['apparent_C']], CB),
              ('TCGA 10-fold CV', [FD['25-gene\noriginal'][1], FD['LASSO\n(26 genes)'][1], M4['tcga']['cv_C']], '#7fb3d5'),
              ('CGGA external', [0.664, 0.726, obs4], CR)]
    xpos = np.arange(3); w = 0.30
    # v6.0 修 4B（两处）：
    #  ① 数值标签过去统一排在"该组最高柱之上"的一条基线上。第一组的 "0.664"
    #     （所属红柱最矮）被抬到 0.716，恰好落进同组浅蓝柱（柱顶 0.796）的柱体里，
    #     白底把柱体挖掉一块；第二组的 "0.726" 也压在浅蓝柱顶。
    #  ② 改为"各标签只放在自己那根柱顶之上"后，第三组又出新问题：0.794 / 0.796
    #     几乎等高，两个标签并排横向相碰（渲染出来是 "0.7940.796"）。
    #  最终方案：第一、三序列贴各自柱顶；**中间那一枚统一抬到同组最高柱之上 0.06**，
    #  于是它必然与左右两个标签纵向错开，两个问题一起解决。
    # v6.1：w 0.26 → 0.30、标签 5.5 → 5.4 pt。这是同一类"标签比柱子还宽"的毛病：
    #  "0.664" 在 5.5 pt 下约 15.6 pt 宽，而 w=0.26 的柱只有 15.0 pt，
    #  标签左右各溢出一点，蹭进邻柱（审计器报「柱/面积 2%」）。柱宽加大即可根治。
    for i, (nm, vals, col) in enumerate(series):
        b.bar(xpos + (i - 1) * w, vals, w, label=nm, color=col, alpha=0.92)
        for j, v in enumerate(vals):
            if i == 1:
                y = max(series[0][1][j], v, series[2][1][j]) + 0.060
            else:
                y = v + 0.014
            b.text(xpos[j] + (i - 1) * w, y, '%.3f' % v, ha='center', va='bottom',
                   fontsize=5.4, color=col if i != 1 else '#3d7ea6', zorder=6)
    b.set_xticks(xpos); b.set_xticklabels(labels, fontsize=7.0)
    b.set_ylabel('C-index', fontsize=7.6); b.set_ylim(0, 1.06)
    b.axhline(0.5, color=CG, ls='--', lw=0.8)
    # v6.0：图例由绘图区内左上移到轴下方——原地正好占着抬高后的中间标签要用的高度。
    b.legend(frameon=False, fontsize=5.6, loc='upper center', bbox_to_anchor=(0.5, -0.22),
             ncol=2, borderaxespad=0.0, handlelength=1.3, columnspacing=1.2, labelspacing=0.25)
    b.grid(axis='y', alpha=0.15, lw=0.5); cleantop(b)
    title(b, 'Discrimination by panel size', pad=26)
    return [b]


def panel_4C(fig, spec):
    """两指数的校正模型森林图，四块并列。

    v5.0：旧稿把 Model A 的校正森林图整块删掉（reviewer 指出信息缺失），
    这里把 Model A 与 MRG-4 放在同一坐标系里、同一 x 轴、同一行样式，
    读者可直接比较"同一个协变量集下两个指数各自的 HR"。
    数值全部来自沉积对象：TCGA 取 tables_v42.json 的 T4/T8，CGGA 取 Table_S14。
    另：旧版"Model A alone, CGGA-325"一行曾误用校正值 lpA_z（1.34），此处不再设"alone"行。
    """
    T4 = json.load(open('tables_v42.json'))['T4']
    T8 = json.load(open('tables_v42.json'))['T8']
    S14 = pd.read_csv('Table_S14_CGGA_adjusted_models.csv', encoding='utf-8-sig')
    mA = S14[S14['Model'] == 'Model A'].set_index('Variable')
    m4 = S14[S14['Model'] == 'MRG-4 cohort-z'].set_index('Variable')

    def t8(blk, key):
        r = T8[blk][[x['var'] for x in T8[blk]].index(key)]
        return ('R', None, r['HR'], r['lo'], r['hi'], r['P'])

    def t4(i, key):
        r = T4[i][[x['var'] for x in T4[i]].index(key)]
        return r['HR'], r['lo'], r['hi'], r['P']

    rows = [('H', 'MRG-4 — TCGA-LGG, multivariable (n = 502)')]
    for lab, k in [('Risk score (per SD)', 'MRG4_score_z'), ('Age (per year)', 'age'), ('Male sex', 'male'),
                   ('IDH-wildtype', 'idhwt'), ('WHO grade 3', 'g3')]:
        r = t8('M2', k); rows.append(('R', lab, r[2], r[3], r[4], r[5]))
    rows.append(('H', 'MRG-4 — CGGA-325 WHO II–III (n = 172)'))
    for lab, k in [('Risk score (per SD)', 'lp4_z'), ('Age (per year)', 'age'), ('Male sex', 'male'),
                   ('IDH-wildtype', 'idhwt'), ('WHO grade III', 'g3')]:
        rows.append(('R', lab, float(m4.loc[k, 'HR']), float(m4.loc[k, 'CI_low']),
                     float(m4.loc[k, 'CI_high']), float(m4.loc[k, 'P'])))
    rows.append(('H', 'Model A — TCGA-LGG, multivariable (n = 502)'))
    for lab, k in [('Risk score (per SD)', 'risk_g_z'), ('Age (per year)', 'age'), ('Male sex', 'male'),
                   ('IDH-wildtype', 'idhwt'), ('WHO grade 3', 'g3')]:
        hr, lo, hi, pv = t4('M2', k); rows.append(('R', lab, hr, lo, hi, pv))
    rows.append(('H', 'Model A — CGGA-325 WHO II–III (n = 172)'))
    for lab, k in [('Risk score (per SD)', 'lpA_z'), ('Age (per year)', 'age'), ('Male sex', 'male'),
                   ('IDH-wildtype', 'idhwt'), ('WHO grade III', 'g3')]:
        rows.append(('R', lab, float(mA.loc[k, 'HR']), float(mA.loc[k, 'CI_low']),
                     float(mA.loc[k, 'CI_high']), float(mA.loc[k, 'P'])))

    gs = spec.subgridspec(1, 3, width_ratios=[1.32, 1.00, 0.78], wspace=0.05)
    aT = fig.add_subplot(gs[0]); aP = fig.add_subplot(gs[1]); aN = fig.add_subplot(gs[2])
    forest(aT, aP, aN, rows, (0.25, 20), [0.3, 0.5, 1, 2, 4, 8], ['0.3', '0.5', '1', '2', '4', '8'],
           'Hazard ratio (log scale)')
    title(aT, '', pad=0)
    return [aT, aP, aN]


# =====================================================================================
# 面板：头对头、校准与决策分析
# =====================================================================================
def panel_5A(fig, spec):
    """时间依赖 AUC，CGGA-325。

    R6：改读 Table_S23（重算所得），不再用 nri_idi_auc.json 里 Model A 均值 0.735 的旧数组。
    v5.0 再改一处口径：IPCW 的删失分布必须取自**被评估队列自身**（CGGA，随访 130.2 月），
    而不是发现队列（TCGA，随访 27.8 月）——后者会把权重算错。更正后
    事件加权均值 0.845 / 0.731（旧口径 0.819 / 0.711）。
    """
    a = fig.add_subplot(spec)
    t = S23['time_months'].values / 12.0
    a4 = S23['AUC_MRG4'].values; a25 = S23['AUC_ModelA'].values
    a.plot(t, a4, '-o', color=CR, ms=3.0, lw=1.5, label='MRG-4 (event-weighted mean 0.845)')
    a.plot(t, a25, '-s', color=CB, ms=3.0, lw=1.5, label='Model A (event-weighted mean 0.731)')
    a.axhline(0.5, color=CG, ls=':', lw=0.9)
    for xi, y4, y25 in zip(t, a4, a25):
        a.plot([xi, xi], [y25, y4], color='#bbbbbb', lw=0.6, zorder=1)
    # 口径说明属方法学，写进图注而不是压在曲线区（旧版此处与图例互撞）
    a.set_xlabel('Time (years)', fontsize=7.6); a.set_ylabel('Time-dependent AUC', fontsize=7.6)
    a.set_ylim(0.45, 1.0); a.set_xticks([1, 2, 4, 6, 8, 10]); a.set_xlim(0.5, 10.5)
    # v6.0 修 5A：图例原本在 'lower left' 且锚在 0.055，正好被 0.5 处的横线
    # （axhline）和网格线从字里穿过（sir 指出的 5a 遮挡）。曲线最高只到 0.89，
    # 顶部 0.90–1.00 那一条带子是空的，图例移到那里即可完全不与任何线相交。
    a.legend(frameon=False, fontsize=6.0, loc='upper left', borderaxespad=0.35,
             labelspacing=0.25, handlelength=1.6)
    a.grid(alpha=0.20, lw=0.5); cleantop(a)
    title(a, 'Discrimination in CGGA-325')
    return [a]


def panel_5B(fig, spec):
    ADD = json.load(open('rev_round1_add.json'))
    L = {e['label']: e for e in ADD['dca']}
    b = fig.add_subplot(spec)
    m4, m25 = L['MRG-4'], L['Model A (25 genes)']
    mcl = L['Clinical model (age, sex, IDH, grade)']
    pt = np.array(m4['pts'], dtype=float); prev = m4['deaths3y'] / m4['n']
    b.plot(pt, m4['nb'], '-o', color=CR, ms=2.6, lw=1.6, label='MRG-4')
    b.plot(pt, m25['nb'], '-s', color=CB, ms=2.6, lw=1.6, label='Model A (25 genes)')
    b.plot(pt, mcl['nb'], '-^', color='#b8860b', ms=2.8, lw=1.7, label='Clinical (age, sex, IDH, grade)')
    b.plot(pt, prev - (1 - prev) * (pt / (1 - pt)), '--', color='#777777', lw=0.9, label='Treat all')
    b.axhline(0.0, color='#444444', ls=':', lw=0.9, label='Treat none')
    # v6.1：把"7% 透明红"改成**预先与白底混好的浅红实色**。
    # 起因：`axvspan(..., color=CR, alpha=0.07)` 在这一版 matplotlib 里
    # 于本图的绘制顺序下会渲染成一整块**饱和正红**（实测像素 (192,57,43,255)，
    # 即 CR 原色、alpha 完全不生效），整块红底把图例文字压在里面。
    # 同一张图上后加的、位置不同的 span 却正常显示为浅色——说明不是 alpha 写法问题，
    # 而是这块 patch 的合成结果不可靠。既然如此就不再依赖 alpha 叠加：
    # 直接给一个已经混好的颜色，任何合成路径下结果都一样。
    # `#f2d7d5` = CR(192,57,43) 按 20% 与白底相混，视觉上仍是"淡红高亮带"。
    b.axvspan(0.24, 0.36, color='#f2d7d5', lw=0, zorder=0)
    b.set_xlabel('Threshold probability (3-year mortality)', fontsize=7.6)
    b.set_ylabel('Net benefit', fontsize=7.6)
    b.set_ylim(-0.06, 0.40); b.set_xlim(0.20, 0.80)      # 抬高上限，给图例留空
    b.legend(frameon=False, fontsize=5.9, loc='upper right', ncol=2, columnspacing=1.0,
             handlelength=1.5, labelspacing=0.30, borderaxespad=0.2)
    b.grid(alpha=0.20, lw=0.5); cleantop(b)
    title(b, 'Decision curve analysis, 3 years')
    return [b]


def panel_5C(fig, spec):
    R1 = json.load(open('rev_round1_add.json'))
    c = fig.add_subplot(spec)
    cal = pd.DataFrame(R1['calibration']['tertiles'])
    c.plot([0, 1], [0, 1], ls='--', lw=0.9, color=CG)
    c.plot(cal['predicted_3y'], cal['observed_3y'], 'o-', color=CB, ms=5.0, lw=1.5)
    for _, r in cal.iterrows():
        # v5.0：右侧三分位的标签改为向左展开，避免 Q1 那行顶到坐标轴右边界
        right = r['predicted_3y'] > 0.70
        c.annotate('%s (n=%d)' % (r['tertile'], r['n']), (r['predicted_3y'], r['observed_3y']),
                   textcoords='offset points', xytext=(-8 if right else 8, -3), fontsize=6.4,
                   ha='right' if right else 'left',
                   bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.2), zorder=6)
    c.set_xlabel('Predicted 3-year OS', fontsize=7.6)
    c.set_ylabel('Observed 3-year OS (KM)', fontsize=7.6)
    c.set_xlim(0, 1); c.set_ylim(0, 1)
    annot(c, 'calibration slope %.3f\n(95%% CI %.3f–%.3f)' % (
        R1['calibration']['slope'], R1['calibration']['lo'], R1['calibration']['hi']), 'nw', fs=6.6)
    c.grid(alpha=0.15, lw=0.5); cleantop(c)
    title(c, 'Calibration in CGGA-325, 3-year OS')
    return [c]


def panel_5D(fig, spec):
    R1 = json.load(open('rev_round1_add.json'))
    d = fig.add_subplot(spec); B = R1['brier']
    d.plot(B['times'], B['mrg4'], 'o-', color=CR, ms=3.4, lw=1.5, label='MRG-4')
    d.plot(B['times'], B['g25'], 's--', color=CB, ms=3.4, lw=1.5, label='Model A (25 genes)')
    d.set_xlabel('Time (months)', fontsize=7.6); d.set_ylabel('Time-dependent Brier score', fontsize=7.6)
    d.legend(frameon=False, fontsize=6.6, loc='upper left')
    # R1：integrated Brier 移出绘图区，改为标题第二行，不再压曲线
    title(d, 'Prediction error (lower is better)\nintegrated Brier score %.3f vs %.3f' % (B['ibs4'], B['ibs25']))
    d.grid(alpha=0.15, lw=0.5); cleantop(d); maxticks(d, 4)
    return [d]


# =====================================================================================
# 面板：单细胞定位
# =====================================================================================
def _dot(ax, Mm, P, order, labels, fs=7.0, rot=0):
    df = pd.read_csv('sl040_expr.csv', index_col=0)
    cols = ['BMP1', 'KIF15', 'TRAF3', 'PDGFA', 'MKI67']
    cmap = LinearSegmentedColormap.from_list('gyr', ['#F2F4F7', '#FDE2E2', '#F5A3A3', '#D7301F', '#7A1004'])
    for i, ct in enumerate(order):
        for j, g in enumerate(cols):
            r = 3 + 100 * np.sqrt(P.iloc[i, j] / 100.0)
            ax.scatter(j, i, s=r, c=[cmap(min(Mm.iloc[i, j] / 4.0, 1.0))],
                       edgecolors='#4A4A4A', linewidths=0.3, zorder=3)
    ax.set_xticks(range(len(cols)))
    # v6.0：45° 斜排在窄面板上仍会互压（B 面板只有约 1.4 in 绘图宽，5 个基因名
    # 斜着排的实际水平footprint超过刻度间距，"TRAF3" 直接糊进 "PDGFA"）。
    # rot >= 60 时改为竖直、居中悬于刻度下方，水平方向只剩一个字符宽，杜绝相撞。
    if rot >= 60:
        # 注意不要用 rotation_mode='anchor'：它先按未旋转状态对齐、再绕锚点旋转，
        # 结果是各标签顶端高低不齐像扇面。不加 rotation_mode 时 matplotlib 按
        # **旋转后的包围盒**对齐，va='top' 就能让 5 个标签顶端严格齐平。
        ax.set_xticklabels(cols, fontsize=fs, style='italic', rotation=rot,
                           ha='center', va='top')
    else:
        ax.set_xticklabels(cols, fontsize=fs, style='italic', rotation=rot,
                           ha='right' if rot else 'center')
    ax.set_yticks(range(len(order))); ax.set_yticklabels(labels, fontsize=7.0)
    ax.invert_yaxis(); ax.set_ylim(len(order) - 0.4, -0.6); ax.set_xlim(-0.6, len(cols) - 0.4)
    for j in range(1, len(cols)):
        ax.axvline(j - 0.5, color='#D9D9D9', lw=0.5, zorder=0)
    for s_ in ('top', 'right'):
        ax.spines[s_].set_visible(False)
    return cmap


def _ctname(c):
    """细胞类型名缩写：只用于 y 轴刻度标签，避免长名越出画布左缘。

    v5.0：'Oligodendrocytes (n=10,895)' 在 7.8 in 画布上曾向左溢出 0.54 in。
    """
    return {'Oligodendrocytes': 'Oligodendro.', 'Macrophages': 'Macrophage'}.get(c, c)


def panel_6A(fig, spec):
    df = pd.read_csv('sl040_expr.csv', index_col=0)
    cols = ['BMP1', 'KIF15', 'TRAF3', 'PDGFA', 'MKI67']
    ORDER = [c for c in ['Tumor', 'Astrocytes', 'OPCs', 'Oligodendrocytes', 'Endothelial',
                         'Pericytes', 'VSMC', 'Fibroblasts', 'Microglia', 'MgTAM', 'MoTAM',
                         'Macrophages', 'Monocytes', 'NK, T-cells', 'B-cells']
             if c in set(df['CellType'])]
    mean = df.groupby('CellType')[cols].mean().loc[ORDER]
    pct = df.groupby('CellType')[cols].apply(lambda d: (d > 0).mean() * 100).loc[ORDER]
    cnt = df['CellType'].value_counts().loc[ORDER]
    a = fig.add_subplot(spec)
    _dot(a, mean, pct, ORDER, ['%s (n=%s)' % (_ctname(c), format(cnt[c], ',')) for c in ORDER])
    title(a, 'All cells, by cell type')
    return [a]


def panel_6B(fig, spec):
    df = pd.read_csv('sl040_expr.csv', index_col=0)
    cols = ['BMP1', 'KIF15', 'TRAF3', 'PDGFA', 'MKI67']
    t_ = df[df['CellType'] == 'Tumor']
    SO = [s for s in ['NPC', 'OPC', 'AC', 'MES'] if s in set(t_['NeftelClass'].dropna())]
    sm = t_.groupby('NeftelClass')[cols].mean().reindex(SO)
    sp = t_.groupby('NeftelClass')[cols].apply(lambda d: (d > 0).mean() * 100).reindex(SO)
    sc = t_['NeftelClass'].value_counts().reindex(SO)
    b = fig.add_subplot(spec)
    cmap = _dot(b, sm, sp, SO, ['%s (n=%s)' % (s, format(sc[s], ',')) for s in SO], fs=6.4, rot=90)
    h = [plt.scatter([], [], s=3 + 100 * np.sqrt(pp / 100.0), c='#C9C9C9', edgecolors='#4A4A4A',
                     linewidths=0.3, label='%d%%' % pp) for pp in (10, 25, 50, 75)]
    b.legend(handles=h, title='% expressing', frameon=False, fontsize=6.0, title_fontsize=6.2,
             loc='upper left', bbox_to_anchor=(1.14, 1.0), labelspacing=0.75, handletextpad=1.0,
             borderpad=0.2)
    sm_ = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 4))
    cb = b.get_figure().colorbar(sm_, ax=b, fraction=0.042, pad=0.20)
    cb.set_label('mean raw UMI count', fontsize=6.4)     # R7：与矩阵实测一致
    cb.ax.tick_params(labelsize=5.9)
    title(b, 'Malignant cells, by state')
    return [b, cb.ax]     # 色条轴必须一并返回，否则导出 B 面板时被 save_panel 隐藏
