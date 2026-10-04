# -*- coding: utf-8 -*-
"""v5.0 拼版：把面板组成 6 张主图 + 3 张补充图，并导出每张面板的矢量文件。

用法：python make_figs_v5.py            # 全部
      python make_figs_v5.py F1 F3      # 只做指定图
输出：
  /tmp/mitoxy/figs_v5/           整图 PNG(600 dpi) + PDF
  /tmp/mitoxy/figs_v5/panels/    逐面板 PDF + PNG（供 Illustrator 拼版用）
"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import fig_panels_v5 as P

FIG = {}
# v5.3（sir 2026-10-04）：**图号按正文首次引用顺序重排**。
#   v5.2 把原 FigS1 / FigS2 并进正文后（→ Figure 7 / Figure 8），编号顺序与首次引用顺序
#   脱钩，实测为 1, 2, 4, 3, 5, 7, 8, 6：Figure 4（外部验证）在 §3.2 的
#   "(Figure 1D, 2A) … (Figure 4C; Table 3)" 处就被引到，早于 §3.3 的 Figure 3；
#   Figure 7/8 在 §3.6/§3.8 被引，早于 §3.10 的 Figure 6（单细胞）。多数期刊要求
#   图号按首次引用递增，故整体重排。映射（旧→新）：
#       1→1  2→2  3→4  4→3  5→5  6→8  7→6  8→7
#   即 3 与 4 对调，6→8 / 7→6 / 8→7 三循环。
#   键 = **新图号**；内容绘制函数不动，因此 BUILDERS 里 F3↔F4、F6↔F7↔F8 交叉指回去。
PLAN = {
    'F1': ('Figure1_发现与ModelA', dict(figsize=(7.4, 7.7))),
    'F2': ('Figure2_ModelA与随机面板基准', dict(figsize=(7.4, 5.5))),
    'F3': ('Figure3_MRG4外部验证', dict(figsize=(7.4, 7.4))),
    'F4': ('Figure4_MRG4推导与稳定性', dict(figsize=(7.4, 5.4))),
    'F5': ('Figure5_头对头与决策分析', dict(figsize=(7.4, 5.5))),
    'F6': ('Figure6_分层生存', dict(figsize=(7.4, 2.6))),
    # v6.1：6.4 → 7.6 in 高。原版面下面板 C 只有 2.30 in 装 10 行（含 Fisher P 表头
    # 共 10.9 个数据单位），1 单位 = 15.2 pt；而 6 pt 的字形墨迹高 4.6 pt，
    # 柱高 0.38 单位只有 5.8 pt——字形比它自己那根柱还高，蓝柱数值标签必然
    # 蹭进相邻红柱 0.9 pt（审计器报 20%）。加高面板 C 到 ~2.97 in 后墨迹完全入柱。
    'F7': ('Figure7_分子相关性', dict(figsize=(7.4, 7.6))),
    'F8': ('Figure8_单细胞定位', dict(figsize=(7.8, 4.9))),
    # v5.0：FigS3 撤下 —— build_S3 只是 P.panel_2C 的放大副本，与 Figure 2C 重复。
}

# ---------------------------------------------------------------------------
# 图幅尺寸表 —— 按**内容**维护，与图号解耦。
# 键 = Figure 名里 "_" 之后的内容串；build_Fx() 正是按内容取尺寸的。
# ⚠️ 不要把它改成 "按图号查"：v5.3 的重排（3↔4、6→8/7→6/8→7）会立刻让
#    「图号」与「内容」错位，而错的只是画布大小——图照样画得出来，
#    只有拿 layout.json 复核 figsize 才看得见。下面的断言就是拦这个的。
# ---------------------------------------------------------------------------
SIZE = {
    '发现与ModelA':        (7.4, 7.7),
    'ModelA与随机面板基准':   (7.4, 5.5),
    'MRG4推导与稳定性':      (7.4, 5.4),
    'MRG4外部验证':         (7.4, 7.4),
    '头对头与决策分析':       (7.4, 5.5),
    '单细胞定位':           (7.8, 4.9),
    '分层生存':            (7.4, 2.6),
    '分子相关性':           (7.4, 7.6),
}
for _k, (_nm, _opt) in PLAN.items():
    _content = _nm.split('_', 1)[1]
    if _content not in SIZE:
        raise SystemExit('!! PLAN[%s] 的内容 %r 不在 SIZE 表里' % (_k, _content))
    if tuple(_opt['figsize']) != SIZE[_content]:
        raise SystemExit('!! PLAN[%s]=%s 的 figsize %s 与 SIZE[%r]=%s 不一致'
                         % (_k, _nm, tuple(_opt['figsize']), _content, SIZE[_content]))


def build_F1():
    fig = plt.figure(figsize=SIZE['发现与ModelA'])
    gs = GridSpec(3, 2, width_ratios=[0.86, 1.30], height_ratios=[0.80, 0.80, 1.35],
                  wspace=0.36, hspace=0.62)
    axs = {}
    axs['A'] = P.panel_1A(fig, gs[:, 0])
    axs['B'] = P.panel_1B(fig, gs[0, 1])
    axs['C'] = P.panel_1C(fig, gs[1, 1])
    axs['D'] = P.panel_1D(fig, gs[2, 1])
    return fig, axs


def build_F2():
    fig = plt.figure(figsize=SIZE['ModelA与随机面板基准'])
    gs = GridSpec(2, 2, wspace=0.36, hspace=0.52)
    axs = {'A': P.panel_2A(fig, gs[0, 0]), 'B': P.panel_2B(fig, gs[0, 1]),
           'C': P.panel_2C(fig, gs[1, 0]), 'D': P.panel_2D(fig, gs[1, 1])}
    return fig, axs


def build_F3():
    fig = plt.figure(figsize=SIZE['MRG4推导与稳定性'])
    # v5.0：wspace 由 0.36 放到 0.58——左下面板带 twinx 右轴，
    # 它的右轴标签原本正压在右下面板的 y 轴标签上。
    # v6.0：hspace 由 0.62 放到 0.74——3A 的 λ 说明改成轴下方图例，需要这一行的间隙。
    gs = GridSpec(2, 2, wspace=0.58, hspace=0.74)
    axs = {'A': P.panel_3A(fig, gs[0, 0]), 'B': P.panel_3B(fig, gs[0, 1]),
           'C': P.panel_3C(fig, gs[1, 0]), 'D': P.panel_3D(fig, gs[1, 1])}
    return fig, axs


def build_F4():
    fig = plt.figure(figsize=SIZE['MRG4外部验证'])
    # v5.0：C 面板并入 Model A 的校正森林图（四块共 20 行），故加高下半区。
    gs = GridSpec(2, 2, height_ratios=[0.86, 1.42], wspace=0.40, hspace=0.38)
    axs = {'A': P.panel_4A(fig, gs[0, 0]), 'B': P.panel_4B(fig, gs[0, 1]),
           'C': P.panel_4C(fig, gs[1, :])}
    return fig, axs


def build_F5():
    fig = plt.figure(figsize=SIZE['头对头与决策分析'])
    gs = GridSpec(2, 2, wspace=0.36, hspace=0.58)
    axs = {'A': P.panel_5A(fig, gs[0, 0]), 'B': P.panel_5B(fig, gs[0, 1]),
           'C': P.panel_5C(fig, gs[1, 0]), 'D': P.panel_5D(fig, gs[1, 1])}
    return fig, axs


def build_F6():
    """单细胞定位。

    v5.0：A 面板 y 轴刻度标签最长到 'Oligodendrocytes (n=10,895)'，
    默认 left=0.125 时向左溢出画布 0.54 in（正文图注里被裁成 'ocytes'）。
    这里把左边距放到 0.185，并把两个最长名做缩写（见 fig_panels_v5._ctname）。

    v6.1（sir 2026-10-04 点名 Fig6 有重叠）：wspace 0.30 → 0.46。
    实测墨迹盒 A 右沿 0.5703 × 7.80 in、B 左沿 0.5494——**A 的最后一个
    斜排刻度标签 MKI67 与 B 的 y 轴标签（NPC/OPC/AC/MES）在水平方向交了
    0.163 in**。check_panel_overlap.py 报出后按需把列间距拉开。
    """
    fig = plt.figure(figsize=SIZE['单细胞定位'])
    gs = GridSpec(1, 2, width_ratios=[1.62, 1.0], wspace=0.54)
    gs.update(left=0.185)
    axs = {'A': P.panel_6A(fig, gs[0]), 'B': P.panel_6B(fig, gs[1])}
    return fig, axs


def build_F7():
    """分层生存三张 KM，用 Model A 风险分组。

    v5.0：标题去掉 "(n = …)"——面板窄、标题长时会顶到左缘，与面板字母互压
    （曾出现 "n = 4B10"）；样本量图例里已经给了。
    v6.0：原 build_S1，图号由 Supplementary FigS1 升为正文 Figure 7。
    """
    fig = plt.figure(figsize=SIZE['分层生存'])
    gs = GridSpec(1, 3, wspace=0.50)
    # v5.0：xlabel 'Time (months)' 在 bottom=0.11 时向下溢出 0.086 in；抬到 0.185。
    gs.update(bottom=0.185)
    A_ = P.M.A_
    sub = A_[A_['idhwt'] == 0]
    sub3 = A_[A_['g3'] == 1]
    a = fig.add_subplot(gs[0])
    pa = P.km(a, sub, 't', 'e', '组', 'TCGA-LGG, IDH-mutant', leg='upper right', note='')
    P.annot(a, 'n = %d\nlog-rank $P$ = %.1e' % (len(sub), pa), 'sw', fs=6.6)
    b = fig.add_subplot(gs[1])
    pb = P.km(b, sub3, 't', 'e', '组', 'TCGA-LGG, WHO grade 3', leg='upper right', note='')
    P.annot(b, 'n = %d\nlog-rank $P$ = %.1e' % (len(sub3), pb), 'sw', fs=6.6)
    c = fig.add_subplot(gs[2])
    pc = P.km(c, P.GGA, 'time', 'ev', 'grp25_c', 'CGGA-325, WHO II–III',
              leg='upper right', note='')
    P.annot(c, 'n = %d\nlog-rank $P$ = %.1e' % (len(P.GGA), pc), 'sw', fs=6.6)
    return fig, {'A': [a], 'B': [b], 'C': [c]}


def build_F8():
    """分子相关性：A PCA、B 免疫（读 Table_S26c）、C 突变。

    v6.0：原 build_S2，图号由 Supplementary FigS2 升为正文 Figure 8。
    """
    import json
    import numpy as np
    from lifelines.statistics import multivariate_logrank_test
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score
    fig = plt.figure(figsize=SIZE['分子相关性'])
    # v6.0：hspace 由 0.50 放到 0.72——A 面板的图例从绘图区内部搬到了轴下方，
    # 需要给这一行腾出高度（sir 指出 S2 有遮挡：白底图例正盖在 PCA 散点上）。
    # v6.1：wspace 0.36 → 0.60，height_ratios 改成 [0.80, 1.20]，hspace 回调到 0.52。
    #   ① wspace：实测 A 的墨迹右沿 0.4687、B 的左沿 0.4257，**A 最右侧的红点
    #      与 B 的 y 轴标签（Interferon signature / Immune checkpoints）在水平
    #      方向交了 0.318 in**——sir 说的 "s2 有重叠" 就是这个。
    #      解得 wspace ≥ 0.53 才分得开，取 0.60 留余量。
    #   ② height_ratios：把省下的高度给下面板 C（见 PLAN 处的说明）。
    gs = GridSpec(2, 2, height_ratios=[0.80, 1.20], hspace=0.52, wspace=0.78)
    gs.update(bottom=0.115)

    X = P.EX[P.G25].values
    p = PCA(n_components=2).fit(X); Z = p.transform(X)
    km = KMeans(n_clusters=2, n_init=50, random_state=42).fit(X)
    lab = km.labels_
    if np.bincount(lab)[0] < np.bincount(lab)[1]:
        lab = 1 - lab
    sil = silhouette_score(X, lab); sz = np.bincount(lab)
    dd = P.COH.set_index('upk').loc[P.EX.index].copy(); dd['cl'] = lab
    lr = multivariate_logrank_test(dd['t'], dd['cl'], dd['e'])
    a = fig.add_subplot(gs[0, 0])
    for cc, col, nm in ((0, P.CB, 'C1'), (1, P.CR, 'C2')):
        m = lab == cc
        a.scatter(Z[m, 0], Z[m, 1], s=3.6, color=col, alpha=0.75, linewidths=0,
                  label='%s (n = %d, %d deaths)' % (nm, sz[cc], int(dd.loc[dd['cl'] == cc, 'e'].sum())))
    a.set_xlabel('PC1 (%.1f%%)' % (p.explained_variance_ratio_[0] * 100), fontsize=7.6)
    a.set_ylabel('PC2 (%.1f%%)' % (p.explained_variance_ratio_[1] * 100), fontsize=7.6)
    # v6.0 修 S2A 遮挡：v5.0 把 silhouette/log-rank 说明放在 'sw'、图例放在 'lower right'，
    # 两处都带白底，而 PCA 散点铺满整个下半区，于是白底把散点挖掉两块
    # （sir 说的 "s2 有遮挡"）。现在：说明性文字并进标题第二行（标题在轴外，不可能压数据），
    # 图例移到轴下方的空白带，绘图区内不再有任何说明文字或白底。
    P.title(a, '$k$-means subtypes of the 25 Model A genes\n'
               'silhouette = %.3f; log-rank $P$ = %.1e' % (sil, lr.p_value))
    a.legend(frameon=False, fontsize=6.2, loc='upper center', bbox_to_anchor=(0.5, -0.22),
             ncol=2, borderaxespad=0.0, handletextpad=0.5, columnspacing=1.6)
    a.grid(alpha=0.12, lw=0.5); P.cleantop(a); P.maxticks(a, 4)

    # B 免疫：改为读 Table_S26c，并与 Table 5 同口径
    b = fig.add_subplot(gs[0, 1])
    I = P.S26C.sort_values('difference_high_minus_low').reset_index(drop=True)
    yy = np.arange(len(I))
    b.barh(yy, I['difference_high_minus_low'], height=0.62, color=P.CR, alpha=0.88)
    for i, r in I.iterrows():
        st = '***' if r['P_MannWhitney'] < 1e-3 else '**' if r['P_MannWhitney'] < 1e-2 else '*'
        b.text(r['difference_high_minus_low'] + 0.010, i, '%.3f %s' % (r['difference_high_minus_low'], st),
               va='center', fontsize=6.2)
    b.set_yticks(yy); b.set_yticklabels(I['Immune_feature'], fontsize=6.2)
    b.set_xlim(0, I['difference_high_minus_low'].max() * 1.36)
    b.set_xlabel('Difference in feature score (high − low); all 11 higher in high-risk', fontsize=6.9)
    b.grid(axis='x', alpha=0.15, lw=0.5); b.set_axisbelow(True); P.cleantop(b)
    P.title(b, 'Immune features, Model A split (n = 255 per group)')

    # C 突变
    S = json.load(open('mut_final.json')); rows = S['rows']
    nh, nl = S['n_hi'], S['n_lo']
    genes = [r['gene'] for r in rows]
    order = sorted(range(len(genes)), key=lambda i: -(rows[i]['ph'] - rows[i]['pl']))
    genes = [genes[i] for i in order]; rows = [rows[i] for i in order]
    hi = np.array([r['ph'] for r in rows]); lo = np.array([r['pl'] for r in rows])
    pv = [r['p'] for r in rows]
    star = lambda p: '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'n.s.'
    c = fig.add_subplot(gs[1, :])
    y = np.arange(len(genes)); h = 0.38
    c.barh(y + h / 2, hi, height=h, color=P.CR, label='High-risk (n = %d)' % nh)
    c.barh(y - h / 2, lo, height=h, color=P.CB, label='Low-risk (n = %d)' % nl)
    for i, (x1, x2) in enumerate(zip(hi, lo)):
        c.text(x1 + 1.4, i + h / 2, '%.1f' % x1, va='center', fontsize=6.0, color=P.CR)
        c.text(x2 + 1.4, i - h / 2, '%.1f' % x2, va='center', fontsize=6.0, color=P.CB)
        c.text(105, i, star(pv[i]), va='center', ha='left', fontsize=6.4,
               color='#222' if pv[i] < 0.05 else P.CG)
    c.text(105, -0.98, 'Fisher $P$', va='center', ha='left', fontsize=6.4)
    c.set_yticks(y); c.set_yticklabels(genes, fontsize=7.2)
    c.set_ylim(len(genes) - 0.5, -1.40); c.set_xlim(0, 118)
    c.set_xticks([0, 20, 40, 60, 80, 100])
    c.set_xlabel('Patients with a non-synonymous mutation (%)', fontsize=7.6)
    c.grid(axis='x', alpha=0.15, lw=0.5); c.set_axisbelow(True)
    c.legend(frameon=False, fontsize=6.6, loc='upper center', bbox_to_anchor=(0.5, -0.145), ncol=2)
    P.cleantop(c)
    P.title(c, 'Somatic mutation frequencies by Model A risk group (n = 510)')
    return fig, {'A': [a], 'B': [b], 'C': [c]}


# v6.0：build_S3（IDH 分层风险评分箱线图）整体删除。
# 它是 P.panel_2C 的放大副本，与 Figure 2C 逐字节相同，v5.0 已从交付集撤下；
# 留在代码里只会让 PLAN['S3'] 成为一个不存在键的悬挂引用。


# v5.3：键 = 新图号，值为**内容不变**的原绘制函数。3↔4、6→8/7→6/8→7 的重排在
# PLAN 的名字里完成，这里显式交叉指回，避免有人误以为 build_F3 画的是新 Figure 3。
BUILDERS = {'F1': build_F1, 'F2': build_F2, 'F3': build_F4, 'F4': build_F3,
            'F5': build_F5, 'F6': build_F7, 'F7': build_F8, 'F8': build_F6}

if __name__ == '__main__':
    import json
    import os
    which = sys.argv[1:] or list(PLAN)
    LAYOUT = {}
    for k in which:
        if k not in BUILDERS:
            print('跳过未知图号', k); continue
        name = PLAN[k][0]
        fig, axs = BUILDERS[k]()
        fig.canvas.draw()
        # 顺序很关键：
        #   ① 先导面板（此时图上还没有面板字母，导出的单幅 PDF 是干净的矢量面板）；
        #   ② 再加面板字母；③ 再导整图。
        # save_panel 用的是整幅画布尺寸，所以每张面板 PDF 与整图逐像素对齐，
        # Illustrator 里把每张都放在原点、100% 缩放即可复原版面。
        for L, ax_list in axs.items():
            P.save_panel(ax_list, '%s_%s' % (name, L))
        for L, ax_list in axs.items():
            P.panellabel(fig, ax_list[0], L)
        fig.canvas.draw()
        # 记录每张面板的坐标轴矩形与面板字母位置（图幅比例），供 .jsx 放活字字母
        rec = {'figsize_in': list(fig.get_size_inches()),
               'file': name, 'panels': {}}
        for L, ax_list in axs.items():
            bb = ax_list[0].get_position()
            rec['panels'][L] = {'x0': bb.x0, 'y0': bb.y0, 'x1': bb.x1, 'y1': bb.y1}
        # 面板字母的 fig.text 是最后加入 fig.texts 的若干个，按加入顺序对应 axs
        # v5.0 关键修正：Illustrator 的 tf.position 是文本框**左上角**，且 y 向下为负。
        # 因此除了 matplotlib 的锚点 y（va='bottom' 时等于 bbox 下沿）之外，
        # 还要记录 bbox **上沿** letter_ytop，.jsx 才能把字母放到与 matplotlib 预览一致的位置。
        letters = [t for t in fig.texts if t.get_text() in list(axs.keys())]
        rend = None
        try:
            rend = fig.canvas.get_renderer()
        except Exception:
            rend = None
        fig_h_in = fig.get_size_inches()[1]
        for t in letters:
            pos = t.get_position()
            k = t.get_text()
            rec['panels'][k]['letter_x'] = pos[0]
            rec['panels'][k]['letter_y'] = pos[1]
            rec['panels'][k]['letter_fs'] = t.get_fontsize()
            ytop = None
            if rend is not None:
                try:
                    tb = t.get_window_extent(rend).transformed(fig.transFigure.inverted())
                    ytop = float(tb.y1)
                    # v6.0：把字母的**完整墨迹矩形**也记下来。Illustrator 只用到
                    # 左上角（letter_x / letter_ytop），而 check_panel_overlap.py
                    # 要判断「面板字母有没有压进邻面板」，必须知道右下沿。
                    rec['panels'][k]['letter_box'] = [float(tb.x0), float(tb.y0),
                                                      float(tb.x1), float(tb.y1)]
                except Exception:
                    ytop = None
            if ytop is None or not (0.0 < ytop <= 1.0):
                # 兜底：按字号估算行高（DejaVu Sans 约 1.19 em）
                ytop = min(pos[1] + 1.19 * t.get_fontsize() / 72.0 / fig_h_in, 0.995)
            rec['panels'][k]['letter_ytop'] = float(ytop)
        LAYOUT[name] = rec
        P.savefig(fig, name)
        plt.close(fig)
    # v6.0 修：过去这里直接覆盖写 layout.json，于是 "python make_figs_v5.py F3"
    # 会把其余 7 张图的版面记录全部抹掉，而 Illustrator 的 .jsx 正是靠它定位面板字母。
    # 现在改为**合并**：只更新本次跑过的图，其余键原样保留。
    LP = '/tmp/mitoxy/figs_v5/layout.json'
    if os.path.exists(LP):
        try:
            old = json.load(open(LP, encoding='utf-8'))
        except Exception:
            old = {}
        old.update(LAYOUT)
        LAYOUT = old
    json.dump(LAYOUT, open(LP, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
    print('\n'.join(P.LOG))
    print('layout -> %s（本次更新 %d 张，文件内共 %d 张）' % (LP, len(which), len(LAYOUT)))
