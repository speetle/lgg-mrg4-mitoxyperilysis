# -*- coding: utf-8 -*-
"""版面遮挡审计器 v6.1 —— 专查「文字压在数据上」。

为什么必须新增这一件
--------------------
v5.0 已有 check_overflow.py，但它只查三种失效模式：
  ① 文字越出画布、② 文字侵入**邻面板**的坐标区、③ 文字与**别的文字**相撞。
sir 2026-10-04 指出的 Fig3A/C/D、Fig5A、Fig6、FigS2 遮挡，**三种一条都不命中**——
因为那批缺陷是「文字压在自己的数据上」：白底注释框把误差棒挖空、
图例横线从字里穿过、数值标签落进同组另一根柱子里。

v6.1 的关键改动：文字盒子改成**字形级**
--------------------------------------
v6.0 用 `Text.get_window_extent()`（排版盒）判定。排版盒含 ascent/descent，
比实际墨迹高约 25%。在密集的 barh 面板里，一个 6 pt 的标签排版盒 7.1 pt、
而它所在的那根柱只有 5.8 pt，于是**每一行都会报 14% 重叠**，其实字形根本没碰。
把阈值调大只是掩盖问题，正确做法是量准：用 `TextPath` 取字形轮廓的真实包围盒，
再把baseline 对齐回排版盒。旋转文本和含 LaTeX 的文本回退到排版盒（偏保守）。

判定口径
--------
· 柱／面积／置信带（有面积的）：**字形盒子**与图元包围盒的交集面积 ÷ 字形盒面积；
· 散点：落在字形盒内的点数；
· 折线：沿路径密采样后落在字形盒内的采样点数（≥2 才算压住）。
网格线**不算**数据（它由 axis 绘制，不在 ax.lines 里），因此不会被误报。

「合法压盖」白名单
------------------
KM 图的置信带铺满整个绘图区，图例只能落在带上——这是期刊常规做法，不是缺陷；
DCA 的决策曲线同理（axvspan 底色 alpha 0.07）。这些**显式登记**在白名单里，
审计时单独列出「已声明豁免」，既不让它们污染结论，也不把问题藏起来。

用法
----
    python check_occlusion.py            # 全部 8 张图
    python check_occlusion.py F3 F5      # 只查指定图
    python check_occlusion.py -v         # 连带打印豁免项与重叠百分比
退出码 0 = 无**未声明**的遮挡；1 = 有。
"""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection, PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath

import make_figs_v5 as B

# ------------------------------------------------------------------ 判定阈值
AREA_FRAC = 0.02      # 有面积图元：交集 ≥ 字形盒面积的 2% 即算压住
N_SCATTER = 1         # 散点：≥1 个点落在字形盒内即算压住
N_LINE = 2            # 折线：≥2 个采样点落在字形盒内即算压住

# ------------------------------------------------- 合法压盖（显式声明，附理由）
# ⚠️ 键用 **内容名**（Figure 名里 "_" 之后那一段），不要用 PLAN 键 'F1'…'F8'。
#    v5.3 按首次引用顺序重排过图号（3↔4、6→8/7→6/8→7），PLAN 键的含义跟着变了：
#    白名单里写 ('F7', …) 原本指「分层生存」，重排后 F7 是「分子相关性」——
#    于是那 4 条豁免**静默失效**，审计器把 10 处本来合法的压盖报成「未声明遮挡」。
#    按内容名索引后，图号怎么重排都不会断。
WHITELIST = [
    # (内容, 面板, 标签前缀, 理由)
    ('发现与ModelA', 'D', '图例:', 'KM 置信带铺满绘图区，图例只能落在带上（期刊常规）'),
    ('发现与ModelA', 'D', 'HR ', '同上：KM 面板的 HR / 中位生存注记'),
    ('MRG4外部验证', 'A', '图例:', 'KM 置信带铺满绘图区，图例只能落在带上（期刊常规）'),
    ('MRG4外部验证', 'A', 'n = ', '同上：KM 面板的风险人数/检验注记必然落在置信带上'),
    ('分层生存', 'A', '图例:', 'KM 置信带铺满绘图区，图例只能落在带上'),
    ('分层生存', 'A', 'n = ', '同上：KM 面板注记'),
    ('分层生存', 'B', '图例:', 'KM 置信带铺满绘图区，图例只能落在带上'),
    ('分层生存', 'B', 'n = ', '同上：KM 面板注记'),
    ('分层生存', 'C', '图例:', 'KM 置信带铺满绘图区，图例只能落在带上'),
    ('分层生存', 'C', 'n = ', '同上：KM 面板注记'),
    ('头对头与决策分析', 'B', '图例:', 'DCA 的 axvspan 是浅红高亮带（#f2d7d5，20% 混合），其下无数据'),
]


def _content(k):
    """PLAN 键 → 内容名（Figure 名里 "_" 之后那一段）。"""
    return B.PLAN[k][0].split('_', 1)[1]


def _wl(k, L, label):
    c = _content(k)
    for cc, p, pre, why in WHITELIST:
        if cc == c and p == L and label.startswith(pre):
            return why
    return None


# ------------------------------------------- 文字 × 文字（含轴标签、刻度标签）
TEXT_FRAC = 0.10      # 交集 ≥ 较小字形盒的 10% 即算相撞


def _tick_labels(axis, lo, hi):
    """只收**落在坐标轴可视区间内**的刻度标签。

    坑：matplotlib 会为区间外的刻度也建 Text 对象，而且 `get_visible()` 是 True。
    图 1A 的 x 轴是 log 刻度，取景 0.35–3.2，但 LogLocator 仍在 0.2 / 0.3 / 4 / 6
    处建了 `$\\mathdefault{3\\times10^{-1}}$` 这类标签——它们压根不会被画出来，
    却和真正的刻度标签 '0.4'、'3.0' 在屏幕坐标里撞在一起。
    不按 tick.get_loc() 过滤，就会凭空报出一批假相撞。
    """
    inv = lo > hi
    a, b = (hi, lo) if inv else (lo, hi)
    out = []
    for tk in list(axis.get_major_ticks()) + list(axis.get_minor_ticks()):
        lbl = tk.label1
        if not lbl.get_text().strip() or not lbl.get_visible():
            continue
        try:
            loc = float(tk.get_loc())
        except (TypeError, ValueError):
            continue
        if a - 1e-12 <= loc <= b + 1e-12:
            out.append(lbl)
    return out


def _text_artists(ax):
    """把一个坐标轴上所有**会显示的**文字收齐。

    注意刻度标签也算——sir 说的「字与字叠在一起」多半就是刻度标签互压，
    而它们既不在 ax.texts 里，也不是 legend。
    """
    out = [t for t in ax.texts if t.get_visible()]
    for t in (ax.title, ax.xaxis.label, ax.yaxis.label):
        if t is not None and t.get_visible():
            out.append(t)
    out += _tick_labels(ax.xaxis, *ax.get_xlim())
    out += _tick_labels(ax.yaxis, *ax.get_ylim())
    leg = ax.get_legend()
    if leg is not None:
        out += list(leg.get_texts())
    return [t for t in out if t.get_text().strip()]


def audit_text_text(ax, fig, rend, k, L, out, wl):
    ts = _text_artists(ax)
    if len(ts) < 2:
        return
    boxes = []
    for t in ts:
        b = ink_bbox(t, rend, fig)
        if b is None or b.width <= 0 or b.height <= 0:
            continue
        boxes.append((t.get_text().strip().replace('\n', ' ⏎ ')[:30], b))
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, bb = boxes[i][1], boxes[j][1]
            ov = _overlap_area(a, bb)
            if ov <= 0:
                continue
            frac = ov / min(a.width * a.height, bb.width * bb.height)
            if frac >= TEXT_FRAC:
                la, lb = boxes[i][0], boxes[j][0]
                why = _wl(k, L, '文字相撞:') or _wl(k, L, la)
                rec = ('%s-%s' % (B.PLAN[k][0], L),
                       '文字相撞: 「%s」 × 「%s」' % (la, lb),
                       '重叠 %.0f%%' % (frac * 100), why)
                (wl if why else out).append(rec)


def ink_bbox(t, rend, fig):
    """文字的**字形级**包围盒（像素，原点左下）。

    做法：`Text.get_window_extent()` 给出排版盒与 baseline 的关系，
    `TextPath` 给出字形相对 baseline 的轮廓范围，两者对接即得墨迹盒。
    旋转文本 / 含 `$...$` 的文本直接回退到排版盒（宁可保守一点）。

    ⚠️ 单位陷阱：`TextPath.get_extents()` 的坐标在**点**（point），要乘 dpi/72；
    而 `renderer.get_text_width_height_descent()` 返回的**已经是像素**，
    不能再乘。v6.1 初版把后者也乘了一次 dpi/72，descent 被放大 4.17 倍，
    baseline 抬到排版盒上方，于是 xlabel 的字形盒飘到刻度标签头上，
    凭空报出一批「Time (years) × 6」之类的假相撞（44 处里大半是这个）。
    """
    bb = t.get_window_extent(rend)
    s = t.get_text()
    if not s.strip():
        return None
    if (t.get_rotation() % 360) != 0 or '$' in s or '\\' in s:
        return bb
    try:
        fp = t.get_fontproperties()
        size = t.get_fontsize()
        _w_px, _h_px, d_px = rend.get_text_width_height_descent(s, fp, False)
        tp = TextPath((0, 0), s, size=size, prop=fp, usetex=False)
        e = tp.get_extents()
        k = fig.dpi / 72.0                 # 点 -> 像素
        base_y = bb.y0 + d_px              # descent 已是像素
        return matplotlib.transforms.Bbox.from_extents(
            bb.x0 + e.x0 * k, base_y + e.y0 * k,
            bb.x0 + e.x1 * k, base_y + e.y1 * k)
    except Exception:
        return bb


def _overlap_area(a, b):
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    return dx * dy if (dx > 0 and dy > 0) else 0.0


def _inside(pts, bb):
    if len(pts) == 0:
        return 0
    return int(np.sum((pts[:, 0] >= bb.x0) & (pts[:, 0] <= bb.x1) &
                      (pts[:, 1] >= bb.y0) & (pts[:, 1] <= bb.y1)))


def audit_axes(ax, fig, rend, k, L, out, wl):
    """检查单个坐标轴内 文字 × 数据 的重叠。"""
    texts = []
    for t in ax.texts:
        if t.get_text().strip():
            texts.append((t.get_text().strip().replace('\n', ' ⏎ ')[:44], t))
    leg = ax.get_legend()
    if leg is not None:
        for t in leg.get_texts():
            if t.get_text().strip():
                texts.append(('图例: ' + t.get_text().strip()[:34], t))
    if not texts:
        return

    # v6.1：数据一律被裁在坐标区里（matplotlib 默认 clip_on=True）。
    # 所以「轴外的文字」不可能被数据压住——3A / 3D 的图例搬到轴下方之后，
    # 数值折线在屏幕坐标里仍会穿到那片区域，若不裁就会误报。
    axbb = ax.get_window_extent(rend)

    rects, scats, polys, lines = [], [], [], []
    for p in ax.patches:
        if isinstance(p, Rectangle) and p.get_visible() and p is not ax.patch:
            rects.append(p)
    for c in ax.collections:
        if not c.get_visible():
            continue
        if isinstance(c, PathCollection):
            scats.append(c)
        elif isinstance(c, PolyCollection):
            polys.append(c)
    for ln in ax.lines:
        if not ln.get_visible():
            continue
        ls = ln.get_linestyle()
        # v6.1：`plot(x, y, 'o')` 这类**只有标记没有连线**的系列，linestyle 是 'None'。
        # 旧版把它当折线做直线插值，于是在两两标记之间凭空造出几十个采样点，
        # 报出 Fig3D「折线 33 采样点」的假遮挡。这类系列只能按点算。
        if ls in ('None', ' ', '', 'none') or ln.get_linewidth() == 0:
            scats.append(ln)
        else:
            lines.append(ln)

    trans = ax.transData

    def _clip(pts):
        m = ((pts[:, 0] >= axbb.x0) & (pts[:, 0] <= axbb.x1) &
             (pts[:, 1] >= axbb.y0) & (pts[:, 1] <= axbb.y1))
        return pts[m]

    for label, t in texts:
        tb = ink_bbox(t, rend, fig)
        if tb is None:
            continue
        tb = matplotlib.transforms.Bbox.intersection(tb, axbb)
        if tb is None or tb.width <= 0 or tb.height <= 0:
            continue
        tarea = tb.width * tb.height
        hits = []

        for p in rects + polys:
            try:
                pb = p.get_window_extent(rend)
            except Exception:
                continue
            ov = _overlap_area(tb, pb) / tarea
            if ov >= AREA_FRAC:
                hits.append('柱/面积 %.0f%%' % (ov * 100))

        for c in scats:
            if isinstance(c, Line2D):
                pts = trans.transform(np.column_stack(
                    [np.atleast_1d(c.get_xdata()), np.atleast_1d(c.get_ydata())]))
            else:
                off = c.get_offsets()
                if off is None or len(off) == 0:
                    continue
                pts = trans.transform(np.asarray(off))
            n = _inside(_clip(pts), tb)
            if n >= N_SCATTER:
                hits.append('散点 %d 点' % n)

        for ln in lines:
            xd, yd = ln.get_xdata(), ln.get_ydata()
            if len(xd) < 2:
                continue
            xs = np.linspace(0, 1, 240)
            try:
                x = np.interp(xs, np.linspace(0, 1, len(xd)), np.asarray(xd, float))
                y = np.interp(xs, np.linspace(0, 1, len(yd)), np.asarray(yd, float))
                pts = trans.transform(np.column_stack([x, y]))
            except Exception:
                continue
            n = _inside(_clip(pts), tb)
            if n >= N_LINE:
                hits.append('折线 %d 采样点' % n)

        if hits:
            why = _wl(k, L, label)
            rec = ('%s-%s' % (B.PLAN[k][0], L), label, '，'.join(hits), why)
            (wl if why else out).append(rec)


def main():
    argv = [a for a in sys.argv[1:] if a != '-v']
    verbose = '-v' in sys.argv
    which = argv or list(B.PLAN)
    found, excused = [], []
    for k in which:
        if k not in B.BUILDERS:
            print('跳过未知图号', k)
            continue
        fig, axs = B.BUILDERS[k]()
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        seen = set()
        for L, ax_list in axs.items():
            for ax in ax_list:
                if id(ax) in seen:
                    continue
                seen.add(id(ax))
                audit_axes(ax, fig, rend, k, L, found, excused)
                audit_text_text(ax, fig, rend, k, L, found, excused)
        plt.close(fig)

    print('=' * 74)
    print('面板内遮挡审计（字形级）：文字 × 数据 ＋ 文字 × 文字')
    print('=' * 74)
    if not found:
        print('通过：没有任何**未声明**的文字压数据。')
    else:
        cur = None
        for tag, label, how, _ in found:
            if tag != cur:
                print('\n[%s]' % tag); cur = tag
            print('   ✗ 「%s」 ← %s' % (label, how))
        print('\n未声明遮挡合计 %d 处。' % len(found))
    if excused and verbose:
        print('\n---- 已声明豁免（合法压盖，不阻断）----')
        cur = None
        for tag, label, how, why in excused:
            if tag != cur:
                print('  [%s]' % tag); cur = tag
            print('     · 「%s」 ← %s   ｜ %s' % (label, how, why))
    elif excused:
        print('\n（另有 %d 处已声明豁免，用 -v 查看）' % len(excused))
    return 1 if found else 0


if __name__ == '__main__':
    sys.exit(main())
