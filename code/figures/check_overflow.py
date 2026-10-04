# -*- coding: utf-8 -*-
"""系统性版面体检：找出所有「越出画布」与「文字互相遮挡」的元素。

对每张图：
  ① 越界：任何 Text / Legend 的窗口外框超出 figure 边界的比例 > tol_in 即报；
  ② 跨面板压字：A 面板的文字落进了 B 面板的坐标轴矩形内（排除该面板自身元素）；
  ③ 同图文字相撞：任意两段文字的窗口外框（各向内缩 22% 抵消行距假重叠）交叠 > 30%。

每处都带 axes 标签，便于直接定位到代码。

用法：python check_overflow.py            # 全部
      python check_overflow.py F1 F6
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import make_figs_v5 as M

TOL_IN = 0.02
SHRINK = 0.22
COLL_TH = 0.30
CROSS_TH = 0.35


def bbox_of(fig, artist, rend):
    try:
        bb = artist.get_window_extent(rend)
    except Exception:
        return None
    if bb is None or bb.width <= 0 or bb.height <= 0:
        return None
    return bb.transformed(fig.transFigure.inverted())


def label_of(t):
    if hasattr(t, 'get_texts'):
        s = 'legend[' + ' | '.join(x.get_text() for x in t.get_texts()) + ']'
    else:
        s = t.get_text() or ''
    s = s.replace('\n', '\\n')
    return (s[:34] + '…') if len(s) > 34 else s


def collect(fig):
    """→ [(artist, ax, tag)]，tag 标明这是哪个轴上的什么东西。"""
    out = []
    for i, ax in enumerate(fig.axes):
        tag = 'ax%d' % i
        if not ax.get_visible():
            continue
        for lbl, lst in (('xlabel', [ax.xaxis.label]), ('ylabel', [ax.yaxis.label]),
                         ('title', [ax.title])):
            for t in lst:
                if t.get_text() and t.get_visible():
                    out.append((t, ax, '%s.%s' % (tag, lbl)))
        for t in ax.texts:
            if t.get_text() and t.get_visible():
                out.append((t, ax, '%s.text' % tag))
        # 刻度：matplotlib 对越出视野的刻度仍报 visible=True，需自己按视野区间过滤
        xlo, xhi = sorted(ax.get_xlim())
        ylo, yhi = sorted(ax.get_ylim())
        if ax.xaxis.get_visible():
            for tk in ax.xaxis.get_major_ticks():
                loc = tk.get_loc()
                if not (xlo - 1e-9 <= loc <= xhi + 1e-9):
                    continue
                if tk.label1.get_text() and tk.label1.get_visible():
                    out.append((tk.label1, ax, '%s.xtick' % tag))
        if ax.yaxis.get_visible():
            for tk in ax.yaxis.get_major_ticks():
                loc = tk.get_loc()
                if not (ylo - 1e-9 <= loc <= yhi + 1e-9):
                    continue
                if tk.label1.get_text() and tk.label1.get_visible():
                    out.append((tk.label1, ax, '%s.ytick' % tag))
        lg = ax.get_legend()
        if lg is not None and lg.get_visible():
            out.append((lg, ax, '%s.legend' % tag))
    for t in fig.texts:
        if t.get_text() and t.get_visible():
            out.append((t, None, 'figtext'))
    return out


def check(key):
    name = M.PLAN[key][0]
    fig, axs = M.BUILDERS[key]()
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W_in, H_in = fig.get_size_inches()
    items = collect(fig)

    full = [Bbox.unit()] + [Bbox.unit()]      # 占位，下面直接用
    panel_own = {}
    panel_box = {}
    for L, al in axs.items():
        panel_own[L] = set(id(x) for x in al if x is not None)
        panel_box[L] = al[0].get_position()

    # ① 越界
    oob = []
    for t, ax, tag in items:
        bb = bbox_of(fig, t, rend)
        if bb is None:
            continue
        dl, dr = -bb.x0 * W_in, (bb.x1 - 1.0) * W_in
        db, dt = -bb.y0 * H_in, (bb.y1 - 1.0) * H_in
        side = []
        if dl > TOL_IN: side.append('左 %.3f in' % dl)
        if dr > TOL_IN: side.append('右 %.3f in' % dr)
        if db > TOL_IN: side.append('下 %.3f in' % db)
        if dt > TOL_IN: side.append('上 %.3f in' % dt)
        if side:
            oob.append((tag, label_of(t), ', '.join(side)))

    # ② 跨面板压字
    cross = []
    for t, ax, tag in items:
        if ax is None:
            continue
        bb = bbox_of(fig, t, rend)
        if bb is None:
            continue
        for L, pb in panel_box.items():
            if id(ax) in panel_own[L] or pb is None:
                continue
            x0, x1 = max(bb.x0, pb.x0), min(bb.x1, pb.x1)
            y0, y1 = max(bb.y0, pb.y0), min(bb.y1, pb.y1)
            if x1 <= x0 or y1 <= y0:
                continue
            frac = (x1 - x0) * (y1 - y0) / max(bb.width * bb.height, 1e-12)
            if frac > CROSS_TH:
                cross.append((tag, label_of(t), L, frac))

    # ③ 同图文字相撞
    def shrink(bb):
        return Bbox.from_extents(bb.x0 + bb.width * SHRINK, bb.y0 + bb.height * SHRINK,
                                 bb.x1 - bb.width * SHRINK, bb.y1 - bb.height * SHRINK)
    boxes = []
    for t, ax, tag in items:
        bb = bbox_of(fig, t, rend)
        if bb is not None:
            boxes.append((shrink(bb), label_of(t), tag))
    coll = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, la, ta = boxes[i]
            b, lb, tb = boxes[j]
            x0, x1 = max(a.x0, b.x0), min(a.x1, b.x1)
            y0, y1 = max(a.y0, b.y0), min(a.y1, b.y1)
            if x1 <= x0 or y1 <= y0:
                continue
            inter = (x1 - x0) * (y1 - y0)
            fr = max(inter / max(a.width * a.height, 1e-12),
                     inter / max(b.width * b.height, 1e-12))
            if fr > COLL_TH:
                coll.append((ta, la, tb, lb, fr))

    print('=' * 78)
    print('%s  %s   %.1f x %.1f in' % (key, name, W_in, H_in))
    if oob:
        print('  ① 越出画布 %d 处：' % len(oob))
        for tag, s, where in oob:
            print('       [%s] "%s"  → %s' % (tag, s, where))
    if cross:
        print('  ② 压到别的面板 %d 处：' % len(cross))
        for tag, s, L, fr in cross:
            print('       [%s] "%s" 落进面板 %s（占其面积 %.0f%%）' % (tag, s, L, fr * 100))
    if coll:
        print('  ③ 文字相撞 %d 对：' % len(coll))
        for ta, la, tb2, lb, fr in coll:
            print('       [%s] "%s"  ×  [%s] "%s"  (%.0f%%)' % (ta, la, tb2, lb, fr * 100))
    if not (oob or cross or coll):
        print('  干净：无越界、无压字、无相撞')
    plt.close(fig)


if __name__ == '__main__':
    which = sys.argv[1:] or list(M.PLAN)
    for k in which:
        if k in M.BUILDERS:
            check(k)
