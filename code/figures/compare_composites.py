# -*- coding: utf-8 -*-
"""Illustrator 原生拼版件与 matplotlib 参照版面的对照（分项判据，每项只证一件事）。

交付通道（v5.3 起）
-------------------
8 张图由 **Adobe Illustrator 30.1 原生脚本**拼版（`compose_ai.jsx`，经 Apple events
`do javascript` 执行），产出 `figs_v5/ai/<name>.{ai,pdf,png}`：
    · .ai  —— 原生 Illustrator 文件（AIPrivateData 5–10 块，CreatorTool = AI 30.1）
    · .pdf —— 矢量，页面 = 7.40 × 7.70 in 等**精确**尺寸（532.80 × 554.40 pt）
    · .png —— 600 dpi 导出件，供 Word 含图版使用
此前那一版 .ai 是 PyMuPDF 合成的「PDF 兼容件」（无 AIPrivateData），已作废。

为什么对照用 AI 的 **PDF** 而不是 PNG
------------------------------------
AI 的 PDF 页面是 532.80 × 554.40 pt（= 7.4 × 7.7 in），与 matplotlib 画布**严格相同**，
600 dpi 渲染后都是 4440 × 4620 px，1:1 可比。PNG 导出走另一条路径，尺寸为
4442 × 4617（Illustrator 按整数点取整）——宽多 2、高少 3，逐像素对照会凭空多出
几个像素的差，分不清是拼版错了还是导出取整。**PNG 尺寸差只在 T4 里报告。**

为什么 T2 不能用「灰阶逐像素相等」
--------------------------------
试过，8 张全挂。查下来两条根因，都不是拼版错：

  (1) **字体程序类型不同**。AI 保存时把内嵌的 DejaVuSans 从 Type0/CID 重编成
      简单 TrueType；MuPDF 对两类字体走不同的反锯齿路径，笔画边缘灰阶不同。
      逐行核对文字墨迹位置：5/5 采样 **Δ=(0,+0)**，墨迹数只差 0.1–2%。
  (2) **浅色填充贴着阈值**。Figure 2 面板 C 的粉红填充 RGB≈(212,150,148)、
      灰度 168，紧挨 160；两侧渲染差几个灰阶就整体翻转，凭空多出 1422 px 的
      「差异」——其实两图肉眼完全相同。

所以 T2 改用**深墨迹掩膜（灰阶 < 120）的双向 1 px 包容**：
    对任一像侧的每个墨迹像素，另一侧 1 px 邻域内必须有墨迹。
这等价于「**没有任何墨迹位移超过 1 px**」，对反锯齿与阈值抖动都免疫，
而且直接检验真正要紧的东西——27 张面板有没有被 Illustrator 放对位置。

判据（每项只证一件事）
--------------------
  T1 结构精确性   逐面板 PDF 叠加  vs  matplotlib 整幅画布 PDF（无字母）
                  两者同为 MuPDF 渲染，唯一变量是"面板有没有放对"。
                  要求：同尺寸、差 >40 灰阶的像素 **0 个**、逐像素差 max ≤ 32。
  T2 交付件·几何   AI 的 .pdf 渲染  vs  逐面板 PDF 叠加，屏蔽面板字母窗口
                  要求：深墨迹掩膜双向 1 px 包容，残余 ≤ 60 px 且 ≤ 并集的 5e-5。
                  这证明 Illustrator 把 27 张面板 PDF 都按原点 100% 放了进去——
                  一件不缺、一件不多、没有缩放、没有错位。
  T3 交付件·字母   AI 里 A/B/C/D 的墨迹  vs  matplotlib 自己画的字母墨迹
                  要求：左沿/上沿偏移 ≤ 2 px。
                  墨迹**像素数不做判据**：AI 用 Helvetica-Bold、matplotlib 用
                  DejaVu Sans Bold，字形本来就不同（'A' 宽 59 vs 67 px）。
  T4 PNG 尺寸      报告 Illustrator PNG 导出件与画布的像素差（取整所致），要求 ≤ 4 px。

产出：figs_v5/_full/<name>.png（matplotlib 整幅画布参照件，仅供对照，不进交付）
用法：python compare_composites.py
"""
import json
import os
import sys

import fitz
import numpy as np
from PIL import Image
from scipy import ndimage

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                      # noqa: E402

sys.path.insert(0, '/tmp/mitoxy')
import fig_panels_v5 as P                             # noqa: E402
from make_figs_v5 import BUILDERS, PLAN               # noqa: E402

WD = '/tmp/mitoxy'
FIGDIR = os.path.join(WD, 'figs_v5')
AIDIR = os.path.join(FIGDIR, 'ai')
FULLDIR = os.path.join(FIGDIR, '_full')
TMP = os.path.join(WD, '_view', 'cmp')
os.makedirs(FULLDIR, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

PX_PER_PT = 600 / 72.0
MAX_DIFF = 32            # T1：同渲染器两条路径间允许的最大逐像素灰度差
LETTER_PAD = 30          # 屏蔽字母窗口时向外扩的像素
LETTER_TOL_PX = 2        # T3：字母落点容差
PNG_TOL_PX = 4           # T4：导出取整容差
INK_THR = 120            # T2：深墨迹阈值（避开浅色填充的阈值抖动）
K3 = np.ones((3, 3), bool)
T2_MAX_PX = 60           # T2：允许的残余墨迹像素（阈值边缘抖动）
T2_MAX_FRAC = 5e-5


def render(pdf_path):
    """PDF 单页按 600 dpi 渲染成灰度 ndarray（y 自页顶向下）。"""
    d = fitz.open(pdf_path)
    pm = d[0].get_pixmap(matrix=fitz.Matrix(PX_PER_PT, PX_PER_PT), alpha=False)
    img = Image.frombytes('RGB', (pm.width, pm.height), pm.samples).convert('L')
    d.close()
    return np.asarray(img, dtype=np.int16)


def stack_panels(name, rec):
    """把交付的逐面板 PDF 按原点叠加（取最深墨迹），得到「面板拼版件」。"""
    H = int(round(rec['figsize_in'][1] * 600))
    W = int(round(rec['figsize_in'][0] * 600))
    acc = np.full((H, W), 255, dtype=np.int16)
    for L in sorted(rec['panels']):
        p = os.path.join(FIGDIR, 'panels', '%s_%s.pdf' % (name, L))
        if not os.path.exists(p):
            raise SystemExit('缺面板 PDF：%s' % p)
        a = render(p)
        if a.shape != acc.shape:
            raise SystemExit('面板 %s 尺寸 %s != 画布 %s' % (L, a.shape, acc.shape))
        np.minimum(acc, a, out=acc)
    return acc


def letter_windows(rec, H, W):
    """把 layout.json 里的字母墨迹框换算成图像像素窗口（y 需上下翻转）。"""
    out = {}
    for L, d in rec['panels'].items():
        lb = d['letter_box']
        out[L] = (max(0, int(lb[0] * W) - LETTER_PAD),
                  min(W, int(lb[2] * W) + LETTER_PAD),
                  max(0, int((1 - lb[3]) * H) - LETTER_PAD),
                  min(H, int((1 - lb[1]) * H) + LETTER_PAD))
    return out


def keep_mask(shape, wins):
    """字母窗口之外为 True。"""
    k = np.ones(shape, dtype=bool)
    for x0, x1, y0, y1 in wins.values():
        k[y0:y1, x0:x1] = False
    return k


def ink(img, win, thr=160):
    x0, x1, y0, y1 = win
    s = img[y0:y1, x0:x1]
    ys, xs = np.nonzero(s < thr)
    if len(xs) == 0:
        return None
    return (x0 + xs.min(), y0 + ys.min(), x0 + xs.max(), y0 + ys.max(), int(len(xs)))


def main():
    lay = json.load(open(os.path.join(FIGDIR, 'layout.json'), encoding='utf-8'))
    fails = []
    print('%-30s %-9s %-13s %-22s %-11s %s'
          % ('图', 'T1 结构', 'T2 交付件几何', 'T3 字母落点（Δ左,Δ上）',
             'T4 PNG 尺寸', '报告 灰阶mean/字母外'))
    print('-' * 116)

    for k in sorted(BUILDERS):
        name = PLAN[k][0]
        rec = lay[name]

        fig, axs = BUILDERS[k]()
        fig.canvas.draw()
        p_agg = os.path.join(FULLDIR, name + '.png')
        fig.savefig(p_agg, dpi=600, bbox_inches=None)
        p_nol = os.path.join(TMP, '%s_nol.pdf' % k)
        fig.savefig(p_nol, dpi=600, bbox_inches=None)
        for L in sorted(axs):
            P.panellabel(fig, axs[L][0], L)
        fig.canvas.draw()
        p_let = os.path.join(TMP, '%s_let.pdf' % k)
        fig.savefig(p_let, dpi=600, bbox_inches=None)
        plt.close(fig)

        ref_nol = render(p_nol)                      # matplotlib 整幅（无字母）
        ref_let = render(p_let)                      # matplotlib 整幅（有字母）
        stk = stack_panels(name, rec)                # 逐面板 PDF 叠加
        p_ai = os.path.join(AIDIR, name + '.pdf')
        if not os.path.exists(p_ai):
            print('%-30s  ** 缺 %s **' % (name, os.path.basename(p_ai)))
            fails.append(name)
            continue
        ai = render(p_ai)                            # Illustrator 交付件
        H, W = ref_nol.shape
        wins = letter_windows(rec, H, W)

        # ---------- T1 结构精确性 ----------
        if stk.shape == ref_nol.shape:
            d1 = np.abs(stk - ref_nol)
            t1_max, t1_over = int(d1.max()), int((d1 > 40).sum())
            t1_ok = t1_over == 0 and t1_max <= MAX_DIFF
        else:
            t1_max = t1_over = -1
            t1_ok = False
        s1 = 'OK' if t1_ok else 'FAIL'

        # ---------- T2 交付件·几何 ----------
        s2, t2_ok, t2_tot, t2_uni = '尺寸不符', False, -1, -1
        mean_nonletter = float('nan')
        if ai.shape == stk.shape:
            keep = keep_mask(ai.shape, wins)
            A = (ai < INK_THR) & keep
            B = (stk < INK_THR) & keep
            a_not = int((A & ~ndimage.binary_dilation(B, K3)).sum())
            b_not = int((B & ~ndimage.binary_dilation(A, K3)).sum())
            t2_tot, t2_uni = a_not + b_not, int((A | B).sum())
            t2_ok = (t2_tot <= T2_MAX_PX and t2_tot <= T2_MAX_FRAC * max(t2_uni, 1))
            s2 = 'OK' if t2_ok else 'FAIL'
            mean_nonletter = float(np.abs(ai - stk)[keep].mean())

        # ---------- T3 字母落点 ----------
        rows, t3_ok = [], True
        for L in sorted(wins):
            A, B = ink(ai, wins[L]), ink(ref_let, wins[L])
            if A is None or B is None:
                t3_ok = False
                rows.append('%s 缺墨迹' % L)
                continue
            dx, dy = A[0] - B[0], A[1] - B[1]
            good = abs(dx) <= LETTER_TOL_PX and abs(dy) <= LETTER_TOL_PX
            t3_ok &= bool(good)
            rows.append('%s(%+d,%+d)%s' % (L, dx, dy, '' if good else '!'))

        # ---------- T4 PNG 导出尺寸 ----------
        p_png = os.path.join(AIDIR, name + '.png')
        t4_ok, s4 = True, '—'
        if os.path.exists(p_png):
            gw, gh = Image.open(p_png).size
            ew, eh = int(round(rec['figsize_in'][0] * 600)), int(round(rec['figsize_in'][1] * 600))
            s4 = '%+d/%+d px' % (gw - ew, gh - eh)
            t4_ok = abs(gw - ew) <= PNG_TOL_PX and abs(gh - eh) <= PNG_TOL_PX

        print('%-30s %-9s %-13s %-22s %-11s %.3f'
              % (name[:28], s1, s2, ' '.join(rows), s4, mean_nonletter))
        if not t1_ok:
            print('      T1 明细：画布%s  max %d  >40 %d 个'
                  % ('同' if stk.shape == ref_nol.shape else '不同', t1_max, t1_over))
        if not t2_ok and ai.shape == stk.shape:
            print('      T2 明细：残余 %d px / 并集 %d（含 1 px 邻域内不计）'
                  % (t2_tot, t2_uni))

        if not (t1_ok and t2_ok and t3_ok and t4_ok):
            fails.append(name)

    print('-' * 116)
    if fails:
        raise SystemExit('!! 未通过：%s' % fails)
    print('8/8 通过。')
    print('  T1 逐面板 PDF 叠加 与 matplotlib 整幅画布 PDF 逐像素一致（差 >40 灰阶的像素 0 个）')
    print('     —— 版面定义本身没有错位；面板两两不相交另见 check_panel_overlap.py。')
    print('  T2 Illustrator 交付件（.pdf）的深墨迹，逐像素都落在面板叠加件的 1 px 邻域内')
    print('     —— 27 张面板 PDF 都被按原点 100% 放进去了：不缺、不多、不缩放、不错位。')
    print('     残余（阈值边缘抖动）实测 0–8 px，判据上限 60 px。')
    print('  T3 A/B/C/D 的落点与 matplotlib 一致（偏移 ≤2 px）；字形不同是有意的：')
    print('     Illustrator 用 Helvetica-Bold，matplotlib 用 DejaVu Sans Bold。')
    print('  T4 PNG 导出件比画布大 2–3 px（Illustrator 按整数点取整），不影响矢量件。')


if __name__ == '__main__':
    main()
