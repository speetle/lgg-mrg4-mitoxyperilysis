# -*- coding: utf-8 -*-
"""拼版兜底通道（无 Illustrator Apple Events 时）：用 PyMuPDF 生成 PDF 兼容 .ai。

为什么需要它
------------
sir 要求「一定要用 ai 来拼接 fig」。本机 Illustrator 30.1 的 Apple Events 通道
被 macOS「自动化」权限挡住——`osascript -e 'tell application "Adobe Illustrator"
to do javascript …'` 一律返回 `-10004 权限违例`，而 `get version` 这类静态属性
仍能读，所以极易被误判成脚本 bug。同一台机上 `tell application "Finder" to
count windows` 也报 -10004，说明是整机授权缺失。

本脚本给出**功能等效**的产物，逐条对应 sir 的诉求：

  · 「每个图都不要和其他的重叠」 —— 每张面板 PDF 本身就是整幅画布尺寸、已在
    原点，逐一叠加即精确复原 matplotlib 版面；重叠与否由 check_panel_overlap.py
    在墨迹级独立判定。
  · 「调整好位置」 —— 位置全部来自 layout.json，与 matplotlib 合成件同源。
  · 「可编辑」 —— 每张面板落在**独立 OCG 图层**上，Illustrator 打开后在「图层」
    面板里可逐个选中／移动／隐藏；面板字母是**真文本**（嵌入 DejaVu Sans Bold），
    可双击改字、改字号。

⚠️ 这不是 Illustrator 原生脚本产出（没有 AIPrivateData / Roundtrip 块）。
授权恢复后，用 `compose_ai.jsx` 重跑即可覆盖同名文件。

产出：figs_v5/ai/<name>.ai / .pdf / .png（600 dpi）
      并把 .ai 的渲染图与 matplotlib 合成图做像素级交叉校验。
用法：python compose_ai_fallback.py
"""
import json
import os
import shutil

import fitz
import numpy as np
from PIL import Image

import matplotlib

WD = '/tmp/mitoxy'
LAY = os.path.join(WD, 'figs_v5', 'layout.json')
PANEL = os.path.join(WD, 'figs_v5', 'panels')
OUT = os.path.join(WD, 'figs_v5', 'ai')
TTF = os.path.join(os.path.dirname(matplotlib.__file__),
                   'mpl-data', 'fonts', 'ttf', 'DejaVuSans-Bold.ttf')
PX_PER_PT = 600 / 72.0           # 600 dpi 渲染

_CAL_FS = 10.4


def calibrate_asc():
    """标定「`letter_ytop` 到基线」的距离，单位 em。

    为什么必须实测而不是查表：`make_figs_v5` 记录的 `letter_ytop` 来自
    `Text.get_window_extent()`，而 matplotlib 给文本用的**不是**字体的 hhea
    ascender（DejaVu Sans Bold = 0.9282 em），也**不是** cap height（0.7269 em），
    而是排版引擎自己算出来的一个小盒子——2026-10-04 实测为 **0.7462 em**。
    前两版分别按 0.729 和 0.9282 写死，字母整体上/下偏 3–15 px，在
    `compare_composites.py` 的 1:1 对照里露了两次。

    所以改成在这里现量：用 matplotlib 自己画一个 'A'，量出 `letter_ytop`
    到字形下沿（= 基线，'A' 无降部）的像素距离，换算成 em。字体、字号或
    matplotlib 版本变了也会自动跟上。
    """
    import io as _io
    import numpy as _np
    from PIL import Image as _Image
    import matplotlib.pyplot as _plt
    fig = _plt.figure(figsize=(3, 1), dpi=600)
    t = fig.text(0.20, 0.30, 'A', fontsize=_CAL_FS, fontweight='bold',
                 ha='left', va='bottom')
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer()).transformed(
        fig.transFigure.inverted())
    H = fig.get_size_inches()[1] * 600
    ytop = (1 - bb.y1) * H
    buf = _io.BytesIO()
    fig.savefig(buf, format='png', dpi=600, bbox_inches=None)
    _plt.close(fig)
    buf.seek(0)
    img = _np.asarray(_Image.open(buf).convert('L'))
    ys, _ = _np.nonzero(img < 128)
    return (ys.max() - ytop) / (_CAL_FS * 600 / 72)


ASC = None      # 首次使用时标定


def build(name, rec):
    global ASC
    if ASC is None:
        ASC = calibrate_asc()
    W_pt = rec['figsize_in'][0] * 72.0
    H_pt = rec['figsize_in'][1] * 72.0
    doc = fitz.open()
    page = doc.new_page(width=W_pt, height=H_pt)

    letters = sorted(rec['panels'])
    drawn = []
    for L in letters:
        p = os.path.join(PANEL, '%s_%s.pdf' % (name, L))
        if not os.path.exists(p):
            raise SystemExit('缺面板 %s' % p)
        oci = doc.add_ocg('面板 %s' % L)
        pd = fitz.open(p)
        page.show_pdf_page(page.rect, pd, 0, keep_proportion=True, oc=oci)
        pd.close()
        drawn.append(L)

    # ⚠️ 坐标系陷阱（第一版就在这里翻了车）：PyMuPDF 的 `insert_text` 用的是
    # **以页顶为原点、y 向下** 的坐标（`insert_text((72, 72))` 落在页面左上角），
    # 与 PDF 内层的 y 向上坐标系相反。而 layout.json 的 `letter_ytop` 是
    # matplotlib 的图幅分数、**y 自下往上**。第一版直接写成 `ytop * H`，结果整排
    # 字母被上下镜像（A/B 跑到页底、C/D 互换观感），差分校验才发现。
    #
    # 正确换算：基线距页顶 = (1 − letter_ytop) × H + ascender/em × 字号。
    page.insert_font(fontname='DJVB', fontfile=TTF)
    for L in letters:
        d = rec['panels'][L]
        fs = d['letter_fs']
        x = d['letter_x'] * W_pt
        y_top_down = (1.0 - d['letter_ytop']) * H_pt
        page.insert_text((x, y_top_down + ASC * fs), L, fontname='DJVB',
                         fontsize=fs, color=(0, 0, 0))
    # OCG 的可见性由 add_layer 默认打开；这里不做额外开关。

    doc.set_metadata({'title': name, 'creator': 'PyMuPDF composite (AI-compatible PDF)',
                      'producer': 'compose_ai_fallback.py'})
    os.makedirs(OUT, exist_ok=True)
    for ext in ('.ai', '.pdf'):
        with open(os.path.join(OUT, name + ext), 'wb') as f:
            f.write(doc.tobytes(garbage=4, deflate=True, clean=True))
    pm = page.get_pixmap(matrix=fitz.Matrix(PX_PER_PT, PX_PER_PT), alpha=False)
    png = os.path.join(OUT, name + '.png')
    pm.save(png)
    doc.close()
    return png, drawn


def _ink_crop(a, thr=250):
    """裁到非白墨迹的外接矩形。"""
    ys, xs = np.nonzero(a < thr)
    if len(xs) == 0:
        return None
    return a[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def cross_check(name, png):
    """把兜底 .ai 的渲染图与 matplotlib 合成图比一遍。

    两者版面同源，但**裁剪口径不同**：本脚本按 Illustrator 惯例导整幅画布
    （7.4 × 7.7 in → 4440 × 4620 px），matplotlib 的 savefig 用的是
    `bbox_inches='tight', pad_inches=0.06`（→ 4038 × 3982 px）。所以不能直接比
    像素尺寸，改用两个与裁剪无关的量：

      ① **墨迹外接矩形的长宽**。同一份版面在两种裁剪下，墨迹矩形是同一个物理
         区域，尺寸必须一致（600 dpi 下容许 ±4 px，即 ±0.007 in）。这一项才是
         「面板位置对不对」的判据——第一版把面板字母上下镜像时，高度差了 96 px，
         一眼就露。
      ② 在 ±8 px 内做一次整数平移搜索后的平均灰度差，用来兜住渲染器抗锯齿与
         字形的亚像素差异。不做平移搜索的话，线稿只差 1 px 就会报出 20% 的
         「差异」，那是假阳性。
    """
    a = np.asarray(Image.open(png).convert('L'), dtype=np.float32)
    b = np.asarray(Image.open(os.path.join(WD, 'figs_v5', name + '.png'))
                   .convert('L'), dtype=np.float32)
    ia, ib = _ink_crop(a), _ink_crop(b)
    if ia is None or ib is None:
        return None, '整幅空白'
    dw = ia.shape[1] - ib.shape[1]
    dh = ia.shape[0] - ib.shape[0]
    # ① 墨迹矩形尺寸
    if abs(dw) > 4 or abs(dh) > 4:
        return None, ('墨迹矩形不符：宽 %d vs %d（差 %d px）、高 %d vs %d（差 %d px）'
                      % (ia.shape[1], ib.shape[1], dw, ia.shape[0], ib.shape[0], dh))
    # ② 在墨迹矩形内做平移搜索后的平均灰度差
    h = min(ia.shape[0], ib.shape[0]); w = min(ia.shape[1], ib.shape[1])
    A = ia[:h, :w]; B = ib[:h, :w]
    best = 255.0
    for dy in range(-8, 9):
        for dx in range(-8, 9):
            y0, y1 = max(0, dy), min(h, h + dy)
            x0, x1 = max(0, dx), min(w, w + dx)
            if y1 - y0 < h - 20 or x1 - x0 < w - 20:
                continue
            m = float(np.abs(A[y0 - dy:y1 - dy, x0 - dx:x1 - dx] - B[y0:y1, x0:x1]).mean())
            best = min(best, m)
    return best, ('墨迹 %d×%d（差 %+d, %+d px）；平移搜索后平均灰度差 %.2f'
                  % (ia.shape[1], ia.shape[0], dw, dh, best))


def main():
    lay = json.load(open(LAY, encoding='utf-8'))
    print('%-34s %-4s %s' % ('图', '面板', '与 matplotlib 合成件比对（裁到墨迹后）'))
    bad = []
    for name in sorted(lay):
        png, drawn = build(name, lay[name])
        frac, info = cross_check(name, png)
        ok = frac is not None and frac < 12.0
        print('%-34s %-4s %s  %s' % (name, ''.join(drawn), 'OK ' if ok else 'FAIL', info))
        if not ok:
            bad.append(name)
    print()
    # 产地说明：这些 .ai 不是 Illustrator 原生脚本产出，写清楚，免得日后混淆。
    with open(os.path.join(OUT, '_PROVENANCE.txt'), 'w', encoding='utf-8') as f:
        f.write(
            '本目录的 .ai / .pdf / .png 由 /tmp/mitoxy/compose_ai_fallback.py 生成\n'
            '（PyMuPDF，%s）。\n\n'
            '性质：PDF 兼容的 Illustrator 文档。每张面板是一个独立对象、各占一个\n'
            'OCG 图层；面板字母是真文本（嵌入 DejaVu Sans Bold），可双击编辑。\n'
            '**没有 Illustrator 私有块（AIPrivateData / Roundtrip）**，即不是\n'
            'Illustrator 原生脚本产出。\n\n'
            '成因：本机 Illustrator 30.1 的 Apple Events 通道被 macOS「自动化」\n'
            '权限挡住（-10004），`do javascript` 无法执行。授权恢复后在 Illustrator\n'
            '里跑 compose_ai.jsx 即可覆盖同名文件，得到原生产物。\n'
            % __import__('datetime').datetime.now().isoformat(timespec='seconds'))
    if bad:
        raise SystemExit('!! 兜底拼版与 matplotlib 合成件不一致：%s' % bad)
    print('8/8 通过：兜底 .ai 与 matplotlib 合成件的墨迹几何一致。')
    print('-> %s' % OUT)


if __name__ == '__main__':
    main()
