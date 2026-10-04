# -*- coding: utf-8 -*-
"""拼版重叠审计器（v6.0 新增）——保证「图与图不重叠」。

为什么需要它
------------
Illustrator 里的每张面板 PDF 都是**整幅画布尺寸、透明底**，全部放在原点。
所以它们在几何上永远是同一块矩形，永远"重叠"——用 AI 的 boundingBox 判断
是查不出问题的。sir 看到的"不同的图出现重叠"，实际是**墨迹互相压到**：
上面板的 x 轴标签伸进下面板的标题带、左面板的 y 轴标签伸进右面板的坐标区。

判定办法
--------
面板 PNG 是透明底、整幅画布尺寸、600 dpi 导出的——正是交给 AI 的那一份。
取每条 **alpha > 0 的像素包围盒**，就得到该面板的**真实墨迹矩形**（图幅比例）。
然后：
  ① 墨迹矩形必须落在 [0,1]² 内（否则出画布）；
  ② 任意两块面板的墨迹矩形**不得相交**（这是"不重叠"的硬保证）；
  ③ 面板字母（由 layout.json 的 letter_box 给出）不得落进**别的**面板的墨迹矩形。

输出最小间隙，便于判断版面是否过挤。

用法：python check_panel_overlap.py            退出码 0 = 无重叠
"""
import json
import os
import re
import sys

import numpy as np
from PIL import Image

PDIR = '/tmp/mitoxy/figs_v5/panels'
LAY = '/tmp/mitoxy/figs_v5/layout.json'

# v5.3：ORDER 从 layout.json 推导（按图号数值排序），不再硬编码。
# 图号重排过一次（v5.2 → v5.3），手抄的名字表会静默变成悬挂引用。
_LAY = json.load(open(LAY, encoding='utf-8'))
ORDER = sorted(_LAY, key=lambda n: int(re.match(r'Figure(\d+)', n).group(1)))


def ink_box(png):
    """返回 alpha>0 像素的包围盒，归一化为 (x0, y0, x1, y1)，原点在**左下**。"""
    im = Image.open(png)
    if im.mode != 'RGBA':
        im = im.convert('RGBA')
    a = np.array(im.getchannel('A'))
    ys, xs = np.nonzero(a > 8)          # 8/255 以下的抗锯齿残影不算墨迹
    if len(xs) == 0:
        return None
    W, H = im.size
    # 图像 y 向下，转成图幅比例时翻过来
    return (xs.min() / W, 1.0 - (ys.max() + 1) / H,
            (xs.max() + 1) / W, 1.0 - ys.min() / H)


def inter(a, b):
    dx = min(a[2], b[2]) - max(a[0], b[0])
    dy = min(a[3], b[3]) - max(a[1], b[1])
    return (dx, dy) if (dx > 0 and dy > 0) else (None, None)


def gap(a, b):
    """两块矩形的最小水平/垂直间隙（负值=重叠）。"""
    gx = max(b[0] - a[2], a[0] - b[2])
    gy = max(b[1] - a[3], a[1] - b[3])
    return gx, gy


def main():
    lay = json.load(open(LAY, encoding='utf-8'))
    problems = []
    print('=' * 78)
    print('面板墨迹重叠审计（Illustrator 拼版等价几何）')
    print('=' * 78)
    for name in ORDER:
        if name not in lay:
            problems.append('%s 不在 layout.json 中' % name)
            continue
        rec = lay[name]
        W_in, H_in = rec['figsize_in']
        boxes = {}
        for L in rec['panels']:
            png = os.path.join(PDIR, '%s_%s.png' % (name, L))
            if not os.path.exists(png):
                problems.append('%s-%s 缺面板 PNG' % (name, L))
                continue
            bb = ink_box(png)
            if bb is None:
                problems.append('%s-%s 面板全透明（无墨迹）' % (name, L))
                continue
            boxes[L] = bb
        keys = sorted(boxes)
        print('\n[%s]  %.2f x %.2f in' % (name, W_in, H_in))
        for L in keys:
            x0, y0, x1, y1 = boxes[L]
            oob = (x0 < -1e-6 or y0 < -1e-6 or x1 > 1 + 1e-6 or y1 > 1 + 1e-6)
            print('   %s  墨迹 x %.4f–%.4f  y %.4f–%.4f  (%.2f x %.2f in)%s'
                  % (L, x0, x1, y0, y1, (x1 - x0) * W_in, (y1 - y0) * H_in,
                     '   ← 越出画布!' if oob else ''))
            if oob:
                problems.append('%s-%s 墨迹越出画布' % (name, L))
        # 两两相交
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                A, Bb = boxes[keys[i]], boxes[keys[j]]
                dx, dy = inter(A, Bb)
                if dx is not None:
                    problems.append('%s：%s 与 %s 墨迹重叠 %.3f x %.3f in'
                                    % (name, keys[i], keys[j],
                                       dx * W_in, dy * H_in))
                    print('   ✗ %s × %s 重叠 %.3f x %.3f in'
                          % (keys[i], keys[j], dx * W_in, dy * H_in))
                else:
                    gx, gy = gap(A, Bb)
                    print('     %s|%s 间隙 水平 %.3f in / 垂直 %.3f in'
                          % (keys[i], keys[j], gx * W_in, gy * H_in))
        # 面板字母是否压进别的面板
        for L, bb in boxes.items():
            lb = rec['panels'][L].get('letter_box')
            if not lb:
                continue
            for M, ob in boxes.items():
                if M == L:
                    continue
                dx, dy = inter(lb, ob)
                if dx is not None:
                    problems.append('%s：字母 %s 压进面板 %s 的墨迹'
                                    % (name, L, M))
                    print('   ✗ 字母 %s 压进面板 %s' % (L, M))

    print('\n' + '=' * 78)
    if problems:
        print('发现 %d 个问题：' % len(problems))
        for p in problems:
            print('  ✗ ' + p)
        return 1
    print('全部通过：8 张图、每一块面板的墨迹都在画布内，且两两不相交，'
          '面板字母也没有压进邻面板。')
    return 0


if __name__ == '__main__':
    sys.exit(main())
