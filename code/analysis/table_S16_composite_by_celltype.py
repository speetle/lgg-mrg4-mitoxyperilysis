# -*- coding: utf-8 -*-
"""Table S16 —— MRG-4 复合评分按细胞类型（第三轮新增沉积）。

背景：审稿人 R2-M3 指出 §3.11 此前只按**单基因**报告表达，从未把稿件自己的
复合评分按细胞类型算过，因此"是恶性细胞信号、不是基质信号"这一说法缺少直接证据。

本脚本用单细胞矩阵 + 稿件冻结的四基因系数，按两种尺度各算一遍：
  口径 A（within-cell z）：先在**每个细胞**上把四个基因各自 z 标准化（跨全部 135,482
                          个细胞），按系数合成评分，再按细胞类型取均值。
  口径 B（type-level z） ：先按细胞类型取四基因均值，再把**每个基因的均值**跨 15 个
                          细胞类型 z 标准化，按系数合成。

两种尺度均为**总体标准差**（ddof = 0）。

输入：sl040_expr.csv（由 census_sl040.py 从 CELLxGENE SL040 解析）
输出：Table_S16_MRG4_score_by_celltype.csv

用法：python table_S16_composite_by_celltype.py
"""
import sys
import numpy as np
import pandas as pd

EXPR = 'sl040_expr.csv'
OUT = 'Table_S16_MRG4_score_by_celltype.csv'

# 稿件冻结的 MRG-4 系数（与 Table 6 / §2.4 一致）
GENES = ['BMP1', 'KIF15', 'TRAF3', 'PDGFA']
COEF = np.array([+0.258, +0.464, -0.238, +0.480])


def main():
    d = pd.read_csv(EXPR)
    missing = [g for g in GENES if g not in d.columns]
    if missing:
        raise SystemExit('表达矩阵缺少基因列：%s' % missing)

    n_cells = len(d)
    n_types = d['CellType'].nunique()

    # ---- 口径 A：逐细胞 z → 合成 → 按类型均值 ----
    zc = (d[GENES] - d[GENES].mean()) / d[GENES].std(ddof=0)
    d = d.assign(_score_cell=zc.values @ COEF)
    a = d.groupby('CellType')['_score_cell'].mean()

    # ---- 口径 B：类型均值 → 每基因跨类型 z → 合成 ----
    m = d.groupby('CellType')[GENES].mean()
    zb = (m - m.mean()) / m.std(ddof=0)
    b = pd.Series(zb.values @ COEF, index=m.index)

    # ---- 四基因原始计数均值 ----
    raw = m  # 已是每类型四基因均值（原始 UMI 计数）

    out = pd.DataFrame({
        'CellType': m.index,
        'n_cells': d.groupby('CellType').size().reindex(m.index).values,
        'MRG4_score_within_cell_z': [round(a[c], 6) for c in m.index],
        'MRG4_score_type_level_z': [round(b[c], 6) for c in m.index],
        'BMP1_mean_counts': [round(raw.loc[c, 'BMP1'], 6) for c in m.index],
        'KIF15_mean_counts': [round(raw.loc[c, 'KIF15'], 6) for c in m.index],
        'TRAF3_mean_counts': [round(raw.loc[c, 'TRAF3'], 6) for c in m.index],
        'PDGFA_mean_counts': [round(raw.loc[c, 'PDGFA'], 6) for c in m.index],
    })
    out = out.sort_values('MRG4_score_within_cell_z', ascending=False).reset_index(drop=True)
    out.insert(1, 'rank_within_cell_z', range(1, len(out) + 1))
    out.to_csv(OUT, index=False, encoding='utf-8-sig')

    print('单细胞总数 %d，细胞类型 %d 种' % (n_cells, n_types))
    print(out[['rank_within_cell_z', 'CellType', 'n_cells',
               'MRG4_score_within_cell_z', 'MRG4_score_type_level_z']].to_string(index=False))
    print()
    print('口径 A 前四：', ', '.join('%s %+.5f' % (r.CellType, r.MRG4_score_within_cell_z)
                                    for r in out.head(4).itertuples()))
    print('口径 A 末位：', ', '.join('%s %+.5f' % (r.CellType, r.MRG4_score_within_cell_z)
                                    for r in out.tail(3).itertuples()))
    print('OK ->', OUT)


if __name__ == '__main__':
    main()
