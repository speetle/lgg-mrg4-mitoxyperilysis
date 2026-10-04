# -*- coding: utf-8 -*-
"""Figure 16 + 补充表 S10：SL040（GBM 单细胞）中 4 个 MRG 的细胞类型定位。"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

df = pd.read_csv('/tmp/mitoxy/sl040_expr.csv', index_col=0)
GENES = ["BMP1", "KIF15", "TRAF3", "PDGFA"]
REF = "MKI67"

ORDER = ["Tumor", "Astrocytes", "OPCs", "Oligodendrocytes", "Endothelial", "Pericytes",
         "VSMC", "Fibroblasts", "Microglia", "MgTAM", "MoTAM", "Macrophages",
         "Monocytes", "NK, T-cells", "B-cells"]
ORDER = [c for c in ORDER if c in set(df['CellType'])]
cols = GENES + [REF]

mean = df.groupby('CellType')[cols].mean().loc[ORDER]
pct = df.groupby('CellType')[cols].apply(lambda d: (d > 0).mean() * 100).loc[ORDER]
cnt = df['CellType'].value_counts().loc[ORDER]

# 补充表 S10
S10 = pd.concat({'n_cells': cnt, 'mean': mean, 'pct_expressing': pct}, axis=1)
S10.to_csv('/tmp/mitoxy/Table_S10_singlecell_celltype.csv')
print(S10.round(3).to_string())

# 恶性细胞状态
t = df[df['CellType'] == 'Tumor']
STATE_ORDER = ["NPC", "OPC", "AC", "MES"]
sm = t.groupby('NeftelClass')[cols].mean().reindex(STATE_ORDER)
sp = t.groupby('NeftelClass')[cols].apply(lambda d: (d > 0).mean() * 100).reindex(STATE_ORDER)
sc = t['NeftelClass'].value_counts().reindex(STATE_ORDER)
ST = pd.concat({'n_cells': sc, 'mean': sm, 'pct_expressing': sp}, axis=1)
ST.to_csv('/tmp/mitoxy/Table_S10b_tumor_states.csv')
print("\n", ST.round(3).to_string())

# ---------------- Figure 16 ----------------
cmap = LinearSegmentedColormap.from_list('gyr', ['#F2F4F7', '#FDE2E2', '#F5A3A3', '#D7301F', '#7A1004'])
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6),
                         gridspec_kw={'width_ratios': [2.55, 1.0]})

def dot(ax, M, P, order, labels, ref=None, rot=0, fs=8.6):
    for i, ct in enumerate(order):
        for j, g in enumerate(cols):
            r = 4 + 118 * np.sqrt(P.iloc[i, j] / 100.0)
            ax.scatter(j, i, s=r, c=[cmap(min(M.iloc[i, j] / 4.0, 1.0))],
                       edgecolors='#4A4A4A', linewidths=0.35, zorder=3)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, fontsize=fs, style='italic',
                       rotation=rot, ha='right' if rot else 'center')
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(labels, fontsize=8.2)
    ax.invert_yaxis()
    ax.set_ylim(len(order) - 0.4, -0.6)
    ax.set_xlim(-0.6, len(cols) - 0.4)
    for j in range(1, len(cols)):
        ax.axvline(j - 0.5, color='#D9D9D9', lw=0.5, zorder=0)
    ax.grid(axis='y', alpha=0.12, lw=0.5)
    ax.set_axisbelow(True)
    for s_ in ('top', 'right'):
        ax.spines[s_].set_visible(False)

dot(axes[0], mean, pct, ORDER, ['%s  (n=%s)' % (c, format(cnt[c], ',')) for c in ORDER])
axes[0].set_title('A  All cells, by cell type', fontsize=9.4, loc='left')
axes[0].set_ylabel('Cell type', fontsize=8.8)

dot(axes[1], sm, sp, STATE_ORDER, ['%s  (n=%s)' % (s, format(sc[s], ',')) for s in STATE_ORDER],
    rot=45, fs=7.6)
axes[1].set_title('B  Malignant cells, by state', fontsize=9.4, loc='left')

# 图例：大小 + 颜色
h = [plt.scatter([], [], s=4 + 118 * np.sqrt(p / 100.0), c='#C9C9C9',
                 edgecolors='#4A4A4A', linewidths=0.35, label='%d%%' % p) for p in (10, 25, 50, 75)]
leg1 = axes[1].legend(handles=h, title='% expressing', frameon=False, fontsize=6.8,
                      title_fontsize=7.0, loc='upper left', bbox_to_anchor=(1.85, 0.98),
                      labelspacing=0.85, handletextpad=1.0, borderpad=0.2)

sm_ = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 4))
cb = fig.colorbar(sm_, ax=axes[1], fraction=0.04, pad=0.22)
cb.set_label('mean log$_1$p expression', fontsize=7.2)
cb.ax.tick_params(labelsize=6.8)

fig.text(0.5, 0.005,
         'SL040 glioblastoma, CZ CELLxGENE Discover (135,482 cells); 11 genes, log$_1$p-normalised. '
         'MKI67 shown as a proliferation reference. Cell-type annotation validated by markers: '
         'CD68 (macrophage), PECAM1 (endothelial), MOG (oligodendrocyte), GFAP (astrocyte/tumour).',
         ha='center', va='bottom', fontsize=6.4, color='#555555')
fig.tight_layout(rect=[0, 0.045, 1, 1])
for ext in ('png', 'pdf'):
    fig.savefig('/tmp/mitoxy/Fig16_单细胞定位.%s' % ext, dpi=300, bbox_inches='tight')
print("\nOK  Fig16_单细胞定位.png/.pdf")
