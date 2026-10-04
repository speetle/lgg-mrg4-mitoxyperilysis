# -*- coding: utf-8 -*-
"""R3-M6 关闭：免疫特征表可由沉积对象重建
用与主模型相同的 cBioPortal z 谱（median_all_sample_Zscores）重算 11 个免疫特征，
沉积标志基因清单 + 逐患者分值。
"""
import json, numpy as np, pandas as pd, urllib.request
from scipy.stats import mannwhitneyu

B = 'https://www.cbioportal.org/api'
def post(u, p):
    r = urllib.request.Request(u, data=json.dumps(p).encode(),
        headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
    return json.load(urllib.request.urlopen(r, timeout=300))

PROF_Z = 'lgg_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_all_sample_Zscores'
SL = 'lgg_tcga_pan_can_atlas_2018_rna_seq_v2_mrna'

PANEL = {
    'CD8+ T cells':        ['CD8A', 'CD8B', 'GZMA', 'GZMK', 'PRF1', 'NKG7'],
    'CD4+ T cells':        ['CD4', 'IL7R', 'CD40LG', 'ICOS', 'CD28', 'TRAT1'],
    'Regulatory T cells':  ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TNFRSF18'],
    'NK cells':            ['KLRD1', 'KLRF1', 'NCR1', 'GNLY', 'KIR2DL3'],
    'B cells':             ['MS4A1', 'CD79A', 'CD79B', 'CD19'],
    'M1 macrophages':      ['NOS2', 'IL1B', 'TNF', 'CXCL9', 'CXCL11'],
    'M2 macrophages':      ['CD163', 'MRC1', 'MSR1', 'IL10', 'TGFB1'],
    'Dendritic cells':     ['ITGAX', 'CLEC9A', 'LAMP3', 'CD1C', 'FCER1A'],
    'Neutrophils':         ['FCGR3B', 'CSF3R', 'S100A8', 'S100A9', 'CXCR2'],
    'Immune checkpoints':  ['PDCD1', 'CD274', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT'],
    'Interferon signature': ['IFNG', 'STAT1', 'CXCL10', 'IDO1', 'HLA-DRA'],
}
allg = sorted({g for v in PANEL.values() for g in v})
resolved = post(B + '/genes/fetch?geneIdType=HUGO_GENE_SYMBOL', allg)
ok = {x['hugoGeneSymbol'] for x in resolved}
unresolved = [g for g in allg if g not in ok]
print('标志基因 %d 个，UniProt/HUGO 未解析：%s' % (len(allg), unresolved or '无'))
missing = {k: [g for g in v if g not in ok] for k, v in PANEL.items()}
missing = {k: v for k, v in missing.items() if v}
if missing: print('各特征缺失基因：', missing)
inv = {x['entrezGeneId']: x['hugoGeneSymbol'] for x in resolved if x['hugoGeneSymbol'] in ok}
ids = sorted(inv)
rows = post(f'{B}/molecular-profiles/{PROF_Z}/molecular-data/fetch?projection=SUMMARY',
            {'entrezGeneIds': ids, 'sampleListId': SL})
Z = pd.DataFrame([{'upk': x.get('uniquePatientKey'), 'gene': inv[x['entrezGeneId']], 'v': x['value']}
                  for x in rows]).pivot_table(index='upk', columns='gene', values='v')
print('z 谱矩阵', Z.shape)

# ---------- 逐患者免疫特征分值 ----------
per = pd.DataFrame(index=Z.index)
for feat, gl in PANEL.items():
    gl = [g for g in gl if g in Z.columns]
    per[feat] = Z[gl].mean(axis=1)
per.to_csv('Table_S26_immune_features_per_patient.csv', encoding='utf-8-sig')

# ---------- 标志基因清单 ----------
mrows = []
for feat, gl in PANEL.items():
    for g in gl:
        mrows.append([feat, g, 'included' if g in Z.columns else 'not resolved/deposited'])
pd.DataFrame(mrows, columns=['Immune_feature', 'Marker_gene', 'Status']).to_csv(
    'Table_S26b_immune_marker_gene_panel.csv', index=False, encoding='utf-8-sig')
print('-> Table_S26 / S26b')

# ---------- 分组比较（Model A 中位分割，n = 510） ----------
G = pd.read_csv('lgg_g1_dataset.csv', encoding='utf-8-sig', low_memory=False)
d = G.merge(per, left_on='pat', right_index=True, how='inner')
print('可评估患者 %d' % len(d))
d['grp'] = np.where(d['risk'] > d['risk'].median(), 'high', 'low')
out = []
for feat in PANEL:
    lo = d[d.grp == 'low'][feat].values; hi = d[d.grp == 'high'][feat].values
    u, p = mannwhitneyu(hi, lo, alternative='two-sided')
    n1, n2 = len(hi), len(lo)
    rb = 2 * u / (n1 * n2) - 1
    out.append([feat, len([g for g in PANEL[feat] if g in Z.columns]),
                round(float(lo.mean()), 4), round(float(hi.mean()), 4),
                round(float(hi.mean() - lo.mean()), 4), float(p), round(float(rb), 3)])
tab = pd.DataFrame(out, columns=['Immune_feature', 'n_marker_genes', 'low_risk_mean_z',
                                 'high_risk_mean_z', 'difference_high_minus_low', 'P_MannWhitney', 'rank_biserial_r'])
tab.to_csv('Table_S26c_immune_features_by_group.csv', index=False, encoding='utf-8-sig')
print(tab.to_string(index=False))
n_up = int((tab['difference_high_minus_low'] > 0).sum())
n_sig = int((tab['P_MannWhitney'] < 0.05).sum())
print('\n高表达方向：%d/%d ；P<0.05：%d/%d' % (n_up, len(tab), n_sig, len(tab)))
json.dump({'n_patients': int(len(d)), 'n_features': len(tab), 'n_higher_in_high_risk': n_up,
           'n_significant': n_sig, 'unresolved_genes': unresolved,
           'difference_range': [float(tab['difference_high_minus_low'].min()),
                                float(tab['difference_high_minus_low'].max())]},
          open('immune_features_result.json', 'w'), ensure_ascii=False, indent=1)
print('-> immune_features_result.json')
