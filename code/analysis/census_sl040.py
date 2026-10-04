# -*- coding: utf-8 -*-
"""SL040（GBM，CELLxGENE Discover，doi 10.1101/2023.09.01.555882）单细胞定位。
检查 BMP1/KIF15/TRAF3/PDGFA 在恶性细胞 vs 免疫/基质细胞中的表达。
数据：CZ CELLxGENE Discover 官方 h5ad（var_names 为 Ensembl ID，符号在 feature_name）。
"""
import numpy as np, pandas as pd
import anndata as ad

H5 = "/tmp/mitoxy/SL040.h5ad"
GENES = ["BMP1", "KIF15", "TRAF3", "PDGFA", "MKI67"]
MARK = ["PTPRC", "CD68", "GFAP", "MOG", "PECAM1", "TOP2A"]

A = ad.read_h5ad(H5, backed='r')
sym = pd.Series(A.var['feature_name'].values, index=np.arange(A.shape[1]))
want = GENES + MARK
pos = {}
for g in want:
    hit = np.where(sym.values == g)[0]
    if len(hit):
        pos[g] = int(hit[0])
missing = [g for g in want if g not in pos]
print("命中:", list(pos.keys()))
print("缺失:", missing)

idx = sorted(pos.values())
names = [sym.values[i] for i in idx]
sub = A[:, idx].to_memory()
X = sub.X
X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
print("矩阵:", X.shape)

df = pd.DataFrame(X, columns=names, index=A.obs_names)
for c in ["CellType", "Zone", "NeftelClass", "donor_id", "sample_id", "assay", "total_genes"]:
    if c in A.obs.columns:
        df[c] = A.obs[c].values
df.to_csv("/tmp/mitoxy/sl040_expr.csv")
print("写出 sl040_expr.csv:", df.shape)

ok = [g for g in GENES + MARK if g in df.columns]
print("\n=== 各细胞类型平均表达（log1p 归一化值）===")
print(df.groupby("CellType")[ok].mean().round(3).to_string())
print("\n=== 表达细胞比例 (%) ===")
print(df.groupby("CellType")[ok].apply(lambda d: (d > 0).mean() * 100).round(1).to_string())
print("\n=== 各细胞类型细胞数 ===")
print(df["CellType"].value_counts().to_string())
print("\n=== 恶性细胞按 Neftel 状态 ===")
t = df[df["CellType"] == "Tumor"]
print(t.groupby("NeftelClass")[ok].mean().round(3).to_string())
print("\n=== 恶性细胞按 Zone ===")
print(t.groupby("Zone")[ok].mean().round(3).to_string())
