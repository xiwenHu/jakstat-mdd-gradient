#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
rho 符号诊断：IL6R 空间表达 × MDD 皮层厚度效应量
================================================
回答"正 rho 到底是什么意思"，并给出可写进论文的判读。

输出：outputs/sign_check/
"""
import sys, warnings
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
import pandas as pd
import nibabel as nib
from pathlib import Path
from scipy.stats import spearmanr, linregress

BASE = Path(r"F:/炎症因子文章")
OUT  = Path(r"C:/Users/ThinkPad/WorkBuddy/2026-10-20-11-06-57/outputs/sign_check")
OUT.mkdir(parents=True, exist_ok=True)

GENE = 'IL6R'
GENES13 = ['IL6R','IL6ST','JAK1','STAT3','STAT4','JAK2','STAT1']

# ---- 表型 ----
d = pd.read_csv(BASE / "traits_schaefer400/trait_mddadult_case-controls_CortThick.csv")
col = [c for c in d.columns if c != 'parcel'][0]
d = d.copy(); d['parcel'] = d['parcel'].astype(int)
print(f"表型列 = {col}")
v = d[col].to_numpy(float)
print(f"  n={len(v)}  min={v.min():+.4f}  max={v.max():+.4f}  mean={v.mean():+.4f}")
print(f"  负值比例 = {(v<0).mean():.3f}   正值比例 = {(v>0).mean():.3f}")
print(f"  >>> 负值 = 皮层更薄（MDD 相对对照）；正值 = 相对保留")
print(f"  10% / 50% / 90% 分位: {np.percentile(v,10):+.4f} / {np.percentile(v,50):+.4f} / {np.percentile(v,90):+.4f}")

# ---- 表达 ----
expr = pd.read_csv(BASE / "expr_schaefer400_allgenes.csv", index_col=0)
expr.index = expr.index.astype(int)
valid = expr.notna().all(axis=1)
expr = expr.loc[valid]
print(f"\n表达矩阵有效脑区 = {expr.shape[0]}")

# ---- 分区名（Schaefer-400 7 网络）----
annot = BASE / "optimized_figures/schaefer_annot/lh.Schaefer2018_400Parcels_7Networks_order.annot"
names = None
try:
    labels, ctab, nm = nib.freesurfer.read_annot(str(annot))
    names = [x.decode() if isinstance(x, bytes) else str(x) for x in nm]
    print(f"annot 分区名 n={len(names)}  例: {names[:3]}")
except Exception as e:
    print(f"[WARN] 读分区名失败: {e}")

def pname(p):
    if names is not None and 0 <= p - 1 < len(names):
        return names[p - 1]
    return f"parcel_{p}"

df = d[d['parcel'].isin(expr.index)].sort_values('parcel').reset_index(drop=True)
df['IL6R'] = expr.loc[df['parcel'], GENE].to_numpy(float)
df['name'] = df['parcel'].apply(pname)

r, p = spearmanr(df['IL6R'], df[col])
lr = linregress(df['IL6R'], df[col])
print("\n" + "=" * 78)
print(f"{GENE} 表达 × 皮层厚度 d（n = {len(df)} 脑区）")
print("=" * 78)
print(f"  Spearman rho = {r:+.4f}   P = {p:.3e}")
print(f"  线性斜率     = {lr.slope:+.4f} d/单位表达   截距 = {lr.intercept:+.4f}")
print(f"  >>> 斜率为{'正' if lr.slope>0 else '负'}：表达越高，d 越{'接近 0 / 越正' if lr.slope>0 else '负'}")
print(f"      即表达越高的脑区，皮层变薄越{'轻（相对保留）' if lr.slope>0 else '重（更易受损）'}")

# |d| 版本
r_abs, p_abs = spearmanr(df['IL6R'], df[col].abs())
print(f"\n  若改用 |d|（变薄幅度）: rho = {r_abs:+.4f}  P = {p_abs:.3e}")
print(f"  >>> |d| 版本的 rho 符号 = {'正' if r_abs>0 else '负'}，"
      f"{'与"高表达→更易变薄"叙事一致' if r_abs>0 else '仍不支持"高表达→更易变薄"'}")

# 分位数组
df['bin'] = pd.qcut(df['IL6R'], 5, labels=['Q1最低','Q2','Q3','Q4','Q5最高'])
print("\n" + "=" * 78)
print("按 IL6R 表达五分位看皮层厚度 d")
print("=" * 78)
g = df.groupby('bin', observed=True)[col].agg(['count','mean','median','min','max'])
g.columns = ['n脑区','d均值','d中位','d最小','d最大']
print(g.round(4).to_string())
print(f"\n  Q5最高表达组 d 均值 = {df[df['bin']=='Q5最高'][col].mean():+.4f}")
print(f"  Q1最低表达组 d 均值 = {df[df['bin']=='Q1最低'][col].mean():+.4f}")
print(f"  差值 = {df[df['bin']=='Q5最高'][col].mean() - df[df['bin']=='Q1最低'][col].mean():+.4f}")
print(f"  >>> 正差值 = 高表达区变薄更轻；负差值 = 高表达区变薄更重")

# 极端脑区
print("\n" + "=" * 78)
print("IL6R 表达最高的 10 个脑区 vs 最低的 10 个")
print("=" * 78)
hi = df.nlargest(10, 'IL6R')[['parcel','name','IL6R',col]]
lo = df.nsmallest(10, 'IL6R')[['parcel','name','IL6R',col]]
print("  【表达最高 10】")
print(hi.round(4).to_string(index=False))
print(f"  -> d 均值 = {hi[col].mean():+.4f}  （负=变薄，正=保留）")
print("\n  【表达最低 10】")
print(lo.round(4).to_string(index=False))
print(f"  -> d 均值 = {lo[col].mean():+.4f}")

# 13 基因一起
print("\n" + "=" * 78)
print("多个 JAK-STAT 基因的符号（正 = 高表达区变薄轻）")
print("=" * 78)
rows = []
for gn in GENES13:
    if gn not in expr.columns:
        continue
    x = expr.loc[df['parcel'], gn].to_numpy(float)
    rr, pp = spearmanr(x, df[col])
    rr2, _ = spearmanr(x, df[col].abs())
    rows.append(dict(gene=gn, rho_vs_d=rr, p=pp, rho_vs_absd=rr2))
    print(f"  {gn:8s} rho(d) = {rr:+.4f}   rho(|d|) = {rr2:+.4f}")
pd.DataFrame(rows).to_csv(OUT / "sign_check_13genes.csv", index=False)
df.to_csv(OUT / "sign_check_per_parcel.csv", index=False)

print("\n" + "=" * 78)
print("判读")
print("=" * 78)
print("""  若 rho(d) > 0 且 rho(|d|) < 0：
      表达越高 → 变薄越轻。这是"相对保留"方向，与"JAK-STAT 驱动萎缩"相反。
      可能的正确解读：(a) 代偿/保护；(b) 该区域本就少受累；(c) 表达梯度反映
      皮层层次/细胞组成而非疾病过程。
  若 rho(d) < 0 且 rho(|d|) > 0：
      表达越高 → 变薄越重。这才支持损伤叙事。
  无论哪种，都必须在论文里明写"方向"，不能只写"空间汇聚"。""")
print(f"\n输出 -> {OUT}")
