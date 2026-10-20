#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
检验 4B'/4C'：空间结果的匹配随机基因集零分布
=============================================
在拿不到球面坐标（无法跑 spin）的情况下，用"匹配随机基因集"作为零分布，
回答两个同样关键、且不依赖几何的问题：

  4B'  JAK-STAT 13 基因作为"基因集"，与 MDD 皮层厚度的空间对齐
       是否强于随机 13 基因集？（主表型）
  4C'  22 表型 x 13 基因 = 286 个相关构成的"网格"，其极端值强度
       是否强于随机 13 基因集在同一批表型上的表现？

关键设计：零分布与观测使用**同一批表型图**，因此表型的空间自相关
在两侧完全相同，比较是公平的。未控制的是基因表达图自身的空间自相关
（spin 版可控制，需 fsaverage5 球面坐标，另附脚本）。

匹配变量：表达量十分位 x 表达图标准差十分位
输出：outputs/test4/test4b_spatial_geneset_null.csv
"""
import sys, time, warnings
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import rankdata

BASE = Path(r"F:/炎症因子文章")
OUT  = Path(r"C:/Users/ThinkPad/WorkBuddy/2026-10-20-11-06-57/outputs/test4")
OUT.mkdir(parents=True, exist_ok=True)

GENES13 = ['IL6R','IL6ST','JAK1','JAK2','JAK3','TYK2',
           'STAT1','STAT3','STAT4','STAT5A','STAT5B','STAT6','SOCS1']
N_SET   = 2000
RS      = 42

def ranknorm_cols(M):
    R = np.apply_along_axis(rankdata, 0, M)
    return (R - R.mean(0)) / R.std(0)

print("读取表达矩阵 ...")
expr = pd.read_csv(BASE / "expr_schaefer400_allgenes.csv", index_col=0)
expr.index = expr.index.astype(str)
n_raw = expr.shape[0]
valid = expr.notna().all(axis=1)
expr = expr.loc[valid]
expr = expr.loc[sorted(expr.index, key=lambda s: int(s))]
n_par, n_genes = expr.shape
print(f"  原始 {n_raw} parcels；有效掩码后 {n_par} parcels x {n_genes} genes")

print("读取 22 个表型 ...")
valid_parcels = set(expr.index)
traits = {}
for f in sorted((BASE / "traits_schaefer400").glob("trait_*.csv")):
    d = pd.read_csv(f)
    col = [c for c in d.columns if c != 'parcel'][0]
    d2 = d[d['parcel'].astype(str).isin(valid_parcels)].sort_values('parcel')
    vec = d2[col].to_numpy(float)
    if len(vec) != n_par:
        raise SystemExit(f"[FAIL] {f.name} 有效脑区数 {len(vec)} != {n_par}")
    traits[f.stem.replace('trait_', '')] = vec
tkeys = list(traits.keys())
Y = np.column_stack([traits[k] for k in tkeys])
Yn = ranknorm_cols(Y)                              # n_par x 22
print(f"  {len(tkeys)} 个表型")

X  = expr.to_numpy(float)                          # n_par x n_genes
Xn = ranknorm_cols(X)

t0 = time.time()
Rho_all = (Xn.T @ Yn) / (n_par - 1)                # n_genes x 22
print(f"  全基因 x 全表型 相关矩阵算完 ({time.time()-t0:.1f}s), shape={Rho_all.shape}")

idx13 = [expr.columns.get_loc(g) for g in GENES13 if g in expr.columns]
g13   = [g for g in GENES13 if g in expr.columns]
print(f"  13 基因命中 {len(g13)} 个")

primary_i = tkeys.index('mddadult_case-controls_CortThick')

# ---- 匹配变量 ----
mean_expr = X.mean(axis=0)
sd_expr   = X.std(axis=0)
def bins(v, k=10):
    return pd.qcut(pd.Series(v).rank(method='first'), k, labels=False, duplicates='drop').to_numpy()
b1, b2 = bins(mean_expr), bins(sd_expr)
key = b1.astype(str) + '_' + b2.astype(str)
keys_obs = key[idx13]

pool_idx = {k: np.where(key == k)[0] for k in np.unique(keys_obs)}
print("  匹配格基因数: " + ", ".join(f"{k}={len(pool_idx[k])}" for k in sorted(pool_idx)))

rng = np.random.default_rng(RS)

def draw_set():
    out = []
    for k in keys_obs:
        p = pool_idx[k]
        out.append(p[rng.integers(len(p))])
    return np.array(out, dtype=int)

# ---- 4B' 基因集级（主表型）----
abs_obs_13 = np.abs(Rho_all[idx13, primary_i])
T_obs_B = float(abs_obs_13.mean())
print("\n" + "=" * 70)
print("4B'  基因集级（主表型 mddadult_case-controls_CortThick）")
print("=" * 70)
print("  13 基因 rho: " + ", ".join(f"{g}={r:+.3f}"
      for g, r in zip(g13, Rho_all[idx13, primary_i])))
print(f"  T_obs = mean|rho| = {T_obs_B:.4f}")

T_null_B = np.empty(N_SET)
for s in range(N_SET):
    T_null_B[s] = np.abs(Rho_all[draw_set(), primary_i]).mean()
p_B = (np.sum(T_null_B >= T_obs_B) + 1) / (N_SET + 1)
print(f"  零分布 = {T_null_B.mean():.4f} ± {T_null_B.std():.4f}")
print(f"  >>> 经验 P = {p_B:.4f}")

# ---- 4C' 网格级 ----
abs_all = np.abs(Rho_all)
obs_grid = abs_all[idx13, :].ravel()
K = 18                                   # DK68 中 p_spin<0.05 的个数
T_top_obs = float(np.sort(obs_grid)[::-1][:K].sum())
T_mean_obs = float(obs_grid.mean())

print("\n" + "=" * 70)
print(f"4C'  网格级（22 表型 x 13 基因 = 286）")
print("=" * 70)
print(f"  观测 |rho|: mean = {T_mean_obs:.4f}  max = {obs_grid.max():.4f}")
print(f"  统计量: top-{K} sum|rho| = {T_top_obs:.4f}   mean|rho| = {T_mean_obs:.4f}")

T_null_top = np.empty(N_SET)
T_null_mean = np.empty(N_SET)
t0 = time.time()
for s in range(N_SET):
    g = abs_all[draw_set(), :].ravel()
    T_null_mean[s] = g.mean()
    T_null_top[s]  = np.sort(g)[::-1][:K].sum()
    if (s + 1) % 500 == 0:
        print(f"    {s+1}/{N_SET} ({time.time()-t0:.1f}s)")
p_C_top  = (np.sum(T_null_top  >= T_top_obs)  + 1) / (N_SET + 1)
p_C_mean = (np.sum(T_null_mean >= T_mean_obs) + 1) / (N_SET + 1)
print(f"  top-{K} : 零分布 = {T_null_top.mean():.4f} ± {T_null_top.std():.4f}   经验 P = {p_C_top:.4f}")
print(f"  mean   : 零分布 = {T_null_mean.mean():.4f} ± {T_null_mean.std():.4f}   经验 P = {p_C_mean:.4f}")

pd.DataFrame([
    dict(part="4B' 基因集 mean|rho| (主表型)", observed=T_obs_B,
         null_mean=T_null_B.mean(), null_sd=T_null_B.std(), emp_p=p_B),
    dict(part=f"4C' 网格 top-{K} sum|rho|", observed=T_top_obs,
         null_mean=T_null_top.mean(), null_sd=T_null_top.std(), emp_p=p_C_top),
    dict(part="4C' 网格 mean|rho|", observed=T_mean_obs,
         null_mean=T_null_mean.mean(), null_sd=T_null_mean.std(), emp_p=p_C_mean),
]).to_csv(OUT / "test4b_spatial_geneset_null.csv", index=False)

# 每基因在 22 表型上的最大 |rho|（供表中定位用）
rows = []
for g, i in zip(g13, idx13):
    a = abs_all[i, :]
    rows.append(dict(gene=g, rho_primary=Rho_all[i, primary_i],
                     absrho_primary=a[primary_i],
                     max_absrho=a.max(), best_trait=tkeys[int(a.argmax())],
                     mean_absrho_22=a.mean()))
pd.DataFrame(rows).to_csv(OUT / "test4b_per_gene_grid.csv", index=False)

print("\n" + "=" * 70)
print("汇总")
print("=" * 70)
print(f"  4B'  基因集 mean|rho| = {T_obs_B:.4f}, 经验 P = {p_B:.4f}")
print(f"  4C'  网格 top-{K} sum|rho| 经验 P = {p_C_top:.4f}")
print(f"  4C'  网格 mean|rho| 经验 P = {p_C_mean:.4f}")
print(f"\n  输出 -> {OUT}")
