#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
决定性检验：空间"汇聚"是否只是皮层感觉-联合轴（S-A axis）的假象
================================================================
假说：AHBA 基因表达的主梯度（PC1）与 MDD 皮层变薄的分布都沿感觉→联合轴，
因此几乎任何有梯度的基因都会与效应量图"显著相关"——这解释了
  (a) 随机基因集 P=0.98；(b) rho 为正（高表达区=感觉区=相对保留）；
  (c) 分区依赖；(d) STAT3/STAT4 等各自"显著"。

检验：控制 PC1 后，IL6R 与 d 的偏相关是否塌陷。
"""
import sys, warnings
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, pearsonr, rankdata

BASE = Path(r"F:/炎症因子文章")
OUT  = Path(r"C:/Users/ThinkPad/WorkBuddy/2026-10-20-11-06-57/outputs/sign_check")
OUT.mkdir(parents=True, exist_ok=True)

GENES = ['IL6R','IL6ST','JAK1','JAK2','STAT3','STAT4','STAT1','TYK2','STAT5A','STAT5B','STAT6','JAK3','SOCS1']

d = pd.read_csv(BASE / "traits_schaefer400/trait_mddadult_case-controls_CortThick.csv")
col = [c for c in d.columns if c != 'parcel'][0]
d['parcel'] = d['parcel'].astype(int)

expr = pd.read_csv(BASE / "expr_schaefer400_allgenes.csv", index_col=0)
expr.index = expr.index.astype(int)
expr = expr.loc[expr.notna().all(axis=1)]
expr = expr.loc[sorted(expr.index)]
common = sorted(set(expr.index) & set(d['parcel']))
expr = expr.loc[common]
dv = d.set_index('parcel').loc[common, col].to_numpy(float)
print(f"有效脑区 n = {len(common)}")

# ---- PC1（S-A 轴代理）----
X = expr.to_numpy(float)
Xz = (X - X.mean(0)) / X.std(0)
U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
pc1 = U[:, 0] * S[0]
print(f"PC1 解释方差比 = {(S[0]**2 / (S**2).sum()):.3f}  （AHBA 典型的第一主梯度）")
pc1 /= pc1.std()

# PC1 与 d
r_pc1_d, p_pc1_d = spearmanr(pc1, dv)
print(f"\nPC1 × d: rho = {r_pc1_d:+.4f}   P = {p_pc1_d:.2e}")

def partial_spearman(x, y, z):
    """控制 z 的 Spearman 偏相关（秩变换后做残差相关）。"""
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def resid(a, b):
        b1 = np.vstack([b, np.ones_like(b)]).T
        beta, *_ = np.linalg.lstsq(b1, a, rcond=None)
        return a - b1 @ beta
    ex = resid(rx, rz); ey = resid(ry, rz)
    return pearsonr(ex, ey)

print("\n" + "=" * 84)
print("JAK-STAT 基因 × d：控制 PC1 前后")
print("=" * 84)
print(f"  {'gene':8s} {'rho(d)':>9s} {'P':>10s} {'rho|PC1':>9s} {'偏相关':>9s} {'P_partial':>11s} {'衰减%':>8s}")
rows = []
for g in GENES:
    if g not in expr.columns:
        continue
    x = expr[g].to_numpy(float)
    r0, p0 = spearmanr(x, dv)
    r_pc, _ = spearmanr(x, pc1)
    rp, pp = partial_spearman(x, dv, pc1)
    drop = (1 - abs(rp) / abs(r0)) * 100 if abs(r0) > 1e-9 else np.nan
    rows.append(dict(gene=g, rho_d=r0, p_d=p0, rho_pc1=r_pc,
                     rho_partial=rp, p_partial=pp, pct_drop=drop))
    print(f"  {g:8s} {r0:+9.4f} {p0:10.2e} {r_pc:+9.4f} {rp:+9.4f} {pp:11.2e} {drop:7.1f}%")
pd.DataFrame(rows).to_csv(OUT / "partial_control_PC1.csv", index=False)

# ---- 再用 PC1-PC5 控制 ----
print("\n" + "=" * 84)
print("控制 PC1–PC5（多重主梯度）后的偏相关")
print("=" * 84)
pcs = (U[:, :5] * S[:5]) / np.std(U[:, :5] * S[:5], axis=0)

def partial_multi(x, y, Z):
    rx, ry = rankdata(x), rankdata(y)
    Zr = np.apply_along_axis(rankdata, 0, Z)
    A = np.hstack([Zr, np.ones((len(ry), 1))])
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None)
    by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    ex = rx - A @ bx; ey = ry - A @ by
    return pearsonr(ex, ey)

print(f"  {'gene':8s} {'rho(d)':>9s} {'偏相关(5PC)':>13s} {'P':>11s} {'衰减%':>8s}")
for r in rows:
    g = r['gene']
    x = expr[g].to_numpy(float)
    rp, pp = partial_multi(x, dv, pcs)
    drop = (1 - abs(rp) / abs(r['rho_d'])) * 100 if abs(r['rho_d']) > 1e-9 else np.nan
    print(f"  {g:8s} {r['rho_d']:+9.4f} {rp:+13.4f} {pp:11.2e} {drop:7.1f}%")

# ---- 反向：PC1 与 d 的关系有多强 ----
print("\n" + "=" * 84)
print("主梯度本身与 MDD 皮层厚度效应量的关系")
print("=" * 84)
for k in range(5):
    r, p = spearmanr(pcs[:, k], dv)
    print(f"  PC{k+1} × d: rho = {r:+.4f}   P = {p:.2e}")
print(f"\n  >>> 若 PC1 与 d 的相关接近 IL6R 与 d 的相关（≈0.3），")
print(f"      则 IL6R 的'汇聚'几乎可完全归因于主梯度，而非基因特异性。")

# 五分位
q = pd.qcut(pc1, 5, labels=['PC1-Q1感觉端','Q2','Q3','Q4','PC1-Q5联合端'])
print("\n按 PC1 五分位看 d：")
for lab in q.cat.categories:
    m = dv[q == lab]
    print(f"  {str(lab):14s} n={len(m):3d}  d均值={m.mean():+.4f}  d中位={np.median(m):+.4f}")

print(f"\n输出 -> {OUT}")
