#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
检验 2（本地可跑部分）：空间细胞类型线索
========================================
审稿人："AHBA 是组织级表达，无法判断信号来自神经元/星形胶质/小胶质/内皮/少突。"
本脚本用已有全基因 Schaefer-400 表达矩阵，在 400 个脑区上检验
JAK-STAT 基因的空间表达模式与哪类细胞标志基因共变。

这不是单细胞证据，而是"区域级细胞组成线索"——可作为补充证据，
并诚实说明其局限（区域混杂、不能定位到细胞）。

数据：F:/炎症因子文章/expr_schaefer400_allgenes.csv
输出：outputs/test2/test2_celltype_spatial.csv 与 .md
"""
import sys, warnings
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, zscore

BASE = Path(r"F:/炎症因子文章")
OUT  = Path(r"C:/Users/ThinkPad/WorkBuddy/2026-10-20-11-06-57/outputs/test2")
OUT.mkdir(parents=True, exist_ok=True)

GENES13 = ['IL6R','IL6ST','JAK1','JAK2','JAK3','TYK2',
           'STAT1','STAT3','STAT4','STAT5A','STAT5B','STAT6','SOCS1']

PANELS = {
    'microglia':   ['AIF1','CSF1R','P2RY12','TMEM119','CX3CR1','ITGAM','CD68','TREM2'],
    'astrocyte':   ['GFAP','AQP4','ALDH1L1','SLC1A3','SLC1A2','GJA1','S100B'],
    'oligo':       ['MBP','PLP1','MOG','MAG','CNP','OLIG1','OLIG2','SOX10'],
    'opc':         ['PDGFRA','CSPG4','OLIG1'],
    'endothelial': ['CLDN5','PECAM1','FLT1','SLC2A1','MFSD2A','RAMP2'],
    'pericyte':    ['PDGFRB','RGS5','ACTA2'],
    'neuron_exc':  ['SLC17A7','SATB2','CUX2','BCL11B','RORB'],
    'neuron_inh':  ['GAD1','GAD2','SLC6A1','PVALB','SST','VIP'],
}

print("读取表达矩阵 ...")
expr = pd.read_csv(BASE / "expr_schaefer400_allgenes.csv", index_col=0)
expr.index = expr.index.astype(str)
n_raw = expr.shape[0]
valid = expr.notna().all(axis=1)
expr = expr.loc[valid]
print(f"  原始 {n_raw} parcels；AHBA 采样有效掩码后 {expr.shape[0]} parcels x {expr.shape[1]} genes"
      f"（每基因固定缺 {n_raw - int(valid.sum())} 个脑区）")

# 面板得分：标志基因的平均 z 分数（跨脑区标准化后平均）
panel_scores = {}
for name, genes in PANELS.items():
    gs = [g for g in genes if g in expr.columns]
    if len(gs) < 2:
        print(f"  [WARN] {name} 可用基因过少 ({gs})")
        continue
    Z = expr[gs].apply(zscore, axis=0)
    panel_scores[name] = Z.mean(axis=1).to_numpy(float)
    print(f"  {name:12s} n_markers={len(gs)}  {gs}")

P = pd.DataFrame(panel_scores)          # 400 x n_panels
print(f"\n面板得分矩阵: {P.shape}")
print("面板间相关（应低度共线，说明可区分）：")
print(P.corr().round(2).to_string())

# 13 基因 x 面板
rows = []
for g in GENES13:
    if g not in expr.columns:
        print(f"  [WARN] {g} 不在表达矩阵")
        continue
    x = expr[g].to_numpy(float)
    rec = {'gene': g}
    best, best_r = None, -2
    for name in P.columns:
        r, p = spearmanr(x, P[name].to_numpy(float))
        rec[f'r_{name}'] = r
        rec[f'p_{name}'] = p
        if r > best_r:
            best_r, best = r, name
    rec['best_panel'] = best
    rec['best_r'] = best_r
    # 神经元 vs 胶质/血管 差值
    neuro = np.nanmean([rec.get('r_neuron_exc', np.nan), rec.get('r_neuron_inh', np.nan)])
    gliov = np.nanmean([rec.get('r_microglia', np.nan), rec.get('r_astrocyte', np.nan),
                        rec.get('r_endothelial', np.nan), rec.get('r_oligo', np.nan)])
    rec['neuro_r_mean'] = neuro
    rec['gliovasc_r_mean'] = gliov
    rec['gliovasc_minus_neuro'] = gliov - neuro
    rows.append(rec)

df = pd.DataFrame(rows)
df.to_csv(OUT / "test2_celltype_spatial.csv", index=False)

print("\n" + "=" * 80)
print("13 基因 x 细胞类型面板（Spearman rho，400 脑区）")
print("=" * 80)
show = ['gene','best_panel','best_r','r_microglia','r_astrocyte','r_endothelial',
        'r_oligo','r_neuron_exc','r_neuron_inh','gliovasc_minus_neuro']
print(df[show].round(3).to_string(index=False))

print("\n" + "=" * 80)
print("关键判读")
print("=" * 80)
for g in ['IL6R','STAT3','STAT4','JAK1','IL6ST']:
    if g in set(df['gene']):
        r = df[df['gene'] == g].iloc[0]
        print(f"  {g:8s} 最高面板 = {r['best_panel']:12s} (rho={r['best_r']:+.3f});"
              f" 胶质/血管 - 神经元 = {r['gliovasc_minus_neuro']:+.3f}")

md = ["# 检验 2（区域级）：JAK-STAT 空间表达的细胞类型线索\n",
      "> 方法：全基因 Schaefer-400 表达矩阵，标志基因面板（文献标准 markers）经跨脑区 z 标准后取均值作为面板得分；",
      "> 13 个 JAK-STAT 基因与各面板做 Spearman 相关。**n = 400 脑区。**\n",
      "> **局限（必须写进论文）**：这是区域级共变，不是单细胞证据；面板得分混入区域间细胞组成差异、",
      "> 血管密度与皮层层次差异。真正的细胞来源仍需 snRNA-seq 或死后脑免疫组化。\n",
      "## 面板\n"]
for name, genes in PANELS.items():
    gs = [g for g in genes if g in expr.columns]
    md.append(f"- **{name}** (n={len(gs)}): {', '.join(gs)}")
md.append("\n## 结果\n")
md.append(df[show].round(3).to_markdown(index=False))
md.append("\n## 面板间相关\n")
md.append(P.corr().round(2).to_markdown())
(OUT / "test2_celltype_spatial.md").write_text("\n".join(md), encoding='utf-8')
print(f"\n输出 -> {OUT}")
