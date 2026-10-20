#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
检验 4：正确零分布（替代"4.6 倍 O/E"与"15/286 稳健性"两处推断）
================================================================
Part 0  表格核对：论文正文数字 vs 实际结果文件
Part A  MAGMA 随机基因集零分布（表达量 + 基因长度 + SNP 数三匹配，10,000 次）
Part B  空间基因集 spin 检验：13 基因 mean|rho| 与 spin 零分布比较
Part C  网格级 spin 检验：286 组合的"极端个数/强度"与 spin 零分布比较

运行环境：C:/Users/ThinkPad/.workbuddy/binaries/python/envs/default
运行方式：python test4_null_distributions.py
"""
import os, sys, time, warnings
warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, rankdata, norm

BASE = Path(r"F:/炎症因子文章")
OUT  = Path(r"C:/Users/ThinkPad/WorkBuddy/2026-10-20-11-06-57/outputs/test4")
OUT.mkdir(parents=True, exist_ok=True)

GENES13 = ['IL6R','IL6ST','JAK1','JAK2','JAK3','TYK2',
           'STAT1','STAT3','STAT4','STAT5A','STAT5B','STAT6','SOCS1']
N_DRAW  = 10000
N_SPIN  = 2000
RS      = 42

def hdr(t):
    print("\n" + "=" * 70)
    print(t)
    print("=" * 70)

# ============================================================================
# Part 0 — 表格核对
# ============================================================================
hdr("Part 0  论文数字 vs 实际结果文件核对")

dk = pd.read_csv(BASE / "results_dk68/dk68_results.csv")
sc = pd.read_csv(BASE / "results_multi_trait/multi_trait_results.csv")

print(f"DK68   组合数 = {len(dk)}  (22 traits x 13 genes = 286)")
print(f"Schaefer 组合数 = {len(sc)}")
for name, df in [("DK68", dk), ("Schaefer400", sc)]:
    n_p = int((df['p_spin'] < 0.05).sum())
    n_q = int((df['q_spin'] < 0.05).sum())
    print(f"  {name:12s}: p_spin<0.05 -> {n_p:3d}   q_spin<0.05 -> {n_q:3d}")
print(f"  论文正文声称: 15 (5.2%)")
print(f"  >>> 结论：三者互不一致，须先确定稿件数字取自哪一次运行")

# 主表型（成人 MDD 全样本 CT）逐基因核对 Table 1
primary = dk[dk['trait'] == 'mddadult_case-controls_CortThick'].copy()
order = {g: i for i, g in enumerate(GENES13)}
primary['ord'] = primary['gene'].map(order)
primary = primary.sort_values('ord')
print("\n主表型 mddadult_case-controls_CortThick（DK68）逐基因：")
print(f"  {'gene':8s} {'rho':>8s} {'p_spin':>10s} {'q_spin':>10s}   论文 Table 1")
manu = {'IL6R':(0.542,0.010),'STAT4':(0.407,0.108),'JAK1':(0.383,0.093),'IL6ST':(0.353,0.106),
        'STAT3':(0.289,0.026),'JAK2':(0.281,0.227),'STAT1':(0.270,0.279),'STAT5B':(0.112,0.489),
        'TYK2':(0.036,0.892),'SOCS1':(0.033,0.913),'STAT5A':(-0.052,0.802),'STAT6':(-0.217,0.306),
        'JAK3':(-0.074,0.611)}
recon_rows = []
for _, r in primary.iterrows():
    g = r['gene']
    m = manu.get(g, (np.nan, np.nan))
    flag = ""
    if abs(r['rho'] - m[0]) > 0.02 or abs(r['p_spin'] - m[1]) > 0.02:
        flag = "  <== 不一致"
    print(f"  {g:8s} {r['rho']:8.4f} {r['p_spin']:10.4f} {r['q_spin']:10.4f}   rho={m[0]:+.3f} p={m[1]:.3f}{flag}")
    recon_rows.append(dict(gene=g, rho_file=r['rho'], p_spin_file=r['p_spin'],
                           q_spin_file=r['q_spin'], rho_manuscript=m[0], p_spin_manuscript=m[1],
                           consistent=abs(r['rho']-m[0])<=0.02 and abs(r['p_spin']-m[1])<=0.02))
pd.DataFrame(recon_rows).to_csv(OUT / "part0_table1_reconciliation.csv", index=False)

# ============================================================================
# Part A — MAGMA 随机基因集零分布
# ============================================================================
hdr("Part A  MAGMA 随机基因集零分布（三匹配：表达量 / 基因长度 / SNP 数）")

mg = pd.read_csv(BASE / "FUMA_JOB757731/magma.genes.out", sep=r"\s+")
mg = mg.dropna(subset=['SYMBOL', 'ZSTAT'])
mg['gene'] = mg['SYMBOL'].astype(str)
mg['len'] = (mg['STOP'] - mg['START']).astype(float)
print(f"  MAGMA 基因数 = {len(mg)}")

hit = mg[mg['gene'].isin(GENES13)].drop_duplicates('gene')
missing = set(GENES13) - set(hit['gene'])
print(f"  命中 13 基因中的 {len(hit)} 个; 缺失: {sorted(missing) if missing else '无'}")
print(hit[['gene','CHR','ZSTAT','P','NSNPS','len']].to_string(index=False))

Z_OBS   = hit['ZSTAT'].to_numpy(float)
N_OBS   = int((Z_OBS > norm.ppf(0.95)).sum())          # Z > 1.645 == P < 0.05
OE_OBS  = N_OBS / (len(Z_OBS) * 0.05)
print(f"\n  观测: mean Z = {Z_OBS.mean():+.4f}   n(Z>1.645) = {N_OBS}   O/E = {OE_OBS:.3f}")
print(f"  论文正文: 3/13 (P<0.05), 4.6 倍富集")

# 表达量（Schaefer-400 全基因均值，作表达水平代理）
print("  读取全基因 Schaefer-400 表达矩阵（匹配用）...")
expr = pd.read_csv(BASE / "expr_schaefer400_allgenes.csv", index_col=0)
expr.index = expr.index.astype(str)
mean_expr = expr.mean(axis=0)
print(f"    表达矩阵: {expr.shape[0]} parcels x {expr.shape[1]} genes")

pool = mg.copy()
pool['log_len']   = np.log10(pool['len'].clip(lower=1))
pool['log_nsnp']  = np.log10(pool['NSNPS'].clip(lower=1))
pool['log_expr']  = np.log10(np.clip(mean_expr.reindex(pool['gene']).to_numpy(float), 1e-3, None))
pool = pool.dropna(subset=['log_expr'])

def decile_bins(s, n_bins=10):
    return pd.qcut(s.rank(method='first'), n_bins, labels=False, duplicates='drop')

for c in ['log_len', 'log_nsnp', 'log_expr']:
    pool[c + '_bin'] = decile_bins(pool[c])
pool['key'] = (pool['log_len_bin'].astype(str) + '_' +
               pool['log_nsnp_bin'].astype(str) + '_' + pool['log_expr_bin'].astype(str))

key_of = dict(zip(pool['gene'], pool['key']))
keys_obs = [key_of.get(g) for g in hit['gene']]
print(f"  观测基因的匹配格: {keys_obs}")
print(f"  匹配格覆盖基因数: " + ", ".join(
    f"{k}={int((pool['key']==k).sum())}" for k in sorted(set(keys_obs)) if k is not None))

rng = np.random.default_rng(RS)
meanZ_null = np.empty(N_DRAW)
count_null = np.empty(N_DRAW, dtype=int)
meanZ_unmatched = np.empty(N_DRAW)
count_unmatched = np.empty(N_DRAW, dtype=int)

pool_z  = pool['ZSTAT'].to_numpy(float)
by_key  = {k: pool.loc[pool['key'] == k, 'ZSTAT'].to_numpy(float) for k in pool['key'].unique()}

for b in range(N_DRAW):
    zs = np.empty(len(Z_OBS))
    for j, k in enumerate(keys_obs):
        if (k is not None) and (k in by_key) and (len(by_key[k]) > 0):
            zs[j] = by_key[k][rng.integers(len(by_key[k]))]
        else:
            zs[j] = pool_z[rng.integers(len(pool_z))]
    meanZ_null[b]  = zs.mean()
    count_null[b]  = int((zs > norm.ppf(0.95)).sum())
    zu = pool_z[rng.integers(0, len(pool_z), size=len(Z_OBS))]
    meanZ_unmatched[b]  = zu.mean()
    count_unmatched[b]  = int((zu > norm.ppf(0.95)).sum())

def emp_p(obs, null, greater=True):
    return (np.sum(null >= obs) + 1) / (len(null) + 1) if greater else (np.sum(null <= obs) + 1) / (len(null) + 1)

p_meanZ = emp_p(Z_OBS.mean(), meanZ_null, True)
p_count = emp_p(N_OBS, count_null, True)
print(f"\n  --- 三匹配置换（10,000 次）---")
print(f"  mean Z   : 观测 {Z_OBS.mean():+.4f}   零分布 {meanZ_null.mean():+.4f} ± {meanZ_null.std():.4f}   经验 P = {p_meanZ:.4f}")
print(f"  n(P<0.05): 观测 {N_OBS}          零分布 {count_null.mean():.2f} ± {count_null.std():.2f}     经验 P = {p_count:.4f}")
print(f"  零分布中 n>=3 的比例 = {np.mean(count_null >= N_OBS):.4f}")
print(f"  --- 未匹配对照 ---")
print(f"  mean Z   : 零分布 {meanZ_unmatched.mean():+.4f} ± {meanZ_unmatched.std():.4f}   经验 P = {emp_p(Z_OBS.mean(), meanZ_unmatched, True):.4f}")
print(f"  n(P<0.05): 零分布 {count_unmatched.mean():.2f} ± {count_unmatched.std():.2f}     经验 P = {emp_p(N_OBS, count_unmatched, True):.4f}")

pd.DataFrame(dict(stat=['mean_Z','n_P_lt_0.05','O_over_E'],
                  observed=[Z_OBS.mean(), N_OBS, OE_OBS],
                  null_mean_matched=[meanZ_null.mean(), count_null.mean(), count_null.mean()/(len(Z_OBS)*0.05)],
                  null_sd_matched=[meanZ_null.std(), count_null.std(), count_null.std()/(len(Z_OBS)*0.05)],
                  emp_p_matched=[p_meanZ, p_count, p_count],
                  null_mean_unmatched=[meanZ_unmatched.mean(), count_unmatched.mean(), count_unmatched.mean()/(len(Z_OBS)*0.05)],
                  emp_p_unmatched=[emp_p(Z_OBS.mean(), meanZ_unmatched, True), emp_p(N_OBS, count_unmatched, True), np.nan])
          ).to_csv(OUT / "partA_magma_geneset_null.csv", index=False)
np.save(OUT / "partA_count_null_matched.npy", count_null)
np.save(OUT / "partA_meanZ_null_matched.npy", meanZ_null)
print(f"\n  >>> 判读: 经验 P 若 > 0.05，则'4.6 倍富集'不成立，必须从摘要撤下")

# ============================================================================
# Part B/C — 空间 spin（依赖 templateflow fsaverage sphere）
# ============================================================================
hdr("Part B/C  空间 spin 检验（基因集级 + 网格级）")

try:
    import nibabel as nib
    from scipy.spatial import cKDTree
except Exception as e:
    print(f"  [FAIL] 缺少 nibabel/scipy: {e}")
    sys.exit(0)

TF = Path(r"C:/Users/ThinkPad/AppData/Local/templateflow/templateflow/Cache/tpl-fsaverage")
SPH = {h: TF / f"tpl-fsaverage_hemi-{h}_den-10k_sphere.surf.gii" for h in ['L', 'R']}
if not all(p.exists() for p in SPH.values()):
    print(f"  [FAIL] 找不到 fsaverage sphere: {list(SPH.values())}")
    sys.exit(0)

coords_l = nib.load(str(SPH['L'])).darrays[0].data
coords_r = nib.load(str(SPH['R'])).darrays[0].data
print(f"  Sphere 顶点数: L={coords_l.shape[0]}  R={coords_r.shape[0]}")

ANNOT_DIR = BASE / "optimized_figures/schaefer_annot"
def load_annot(path):
    labels, ctab, names = nib.freesurfer.read_annot(str(path))
    mapping = {int(lb): i for i, lb in enumerate(sorted(set(labels)))}
    return np.array([mapping[int(lb)] for l in labels], dtype=np.int32)

lh_lab = load_annot(ANNOT_DIR / "lh.Schaefer2018_400Parcels_7Networks_order.annot")
rh_lab = load_annot(ANNOT_DIR / "rh.Schaefer2018_400Parcels_7Networks_order.annot")
print(f"  Annot 顶点数: L={len(lh_lab)} R={len(rh_lab)}")
if len(lh_lab) != coords_l.shape[0] or len(rh_lab) != coords_r.shape[0]:
    print("  [FAIL] annot 与 sphere 顶点数不匹配")
    sys.exit(0)

labels_all = np.hstack([lh_lab, rh_lab])
n_par = 400
counts = np.bincount(labels_all, minlength=n_par).astype(float)

def surf2parc(v):
    s = np.zeros(n_par); np.add.at(s, labels_all, v)
    return s / counts

def parc2surf(p):
    return p[labels_all]

# ---- Freedman–Perry spin：随机旋转球面 + 最近邻回贴（自实现，无 brainspace 依赖）----
print(f"  生成 spin 索引矩阵 (n_rep={N_SPIN}, 球面顶点 {coords_l.shape[0]+coords_r.shape[0]})...")
t0 = time.time()
coords_all = np.vstack([coords_l, coords_r])           # (n_vert, 3)
tree_l = cKDTree(coords_l)
tree_r = cKDTree(coords_r)
n_l = coords_l.shape[0]
rng_spin = np.random.default_rng(RS)

def rand_rot(rng):
    A = rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q

# 标准 spin：同一旋转矩阵作用于两半球，但各自在本半球内最近邻回贴（L 不跨到 R）
IDX = np.empty((N_SPIN, coords_all.shape[0]), dtype=np.int32)
for b in range(N_SPIN):
    Q = rand_rot(rng_spin)
    _, idx_l = tree_l.query(coords_l @ Q.T, k=1)
    _, idx_r = tree_r.query(coords_r @ Q.T, k=1)
    IDX[b, :n_l] = idx_l
    IDX[b, n_l:] = idx_r + n_l
    if (b + 1) % 500 == 0:
        print(f"    {b+1}/{N_SPIN} ({time.time()-t0:.1f}s)")
print(f"  索引就绪 ({time.time()-t0:.1f}s)")

parc_cols = [np.where(labels_all == p)[0] for p in range(n_par)]

# ---- 数据 ----
expr13 = expr.loc[:, [g for g in GENES13 if g in expr.columns]]
X = expr13.to_numpy(float)                      # 400 x 13
Xr = np.apply_along_axis(rankdata, 0, X)
genes_used = list(expr13.columns)

trait_files = sorted((BASE / "traits_schaefer400").glob("trait_*.csv"))
traits = {}
for f in trait_files:
    df = pd.read_csv(f)
    col = [c for c in df.columns if c != 'parcel'][0]
    traits[f.stem.replace('trait_', '')] = df[col].to_numpy(float)
print(f"  载入 {len(traits)} 个表型 x {X.shape[1]} 个基因 = {len(traits)*X.shape[1]} 个组合")

n_par = 400
Xr_n = (Xr - Xr.mean(0)) / Xr.std(0)

def rank_norm(v):
    r = rankdata(v)
    return (r - r.mean()) / r.std()

def trait_spin_parcels(yvec):
    """对一个表型做 spin，返回 (400, n_rep) 的置换表型矩阵。"""
    y_vert = yvec[labels_all]                      # (n_vert,)
    V = y_vert[IDX]                                # (n_rep, n_vert)
    M = np.empty((n_par, V.shape[0]), float)
    for p in range(n_par):
        M[p] = V[:, parc_cols[p]].mean(axis=1)
    return M

# ---- Part B: 基因集级统计量（主表型）----
hdr("Part B  基因集级 spin：13 基因 mean|rho|")
key_primary = 'mddadult_case-controls_CortThick'
Yp = traits[key_primary]
Yp_rn = rank_norm(Yp)
obs_rho = (Xr_n.T @ Yp_rn) / (n_par - 1)
T_obs = float(np.mean(np.abs(obs_rho)))
print(f"  主表型 = {key_primary}")
print(f"  13 基因 rho: " + ", ".join(f"{g}={r:+.3f}" for g, r in zip(genes_used, obs_rho)))
print(f"  T_obs = mean|rho| = {T_obs:.4f}")

print("  生成主表型 spin 置换 ...")
t0 = time.time()
Yperm_primary = trait_spin_parcels(Yp)
print(f"  完成，shape={Yperm_primary.shape} ({time.time()-t0:.1f}s)")
Yp_r = np.apply_along_axis(rankdata, 0, Yperm_primary)
Yp_rn = (Yp_r - Yp_r.mean(0)) / Yp_r.std(0)
C = Xr_n.T @ Yp_rn / (n_par - 1)                 # 13 x n_rep
T_null = np.mean(np.abs(C), axis=0)
p_set = (np.sum(T_null >= T_obs) + 1) / (len(T_null) + 1)
print(f"  零分布 mean|rho| = {T_null.mean():.4f} ± {T_null.std():.4f}")
print(f"  >>> 基因集级 spin 经验 P = {p_set:.4f}")
pd.DataFrame(dict(statistic=['mean_abs_rho_13genes'], observed=[T_obs],
                  null_mean=[T_null.mean()], null_sd=[T_null.std()], emp_p=[p_set])
              ).to_csv(OUT / "partB_geneset_spin.csv", index=False)
np.save(OUT / "partB_T_null.npy", T_null)
np.save(OUT / "partB_rho13_primary.npy", obs_rho)

# ---- Part C: 网格级 spin（286 组合）----
hdr("Part C  网格级 spin：286 组合的极端性")
tkeys = list(traits.keys())
print(f"  逐表型 spin（{len(tkeys)} 个）...")
t0 = time.time()
Rho_obs = np.empty((X.shape[1], len(tkeys)), float)
C_all = np.empty((X.shape[1], len(tkeys), N_SPIN), float)
for ti, k in enumerate(tkeys):
    y = traits[k]
    Rho_obs[:, ti] = (Xr_n.T @ rank_norm(y)) / (n_par - 1)
    M = trait_spin_parcels(y)
    Mr = np.apply_along_axis(rankdata, 0, M)
    Mrn = (Mr - Mr.mean(0)) / Mr.std(0)
    n_act = min(N_SPIN, Mrn.shape[1])
    C_all[:, ti, :n_act] = (Xr_n.T @ Mrn[:, :n_act]) / (n_par - 1)
    if (ti + 1) % 5 == 0:
        print(f"    {ti+1}/{len(tkeys)} ({time.time()-t0:.1f}s)")
n_act = min(N_SPIN, C_all.shape[2])
C_all = C_all[:, :, :n_act]
print(f"  完成 ({time.time()-t0:.1f}s), n_rep={n_act}")

abs_obs = np.abs(Rho_obs).ravel()
abs_obs_sorted = np.sort(abs_obs)[::-1]
K = int((dk['p_spin'] < 0.05).sum())
print(f"  观测 286 组合 |rho|: mean = {abs_obs.mean():.4f}  max = {abs_obs.max():.4f}")
print(f"  取 k = {K}（DK68 中 p_spin<0.05 的个数）")
T_top_obs = float(abs_obs_sorted[:K].sum())

A = np.abs(C_all).reshape(-1, n_act)            # 286 x n_rep
T_mean_null = A.mean(axis=0)
T_top_null  = np.sort(A, axis=0)[::-1, :][:K, :].sum(axis=0)

p_top  = (np.sum(T_top_null >= T_top_obs) + 1) / (len(T_top_null) + 1)
p_mean = (np.sum(T_mean_null >= abs_obs.mean()) + 1) / (len(T_mean_null) + 1)
print(f"  top-{K} sum|rho|: 观测 {T_top_obs:.4f}   零分布 {T_top_null.mean():.4f} ± {T_top_null.std():.4f}   经验 P = {p_top:.4f}")
print(f"  mean  |rho|    : 观测 {abs_obs.mean():.4f}   零分布 {T_mean_null.mean():.4f} ± {T_mean_null.std():.4f}   经验 P = {p_mean:.4f}")

pd.DataFrame(dict(statistic=[f'top_{K}_sum_abs_rho', 'mean_abs_rho_grid'],
                  observed=[T_top_obs, abs_obs.mean()],
                  null_mean=[T_top_null.mean(), T_mean_null.mean()],
                  null_sd=[T_top_null.std(), T_mean_null.std()],
                  emp_p=[p_top, p_mean])
              ).to_csv(OUT / "partC_grid_spin.csv", index=False)
np.save(OUT / "partC_T_top_null.npy", T_top_null)
np.save(OUT / "partC_T_mean_null.npy", T_mean_null)

hdr("汇总")
print(f"  Part 0  论文 15 / DK68 {int((dk['p_spin']<0.05).sum())} / Schaefer {int((sc['p_spin']<0.05).sum())}  —— 三者不一致")
print(f"  Part A  MAGMA 13 基因 n(P<0.05)={N_OBS} vs 匹配随机集 {count_null.mean():.2f}±{count_null.std():.2f}, 经验 P={p_count:.4f}")
print(f"  Part B  基因集 mean|rho|={T_obs:.4f}, 经验 P={p_set:.4f}")
print(f"  Part C  网格 top-{K} sum|rho| 经验 P={p_top:.4f}; mean|rho| 经验 P={p_mean:.4f}")
print(f"\n  全部结果 -> {OUT}")
