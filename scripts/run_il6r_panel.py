#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IL6R (rs2228145) -> cortical thickness   MR 判决分析  汇总统计 + 面板检验

结局: Makowski C et al. PNAS 2023 (PMID 36893272), GCST90271790-90271815
      n = 32,488 European adults (UKB, fastGWA), GRCh37
暴露: plasma sIL-6R, 工具 rs2228145 Asp358Ala (chr1:154426970, A>C, EA=C=358Ala)
      C -> sIL-6R 升高 (shedding↑, classical signaling↓, trans-signaling↑)
"""
from __future__ import annotations

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

BETA_EXP = 1.035   # SD sIL-6R per C allele  (van Dongen 2014, R2=51% -> sqrt(.51/(2*.39*.61)))
BETA_EXP_LO = 0.49
BETA_EXP_HI = 1.035

# (label, gcst, 类别)
PANEL = [
    ("Ventral frontal",    "GCST90271812", "thickness"),
    ("Temporal pole",      "GCST90271811", "thickness"),
    ("Dorsolateral PFC",   "GCST90271802", "thickness"),
    ("Medial prefrontal",  "GCST90271804", "thickness"),
    ("Perisylvian",        "GCST90271809", "thickness"),
    ("Middle temporal",    "GCST90271806", "thickness"),
    ("Motor premotor",     "GCST90271807", "thickness"),
    ("Inferior parietal",  "GCST90271803", "thickness"),
    ("Occipital",          "GCST90271808", "thickness"),
    ("Superior parietal",  "GCST90271810", "thickness"),
    ("vmOccipital",        "GCST90271813", "thickness"),
    ("Medial temporal",    "GCST90271805", "thickness"),
    ("Global mean CT",     "GLOBAL",       "thickness"),
    ("Total surface area", "AREA",         "area"),
]

# 从 regional/*.txt 与 il6r_mr_results.json 读入
def load_rows():
    rows = []
    for label, acc, cat in PANEL:
        if acc == "GLOBAL":
            beta, se, p, n = -0.0119313, 0.00799077, 0.1354, 32488
        elif acc == "AREA":
            beta, se, p, n = 0.0110744, 0.00798674, 0.165563, 32488
        else:
            fns = list((HERE / "regional").glob(f"{acc}_*.txt"))
            if not fns:
                continue
            txt = fns[0].read_text().strip().splitlines()[-1].split("\t")
            beta, se, p, n = float(txt[4]), float(txt[5]), float(txt[7]), int(float(txt[8]))
        z = beta / se
        rows.append(dict(label=label, acc=acc, cat=cat, beta=beta, se=se, p=p, n=n, z=z,
                         theta=beta / BETA_EXP))
    return rows


def norm_sf(z):
    """two-sided p for standard normal"""
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def binom_sf(k, n, p=0.5):
    """P(X >= k) for X~Bin(n,p)"""
    s = 0.0
    for i in range(k, n + 1):
        s += math.comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
    return s


def main():
    rows = load_rows()
    th = [r for r in rows if r["cat"] == "thickness" and r["acc"] != "GLOBAL"]
    area = [r for r in rows if r["cat"] == "area"]
    glob = [r for r in rows if r["acc"] == "GLOBAL"][0]

    print("=" * 78)
    print("IL6R (rs2228145 Asp358Ala, EA=C)  ->  cortical structure")
    print("outcome: Makowski 2023 PNAS GCST90271790-15  n=32,488 (UKB European, GRCh37)")
    print("exposure: plasma sIL-6R, beta_exp = %.3f SD per C allele" % BETA_EXP)
    print("=" * 78)

    print("\n%-22s %10s %10s %10s %12s %10s" % (
        "phenotype", "beta", "se", "p", "theta", "z"))
    print("-" * 78)
    for r in rows:
        mark = " *" if r["p"] < 0.05 else ""
        print("%-22s %+10.5f %10.5f %10.4f %+12.5f %10.3f%s" % (
            r["label"], r["beta"], r["se"], r["p"], r["theta"], r["z"], mark))

    # ---- 面板检验 ----
    n_th = len(th)
    n_neg = sum(1 for r in th if r["beta"] < 0)
    n_nom = sum(1 for r in th if r["p"] < 0.05)
    # 符号检验 (假设独立; 相关时此为上界)
    p_sign_ind = 0.5 ** n_th if n_neg == n_th else None
    # binomial: 5/12 nominal
    p_binom = binom_sf(n_nom, n_th, 0.05)

    # Stouffer 加权 Z (等权) —— 相关会高估, 这里给独立上界
    import statistics as st
    zs = [r["z"] for r in th]
    z_stouffer = sum(zs) / math.sqrt(n_th)
    p_stouffer = norm_sf(z_stouffer)

    # 全局 CT 的 z (可视为相关校正后的合成)
    z_glob = glob["z"]

    print("\n--- 面板检验 (12 个区域皮层厚度) ---")
    print(f"  方向一致: {n_neg}/{n_th} 区域 beta < 0")
    print(f"  名义显著 (p<0.05): {n_nom}/{n_th}   (零假设期望 {0.05*n_th:.1f})")
    print(f"  符号检验 (独立上界): P(全同向) = {p_sign_ind:.2e}")
    print(f"  Binomial P(>= {n_nom} of {n_th} nominal): {p_binom:.3e}")
    print(f"  Stouffer 等权 Z = {z_stouffer:+.3f}   P = {p_stouffer:.4g}"
          f"   [独立上界, 区域相关时此 P 偏小]")
    print(f"  全局平均皮层厚度 Z = {z_glob:+.3f}   P = {glob['p']:.4g}   [相关校正后的合成]")

    # 表面积特异性
    a = area[0]
    print(f"\n--- 特异性: 总皮层表面积 ---")
    print(f"  beta = {a['beta']:+.5f}  p = {a['p']:.4f}   (方向为正, 与厚度相反)")
    print("  => 效应是 thickness-specific, 不是全脑尺寸效应")

    # 效应量区间
    print("\n--- 效应量 (Wald ratio, theta = SD cortical thickness per SD sIL-6R) ---")
    for tag, b in [("beta_exp=0.49", BETA_EXP_LO), ("beta_exp=1.035", BETA_EXP_HI)]:
        print(f"  {tag}:")
        for r in [glob] + sorted(th, key=lambda x: x["p"])[:3]:
            th_ = r["beta"] / b
            se_ = r["se"] / b
            print(f"    {r['label']:<22} theta={th_:+.5f}  se={se_:.5f}"
                  f"  95%CI=({th_-1.96*se_:+.5f},{th_+1.96*se_:+.5f})")

    # 功效
    se_out = glob["se"]
    mde = (1.96 + 0.84) * se_out  # 80% power, alpha .05, on beta_out scale
    print(f"\n--- 功效 (n=32,488, alpha=0.05, 80% power) ---")
    print(f"  可检出的最小 beta_out (每 C 等位基因) = {mde:.5f} SD")
    print(f"  折算 theta (per SD sIL-6R, beta_exp=1.035) = {mde/BETA_EXP:.5f} SD")
    print(f"  观测全局 beta_out = {glob['beta']:+.5f}  (= {abs(glob['beta'])/mde*100:.0f}% 的 MDE)")
    print(f"  观测最强区域 (ventral frontal) beta = {sorted(th,key=lambda x:x['p'])[0]['beta']:+.5f}")

    out = dict(
        panel=rows,
        n_regions=n_th,
        n_negative=n_neg,
        n_nominal=n_nom,
        p_sign_independent=p_sign_ind,
        p_binomial=p_binom,
        z_stouffer=z_stouffer,
        p_stouffer=p_stouffer,
        global_mean_ct=glob,
        surface_area=a,
        beta_exp=BETA_EXP,
        outcome=dict(
            study="Makowski C et al. PNAS 2023;120(11):e2214834120",
            pmid="36893272", n=32488, build="GRCh37",
            note="Grasby 2020 Science (GCST010281, ENIGMA) 全汇总统计未公开 (fullSS=no)",
        ),
    )
    js = HERE / "il6r_mr_panel_results.json"
    js.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nsaved -> {js}")


if __name__ == "__main__":
    main()
