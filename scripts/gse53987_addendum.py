#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GSE53987 补充分析：SCZ 明细 / 功效(MDE) / 基因集方向一致 / 判决汇总"""
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

HERE = Path(__file__).resolve().parent
res = pd.read_csv(HERE / "gse53987_de_results.tsv", sep="\t")
best = res[res["is_best_probe"]].copy()

GENES13 = ["IL6R", "IL6ST", "JAK1", "JAK2", "JAK3", "TYK2",
           "STAT1", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "SOCS1"]
REGIONS = ["PFC", "HIP", "STR"]

print("=" * 78)
print("A. SCZ vs HC 明细（基因级代表探针）")
print("=" * 78)
for reg in REGIONS:
    sub = best[(best.contrast == "SCZ_vs_HC") & (best.region == reg)].sort_values("p_ttest")
    print(f"\n--- {reg} ---")
    print(f"{'gene':<8} {'logFC':>9} {'p_ttest':>10} {'p_mwu':>10} {'p_adjOLS':>10} {'FDR13':>9}")
    for r in sub.itertuples():
        star = " *" if r.p_ttest < 0.05 else ""
        print(f"{r.gene:<8} {r.logFC:>+9.4f} {r.p_ttest:>10.4g} {r.p_mwu:>10.4g} "
              f"{r.p_ttest_adj:>10.4g} {r.fdr_within_13:>9.3g}{star}")

print("\n" + "=" * 78)
print("B. 功效 / 最小可检出效应 (MDE), alpha=0.05, 80% power")
print("=" * 78)
print(f"{'region':<6} {'contrast':<11} {'n_case':>6} {'n_ctrl':>7} {'pooled_SD':>10} {'MDE_log2':>9} {'MDE_fold':>9}")
for cn in ["MDD_vs_HC", "SCZ_vs_HC", "BD_vs_HC"]:
    for reg in REGIONS:
        sub = best[(best.contrast == cn) & (best.region == reg)]
        if len(sub) == 0:
            continue
        n1, n2 = int(sub.n_case.iloc[0]), int(sub.n_ctrl.iloc[0])
        # pooled SD approximated from t and logFC:  se = sd*sqrt(1/n1+1/n2) -> sd = se / sqrt(1/n1+1/n2)
        tvals = pd.to_numeric(sub.t, errors="coerce").replace(0, np.nan)
        se = np.abs(sub.logFC.values) / np.abs(tvals.values)
        se = se[np.isfinite(se)]
        sd = float(np.median(se) / np.sqrt(1 / n1 + 1 / n2)) if len(se) else float("nan")
        mde = (1.96 + 0.84) * sd * np.sqrt(1 / n1 + 1 / n2)
        print(f"{reg:<6} {cn:<11} {n1:>6d} {n2:>7d} {sd:>10.3f} {mde:>9.3f} {2**mde:>9.3f}")

print("\n" + "=" * 78)
print("C. 基因集方向一致（13 基因，MDD vs HC）")
print("=" * 78)
for reg in REGIONS:
    sub = best[(best.contrast == "MDD_vs_HC") & (best.region == reg)]
    if len(sub) == 0:
        continue
    up = int((sub.logFC > 0).sum())
    dn = int((sub.logFC < 0).sum())
    mean_lfc = sub.logFC.mean()
    # 符号检验
    p_sign = stats.binomtest(up, up + dn, 0.5, alternative="two-sided").pvalue
    # 与 SCZ 对照
    sub_s = best[(best.contrast == "SCZ_vs_HC") & (best.region == reg)]
    up_s = int((sub_s.logFC > 0).sum())
    print(f"{reg}: MDD up={up}/13 down={dn}/13  mean logFC={mean_lfc:+.4f}  sign-test p={p_sign:.3f}"
          f"   |  SCZ up={up_s}/13")

print("\n" + "=" * 78)
print("D. 判决汇总")
print("=" * 78)
rows = []
for cn, lab in [("MDD_vs_HC", "MDD"), ("SCZ_vs_HC", "SCZ"), ("BD_vs_HC", "BD")]:
    for reg in REGIONS:
        sub = best[(best.contrast == cn) & (best.region == reg)]
        if len(sub) == 0:
            continue
        rows.append(dict(
            group=lab, region=reg,
            n_nominal=int((sub.p_ttest < 0.05).sum()),
            n_nominal_adj=int((sub.p_ttest_adj < 0.05).sum()),
            n_fdr13=int((sub.fdr_within_13 < 0.05).sum()),
            min_p=float(sub.p_ttest.min()),
            min_gene=sub.loc[sub.p_ttest.idxmin(), "gene"],
            min_logFC=float(sub.loc[sub.p_ttest.idxmin(), "logFC"]),
            max_abs_logFC=float(sub.logFC.abs().max()),
        ))
df = pd.DataFrame(rows)
print(df.to_string(index=False))

print("\n--- 一句话判决 ---")
mdd = df[df.group == "MDD"]
print(f"MDD: 3 区域 × 13 基因 = {len(best[best.contrast=='MDD_vs_HC'])} 个基因级检验")
print(f"     名义 p<0.05: {int(mdd.n_nominal.sum())} 个；协变量校正后 p<0.05: {int(mdd.n_nominal_adj.sum())} 个")
print(f"     FDR<0.05: {int(mdd.n_fdr13.sum())} 个；最小 p = {mdd.min_p.min():.4f} ({mdd.loc[mdd.min_p.idxmin(),'min_gene']} / {mdd.loc[mdd.min_p.idxmin(),'region']})")
print(f"     最大 |logFC| = {mdd.max_abs_logFC.max():.3f} log2 (= {2**mdd.max_abs_logFC.max():.2f}-fold)")
scz = df[df.group == "SCZ"]
print(f"SCZ: 名义 {int(scz.n_nominal.sum())} 个；FDR<0.05: {int(scz.n_fdr13.sum())} 个；最小 p = {scz.min_p.min():.2e}")
