#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
检验 1 — GSE53987 人死后脑 JAK-STAT 13 基因差异表达
====================================================
数据: GSE53987 (PMID 31123247)  Affymetrix HG-U133 Plus 2 (GPL570), RMA 归一化
      三脑区: PFC(BA46) / 海马 / 纹状体
      四组: schizophrenia / bipolar / MDD / control
主要对比:  MDD vs control（判决性）
参照对比:  SCZ vs control、BD vs control（疾病特异性）

13 基因（先验选定，来自 JAK-STAT 通路）:
  IL6R IL6ST JAK1 JAK2 JAK3 TYK2 STAT1 STAT3 STAT4 STAT5A STAT5B STAT6 SOCS1

方法:
  - **GSE53987 矩阵为线性强度（4.59–25,856），已 log2 变换后再检验**
    （RMA 后未取对数；微阵列标准做法是 log2，否则方差随均值膨胀、t 检验失真）
  - logFC = mean(case) - mean(control)  (log2 单位，即 log2 倍数变化)
  - Welch t 检验 + Mann-Whitney U
  - 协变量校正: OLS  ~ group + age + gender + pmi + ph + rin
  - 多探针处理（预设规则，避免挑 p 最小的）: 取该基因在全样本中**平均表达最高**的探针
    作为基因级代表，同时报告全部探针
  - 多重校正: BH-FDR，跨 13 基因（每个 region×contrast）+ 跨 13×3=39 个主检验

输出: gse53987_de_results.tsv / .json + 控制变量平衡表
"""
from __future__ import annotations

import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
MAT = HERE / "GSE53987_series_matrix.txt.gz"
PLAT = HERE / "GPL570_platform_soft.txt"

GENES13 = ["IL6R", "IL6ST", "JAK1", "JAK2", "JAK3", "TYK2",
           "STAT1", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "SOCS1"]

REGION_MAP = {
    "prefrontal cortex": "PFC", "prefrontal": "PFC", "pfc": "PFC", "ba46": "PFC",
    "hippocampus": "HIP", "hip": "HIP",
    "striatum": "STR", "str": "STR",
}
DISEASE_MAP = {
    "control": "HC", "healthy control": "HC", "normal": "HC", "unaffected": "HC",
    "mdd": "MDD", "major depressive disorder": "MDD", "depression": "MDD",
    "bipolar disorder": "BD", "bipolar": "BD", "bpd": "BD",
    "schizophrenia": "SCZ", "scz": "SCZ",
}


# ---------------------------------------------------------------- 平台注释
def load_platform() -> pd.DataFrame:
    """probe -> gene symbols (/// 分隔)"""
    rows = []
    with open(PLAT, "r", encoding="utf-8", errors="ignore") as f:
        in_tab = False
        hdr = None
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_tab = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_tab:
                continue
            c = line.split("\t")
            if hdr is None:
                hdr = c
                idx = {k: i for i, k in enumerate(hdr)}
                continue
            sym = c[idx["Gene Symbol"]] if len(c) > idx["Gene Symbol"] else ""
            rows.append((c[idx["ID"]], sym))
    df = pd.DataFrame(rows, columns=["probe", "gene_symbol"])
    return df


# ---------------------------------------------------------------- 矩阵 + 元数据
def load_matrix():
    meta = {}
    accs, titles = None, None
    chars = defaultdict(list)
    char_i = 0
    expr_rows = []
    hdr = None
    with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
        in_tab = False
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_tab = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_tab:
                c = line.split("\t")
                if hdr is None:
                    hdr = [x.strip('"') for x in c]
                    continue
                expr_rows.append(c)
                continue
            if line.startswith("!Sample_"):
                parts = line.split("\t")
                tag = parts[0]
                vals = [p.strip('"') for p in parts[1:]]
                if tag == "!Sample_geo_accession":
                    accs = vals
                elif tag == "!Sample_title":
                    titles = vals
                elif tag == "!Sample_characteristics_ch1":
                    chars[char_i] = vals
                    char_i += 1
                else:
                    meta[tag] = vals
    samples = hdr[1:]
    expr = pd.DataFrame(expr_rows)
    probes = expr.iloc[:, 0].str.strip('"')
    vals = expr.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    vals.index = probes.values
    vals.columns = samples
    # 重要：GSE53987 矩阵为**线性强度**（实测 4.59–25,856），非 log2。
    # 微阵列标准做法：先 log2 变换再做检验。此后所有 logFC 单位 = log2。
    n_nonpos = int((vals <= 0).sum().sum())
    if n_nonpos:
        print(f"    [warn] {n_nonpos} 个非正值，已裁到 1e-6 后取 log2")
        vals = vals.clip(lower=1e-6)
    vals = np.log2(vals)
    print(f"    [log2] 变换后范围 {vals.min().min():.2f} ~ {vals.max().max():.2f} (log2)")
    return vals, samples, accs, titles, chars, meta


def parse_chars(chars):
    """characteristics 行 -> 列名/值 DataFrame"""
    # 每行形如 'age: 52'
    out = {}
    for i, row in chars.items():
        if not row:
            continue
        key = row[0].split(":")[0].strip().lower()
        vals = []
        for v in row:
            if ":" in v:
                vals.append(v.split(":", 1)[1].strip())
            else:
                vals.append(v)
        out[key] = vals
    return pd.DataFrame(out)


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out


def main():
    print("=" * 78)
    print("检验 1 — GSE53987 死后脑 JAK-STAT 13 基因差异表达")
    print("=" * 78)

    plat = load_platform()
    print(f"[1] GPL570 注释: {len(plat)} 探针")
    # 建立 gene -> probes
    gene2probe = defaultdict(list)
    for pr, sym in plat.itertuples(index=False):
        if not sym or sym == "---":
            continue
        for g in sym.split("///"):
            g = g.strip()
            if g:
                gene2probe[g].append(pr)
    missing = [g for g in GENES13 if g not in gene2probe]
    print(f"    13 基因中缺失: {missing if missing else '无'}")
    for g in GENES13:
        print(f"    {g:8s} probes={len(gene2probe.get(g, []))}  {gene2probe.get(g, [])[:6]}")

    expr, samples, accs, titles, chars, meta = load_matrix()
    print(f"\n[2] 表达矩阵: {expr.shape[0]} probes x {expr.shape[1]} samples")

    ch = parse_chars(chars)
    print(f"    characteristics 列: {list(ch.columns)}")
    pheno = pd.DataFrame({
        "gsm": accs,
        "title": titles,
        "disease_raw": ch.get("disease state", pd.Series([""] * len(accs))).values,
        "tissue_raw": ch.get("tissue", pd.Series([""] * len(accs))).values,
    })
    for k in ["age", "gender", "race", "pmi", "ph", "rin"]:
        if k in ch.columns:
            pheno[k] = ch[k].values
    pheno["region"] = pheno["tissue_raw"].str.strip().str.lower().map(REGION_MAP)
    pheno["group"] = pheno["disease_raw"].str.strip().str.lower().map(DISEASE_MAP)
    # 兜底: 从 title 解析 (形如 'bipolar_hip_10')
    if pheno["group"].isna().any() or pheno["region"].isna().any():
        for i, t in enumerate(pheno["title"]):
            tl = t.lower()
            if pd.isna(pheno.loc[i, "group"]):
                for k, v in DISEASE_MAP.items():
                    if k.split()[0] in tl:
                        pheno.loc[i, "group"] = v
                        break
            if pd.isna(pheno.loc[i, "region"]):
                if "hip" in tl:
                    pheno.loc[i, "region"] = "HIP"
                elif "str" in tl:
                    pheno.loc[i, "region"] = "STR"
                elif "pfc" in tl or "pre" in tl:
                    pheno.loc[i, "region"] = "PFC"

    print(f"\n[3] 表型分组")
    print(pd.crosstab(pheno["group"], pheno["region"], margins=True).to_string())

    # 协变量平衡
    print("\n[4] 协变量平衡 (按 group)")
    for cov in ["age", "pmi", "ph", "rin"]:
        if cov not in pheno.columns:
            continue
        pheno[cov] = pd.to_numeric(pheno[cov], errors="coerce")
    bal = pheno.groupby("group")[["age", "pmi", "ph", "rin"]].agg(["count", "mean", "std"])
    print(bal.round(2).to_string())
    pheno["male"] = (pheno["gender"].str.upper() == "M").astype(int)
    print("\n  gender (male%):")
    print(pheno.groupby("group")["male"].agg(["count", "mean"]).round(3).to_string())

    # 数值矩阵对齐
    common = [s for s in pheno["gsm"] if s in expr.columns]
    pheno = pheno.set_index("gsm").loc[common].reset_index()
    expr = expr[common]
    print(f"\n[5] 对齐后样本 {len(common)}")

    # 探针 -> 基因，取平均表达最高的探针为基因级代表
    gene_probe_expr = {}
    for g in GENES13:
        prs = [p for p in gene2probe.get(g, []) if p in expr.index]
        if not prs:
            continue
        sub = expr.loc[prs]
        means = sub.mean(axis=1)
        best = means.idxmax()
        gene_probe_expr[g] = dict(probes=prs, best=best, best_mean=float(means[best]))

    print("\n[6] 基因级代表探针（预设规则：全样本平均表达最高）")
    for g, d in gene_probe_expr.items():
        print(f"    {g:8s} best={d['best']:12s} mean={d['best_mean']:.2f}  n_probes={len(d['probes'])}")

    # ---------------- DE 主分析 ----------------
    CONTRASTS = [
        ("MDD", "HC", "MDD_vs_HC"),
        ("SCZ", "HC", "SCZ_vs_HC"),
        ("BD", "HC", "BD_vs_HC"),
    ]
    REGIONS = ["PFC", "HIP", "STR"]

    records = []
    for reg in REGIONS:
        for case, ctrl, cname in CONTRASTS:
            m_case = (pheno["region"] == reg) & (pheno["group"] == case)
            m_ctrl = (pheno["region"] == reg) & (pheno["group"] == ctrl)
            n_case, n_ctrl = int(m_case.sum()), int(m_ctrl.sum())
            if n_case < 3 or n_ctrl < 3:
                print(f"    [skip] {reg} {cname}: n={n_case}/{n_ctrl}")
                continue
            idx_case = pheno.index[m_case]
            idx_ctrl = pheno.index[m_ctrl]
            gsm_case = pheno.loc[idx_case, "gsm"].tolist()
            gsm_ctrl = pheno.loc[idx_ctrl, "gsm"].tolist()

            for g, d in gene_probe_expr.items():
                for pr in d["probes"]:
                    y = expr.loc[pr]
                    a = y[gsm_case].astype(float).values
                    b = y[gsm_ctrl].astype(float).values
                    lfc = float(np.mean(a) - np.mean(b))
                    # Welch t
                    t, pt = stats.ttest_ind(a, b, equal_var=False)
                    # Mann-Whitney
                    try:
                        u, pu = stats.mannwhitneyu(a, b, alternative="two-sided")
                    except Exception:
                        u, pu = float("nan"), float("nan")
                    # 协变量校正 OLS（组 + 年龄 + 性别 + pmi + rin + ph）
                    padj_p = float("nan")
                    try:
                        sub = pheno.loc[idx_case.tolist() + idx_ctrl.tolist()].copy()
                        yv = np.concatenate([a, b])
                        grp = np.concatenate([np.ones(len(a)), np.zeros(len(b))])
                        covs = ["age", "pmi", "rin", "ph"]
                        X = [grp]
                        for cv in covs:
                            if cv in sub.columns and sub[cv].notna().sum() > 10:
                                X.append(pd.to_numeric(sub[cv], errors="coerce").fillna(
                                    pd.to_numeric(sub[cv], errors="coerce").mean()).values)
                        if "male" in sub.columns:
                            X.append(sub["male"].astype(float).values)
                        X = np.column_stack([np.ones(len(yv))] + X)
                        beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
                        resid = yv - X @ beta
                        dof = len(yv) - X.shape[1]
                        sigma2 = float(resid @ resid) / max(dof, 1)
                        xtx_inv = np.linalg.pinv(X.T @ X)
                        se = math.sqrt(max(sigma2 * xtx_inv[1, 1], 1e-12))
                        t_adj = beta[1] / se
                        padj_p = float(2 * (1 - 0.5 * (1 + math.erf(abs(t_adj) / math.sqrt(2)))))
                        beta_adj = float(beta[1])
                    except Exception:
                        beta_adj = float("nan")

                    records.append(dict(
                        region=reg, contrast=cname, case=case, ctrl=ctrl,
                        n_case=n_case, n_ctrl=n_ctrl,
                        gene=g, probe=pr, is_best_probe=(pr == d["best"]),
                        logFC=lfc, t=float(t) if t == t else float("nan"),
                        p_ttest=float(pt), p_mwu=float(pu),
                        logFC_adj=beta_adj, p_ttest_adj=padj_p,
                    ))

    res = pd.DataFrame(records)
    print(f"\n[7] 共 {len(res)} 条 gene×probe×region×contrast 检验")

    # BH-FDR
    # (a) 基因级代表探针，MDD_vs_HC，跨 13 基因 × 3 区域 = 39
    m = res["is_best_probe"] & (res["contrast"] == "MDD_vs_HC")
    res.loc[m, "fdr_39"] = bh_fdr(res.loc[m, "p_ttest"].values)
    # (b) 基因级代表探针，每个 region×contrast 内跨 13 基因
    res["fdr_within_13"] = np.nan
    for (reg, cn), sub in res[res["is_best_probe"]].groupby(["region", "contrast"]):
        res.loc[sub.index, "fdr_within_13"] = bh_fdr(sub["p_ttest"].values)
    # (c) 全探针全对比
    res["fdr_all"] = bh_fdr(res["p_ttest"].values)

    # ---------------- 输出 ----------------
    out_tsv = HERE / "gse53987_de_results.tsv"
    res.to_csv(out_tsv, sep="\t", index=False)
    print(f"[8] saved -> {out_tsv}")

    # 判决表：基因级代表探针 × MDD_vs_HC
    print("\n" + "=" * 78)
    print("判决表：MDD vs HC（基因级代表探针）")
    print("=" * 78)
    prim = res[(res["is_best_probe"]) & (res["contrast"] == "MDD_vs_HC")].copy()
    prim = prim.sort_values(["region", "p_ttest"])
    hdr = f"{'region':<5} {'gene':<8} {'probe':<12} {'n_case':>6} {'n_ctrl':>7} {'logFC':>9} {'p_ttest':>10} {'p_mwu':>10} {'p_adjOLS':>10} {'FDR13':>9} {'FDR39':>9}"
    print(hdr)
    print("-" * len(hdr))
    for r in prim.itertuples():
        star = " *" if r.p_ttest < 0.05 else ("†" if r.p_ttest < 0.10 else "")
        print(f"{r.region:<5} {r.gene:<8} {r.probe:<12} {r.n_case:>6d} {r.n_ctrl:>7d} "
              f"{r.logFC:>+9.4f} {r.p_ttest:>10.4g} {r.p_mwu:>10.4g} {r.p_ttest_adj:>10.4g} "
              f"{r.fdr_within_13:>9.3g} {r.fdr_39:>9.3g}{star}")

    # 参照对比摘要
    print("\n" + "=" * 78)
    print("疾病特异性摘要（基因级代表探针，p_ttest < 0.05 的个数 / 13）")
    print("=" * 78)
    for cn in ["MDD_vs_HC", "SCZ_vs_HC", "BD_vs_HC"]:
        line = []
        for reg in REGIONS:
            sub = res[(res["is_best_probe"]) & (res["contrast"] == cn) & (res["region"] == reg)]
            if len(sub) == 0:
                continue
            n05 = int((sub["p_ttest"] < 0.05).sum())
            n05n = int((sub["p_ttest_adj"] < 0.05).sum()) if sub["p_ttest_adj"].notna().any() else 0
            nfdr = int((sub["fdr_within_13"] < 0.05).sum())
            line.append(f"{reg}: {n05}/13 nominal, {n05n}/13 协变量校正, {nfdr}/13 FDR")
        print(f"  {cn:12s}  " + " | ".join(line))

    # 面板得分（13 基因平均 z，作为整体通路读数）
    print("\n" + "=" * 78)
    print("JAK-STAT 13 基因面板得分（各基因 z 分数均值）")
    print("=" * 78)
    for reg in REGIONS:
        rows = []
        for g, d in gene_probe_expr.items():
            y = expr.loc[d["best"]]
            for grp in ["HC", "MDD", "SCZ", "BD"]:
                msk = (pheno["region"] == reg) & (pheno["group"] == grp)
                if msk.sum() < 3:
                    continue
                rows.append((g, grp, float(y[pheno.loc[msk, "gsm"]].mean())))
        if not rows:
            continue
        dfp = pd.DataFrame(rows, columns=["gene", "group", "expr"])
        # z by gene
        z = []
        for g, sub in dfp.groupby("gene"):
            mu, sd = sub["expr"].mean(), sub["expr"].std(ddof=1) or 1.0
            for r in sub.itertuples():
                z.append((r.gene, r.group, (r.expr - mu) / sd))
        dz = pd.DataFrame(z, columns=["gene", "group", "z"])
        piv = dz.groupby("group")["z"].mean().round(4)
        # t-test 面板 MDD vs HC
        a = dz[(dz.group == "MDD")]["z"].values
        b = dz[(dz.group == "HC")]["z"].values
        if len(a) > 2 and len(b) > 2:
            t, p = stats.ttest_ind(a, b, equal_var=False)
        else:
            t, p = float("nan"), float("nan")
        print(f"  {reg}:  " + "  ".join(f"{k}={v:+.3f}" for k, v in piv.items())
              + f"   [MDD vs HC panel z: t={t:+.2f}, p={p:.4g}]")

    # JSON
    out_json = HERE / "gse53987_de_results.json"
    summary = {
        "dataset": "GSE53987 (PMID 31123247), Affymetrix HG-U133 Plus 2 (GPL570), RMA",
        "design": "PFC(BA46)/HIP/STR x SCZ/BD/MDD/HC, postmortem human brain",
        "genes": GENES13,
        "n_probes_per_gene": {g: len(gene_probe_expr.get(g, {}).get("probes", [])) for g in GENES13},
        "primary_contrast": "MDD vs HC (per region)",
        "n_tests_primary": int((res["is_best_probe"] & (res["contrast"] == "MDD_vs_HC")).sum()),
        "group_region_counts": pd.crosstab(pheno["group"], pheno["region"]).to_dict(),
        "results_file": str(out_tsv),
    }
    out_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n[9] saved -> {out_json}")
    print("\nDONE")


if __name__ == "__main__":
    main()
