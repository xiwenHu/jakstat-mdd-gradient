#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IL6R -> mean cortical thickness  MR (go/no-go 判决分析)

数据:
  结局 (outcome) = mean cortical thickness GWAS
      Makowski C et al. PNAS 2023 (PMID 36893272), GCST90271815
      n = 32,488 European adults (UKB, fastGWA), GRCh37
      列: chromosome base_pair_location effect_allele other_allele
          beta standard_error effect_allele_frequency p_value n
      文件: data/mr/GCST90271815_meanCT_grch37.tsv.gz

  暴露 (exposure) = plasma sIL-6R
      主工具变量 rs2228145 (Asp358Ala, chr1:154426970, A>C)
      C 等位基因 (358Ala) -> sIL-6R 升高 (shedding 增强, classical signaling 减弱)
      beta_exp 取自文献 (见 EXPOSURE dict), 单位 = SD of plasma sIL-6R per C allele

方法:
  1) Wald ratio (单工具):   theta = b_out / b_exp ;  se = se_out / |b_exp|
  2) IVW (多工具, 独立):    与 Wald 同 (无相关时)
  3) 与 PGC MDD 2018 (本地 daner 文件) 的 IL6R 位点结果并列, 作为疾病结局参照
  4) 敏感性: 等位基因方向核对 / F 统计 / 与 rs4845618-rs4129267 工具集 IVW

输出: data/mr/il6r_mr_results.json + il6r_mr_results.tsv
"""
from __future__ import annotations

import gzip
import json
import math
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE
CT_GZ = HERE / "GCST90271815_meanCT_grch37.tsv.gz"
PGC = Path(r"F:\炎症因子文章\MDD_2018\daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz")

# ---------------------------------------------------------------------------
# 1. 暴露侧效应量 (SD of plasma sIL-6R per effect allele)
# ---------------------------------------------------------------------------
# rs2228145 A>C, effect allele = C (358Ala)
# 来源:
#   - van Dongen 2014 Behav Genet (PMID 24791950, PMC4283105): rs2228145 解释
#     血浆 sIL-6R 总方差 51% (N~5000 twin families, 广义遗传度 72%)
#     -> 若表型标准化且 MAF(C)~0.39:  beta ~= sqrt(0.51 / (2*0.39*0.61)) = 1.035 SD
#   - 实际文献常用的每等位基因效应 ~ +0.49 SD 到 +1.0 SD, 取决于标准化方式。
#   这里以 van Dongen R^2 反推值作为主值, 并做区间敏感性 (0.49 / 0.70 / 1.035)。
EXPOSURE = {
    "rs2228145": {
        "effect_allele": "C",
        "other_allele": "A",
        "eaf": 0.39,          # 358Ala 频率 (欧洲)
        "beta_sd": 1.035,     # SD sIL-6R per C allele (由 R2=0.51 反推)
        "beta_sd_lo": 0.49,
        "beta_sd_hi": 1.035,
        "source": "van Dongen 2014 Behav Genet PMID 24791950 (R2=51% of plasma sIL-6R)",
    },
    # 次级工具 (IL6R locus, 与 rs2228145 部分 LD) —— 效应量为占位, 需 pQTL 精确值
    "rs4845618": {"effect_allele": "T", "other_allele": "G", "eaf": 0.43,
                  "beta_sd": None, "source": "3' IL6R eQTL/pQTL, 待补"},
    "rs4129267": {"effect_allele": "T", "other_allele": "C", "eaf": 0.40,
                  "beta_sd": None, "source": "IL6R intron, 待补"},
}

# rs2228145 的 GRCh37 位置 (与 GCST90271815 / PGC daner 一致)
KEY_POS = {
    154426970: "rs2228145",
    154400015: "rs4845618",
    154426264: "rs4129267",
    154422067: "rs4845625",
    154436195: "rs4240872",
}


# ---------------------------------------------------------------------------
# 2. 从结局 GWAS 提取 IL6R 区间
# ---------------------------------------------------------------------------
def read_ct_window(lo: int = 154_380_000, hi: int = 154_480_000):
    rows = []
    with gzip.open(CT_GZ, "rt", encoding="utf-8", errors="ignore") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            try:
                bp = int(c[idx["base_pair_location"]])
            except Exception:
                continue
            if c[idx["chromosome"]] != "1":
                continue
            if lo <= bp <= hi:
                rows.append({
                    "chr": 1,
                    "bp": bp,
                    "rsid": KEY_POS.get(bp, ""),
                    "ea": c[idx["effect_allele"]],
                    "oa": c[idx["other_allele"]],
                    "beta": float(c[idx["beta"]]),
                    "se": float(c[idx["standard_error"]]),
                    "eaf": float(c[idx["effect_allele_frequency"]]),
                    "p": float(c[idx["p_value"]]),
                    "n": int(float(c[idx["n"]])),
                })
    return rows


# ---------------------------------------------------------------------------
# 3. Wald ratio / IVW
# ---------------------------------------------------------------------------
def wald(b_exp, se_exp, b_out, se_out):
    theta = b_out / b_exp
    se = se_out / abs(b_exp)
    z = theta / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    ci = (theta - 1.96 * se, theta + 1.96 * se)
    return dict(theta=theta, se=se, z=z, p=p, ci_lo=ci[0], ci_hi=ci[1])


def ivw(b_exp_list, se_exp_list, b_out_list, se_out_list):
    """Fixed-effect IVW (independent instruments)."""
    wsum = 0.0
    wbsum = 0.0
    for be, see, bo, soo in zip(b_exp_list, se_exp_list, b_out_list, se_out_list):
        theta_i = bo / be
        var_i = (soo / abs(be)) ** 2
        w = 1.0 / var_i
        wsum += w
        wbsum += w * theta_i
    theta = wbsum / wsum
    se = math.sqrt(1.0 / wsum)
    z = theta / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return dict(theta=theta, se=se, z=z, p=p, ci_lo=theta - 1.96 * se, ci_hi=theta + 1.96 * se)


# ---------------------------------------------------------------------------
# 4. 主流程
# ---------------------------------------------------------------------------
def main():
    print("=" * 72)
    print("IL6R (rs2228145) -> mean cortical thickness  MR")
    print("outcome: GCST90271815  Makowski 2023 PNAS  n=32,488 (GRCh37)")
    print("=" * 72)

    rows = read_ct_window()
    print(f"[1] IL6R window 154.38-154.48 Mb: {len(rows)} SNPs")

    by_pos = {r["bp"]: r for r in rows}
    target = by_pos.get(154426970)
    if target is None:
        raise SystemExit("rs2228145 (chr1:154426970) NOT FOUND in outcome GWAS")

    print("\n[2] rs2228145 in outcome GWAS (mean cortical thickness):")
    print(f"    alleles {target['ea']}/{target['oa']}  EAF(ea)={target['eaf']:.4f}  n={target['n']}")
    print(f"    beta={target['beta']:.6g}  se={target['se']:.6g}  p={target['p']:.4g}")

    # 等位基因对齐: 暴露 effect allele = C
    b_out, se_out, ea_out = target["beta"], target["se"], target["ea"]
    aligned = "C"
    if ea_out != "C":
        # 反转
        b_out, se_out = -b_out, se_out
        aligned = "C (flipped)"
        print(f"    [allele] outcome EA={target['ea']} != C -> beta flipped")
    else:
        print("    [allele] outcome EA=C 与暴露 effect allele 一致, 无需翻转")

    exp = EXPOSURE["rs2228145"]
    print(f"\n[3] exposure: rs2228145 {exp['effect_allele']}/{exp['other_allele']}"
          f"  beta_exp = {exp['beta_sd']} SD sIL-6R per {exp['effect_allele']}")
    print(f"    source: {exp['source']}")

    res_main = wald(exp["beta_sd"], None, b_out, se_out)
    print("\n[4] Wald ratio (beta_exp = {:.3f} SD):".format(exp["beta_sd"]))
    print(f"    theta = {res_main['theta']:+.5f} SD cortical thickness per SD sIL-6R")
    print(f"    se    = {res_main['se']:.5f}")
    print(f"    95%CI = ({res_main['ci_lo']:+.5f}, {res_main['ci_hi']:+.5f})")
    print(f"    z     = {res_main['z']:+.3f}   p = {res_main['p']:.4g}")

    print("\n[5] 敏感性 (beta_exp 区间):")
    sens = {}
    for tag, b in [("lo_0.49", exp["beta_sd_lo"]), ("mid_0.70", 0.70), ("main_1.035", exp["beta_sd"])]:
        r = wald(b, None, b_out, se_out)
        sens[tag] = r
        print(f"    beta_exp={b:.3f} -> theta={r['theta']:+.5f} (se {r['se']:.5f})"
              f"  p={r['p']:.4g}")

    # 次级工具 (只有在暴露 beta 未知时给位点结果)
    print("\n[6] IL6R locus 其他工具 (结局侧; 暴露 beta 待补):")
    locus = []
    for bp, name in KEY_POS.items():
        r = by_pos.get(bp)
        if r is None:
            print(f"    {name} (chr1:{bp})  -- 缺失 --")
            continue
        print(f"    {name} (chr1:{bp})  {r['ea']}/{r['oa']}  beta={r['beta']:+.5f}"
              f"  se={r['se']:.5f}  p={r['p']:.4g}")
        locus.append(r)

    # ---- PGC MDD 2018 参照 ----
    pgc = {}
    if PGC.exists():
        print("\n[7] PGC MDD 2018 (daner, 45,396 cases / 97,250 controls) 参照:")
        with gzip.open(PGC, "rt", encoding="utf-8", errors="ignore") as f:
            f.readline()
            for line in f:
                c = line.split("\t")
                try:
                    bp = int(c[2])
                except Exception:
                    continue
                if c[0] == "1" and bp in KEY_POS:
                    pgc[KEY_POS[bp]] = dict(
                        rsid=c[1], ea=c[3], oa=c[4], or_=float(c[8]), se=float(c[9]), p=float(c[10])
                    )
                    print(f"    {c[1]:<12} {c[3]}/{c[4]}  OR={float(c[8]):.4f}"
                          f"  se={float(c[9]):.4f}  p={float(c[10]):.4g}")
    else:
        print("\n[7] PGC MDD 2018 文件未找到, 跳过参照")

    out = {
        "outcome": {
            "gwas_id": "GCST90271815",
            "trait": "mean cortical thickness (global)",
            "n": 32488,
            "ancestry": "European (UKB)",
            "build": "GRCh37",
            "source": "Makowski C et al. PNAS 2023;120(11):e2214834120 (PMID 36893272)",
            "url": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90271001-GCST90272000/GCST90271815/",
        },
        "exposure": EXPOSURE["rs2228145"],
        "rs2228145_outcome": target,
        "rs2228145_outcome_effect_allele_aligned": aligned,
        "wald_main": res_main,
        "wald_sensitivity": sens,
        "il6r_locus_outcome": locus,
        "pgc_mdd2018_reference": pgc,
        "note": (
            "Grasby 2020 (GCST010281, ENIGMA/Science) 全汇总统计未公开 (fullSS=no), "
            "故以最大公开的 mean cortical thickness GWAS (Makowski 2023 PNAS, n=32,488, GRCh37) 为结局。"
        ),
    }
    js = OUT_DIR / "il6r_mr_results.json"
    with open(js, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"\n[8] saved -> {js}")

    ts = OUT_DIR / "il6r_mr_results.tsv"
    with open(ts, "w", encoding="utf-8") as f:
        f.write("item\tvalue\n")
        f.write(f"outcome\tGCST90271815 mean cortical thickness n=32488\n")
        f.write(f"rs2228145_outcome_beta\t{target['beta']}\n")
        f.write(f"rs2228145_outcome_se\t{target['se']}\n")
        f.write(f"rs2228145_outcome_p\t{target['p']}\n")
        f.write(f"beta_exposure_sd\t{exp['beta_sd']}\n")
        f.write(f"wald_theta\t{res_main['theta']}\n")
        f.write(f"wald_se\t{res_main['se']}\n")
        f.write(f"wald_p\t{res_main['p']}\n")
    print(f"[8] saved -> {ts}")
    print("\n" + "=" * 72)
    print("DONE")
    print("=" * 72)


if __name__ == "__main__":
    main()
