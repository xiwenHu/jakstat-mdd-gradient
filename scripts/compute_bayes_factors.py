# -*- coding: utf-8 -*-
"""
Bayes factor calculations for Path-1 manuscript.
A. BIC-approximation BF01 (Wagenmakers 2007) for correlations and two-sample t-tests
B. JZS BF10 (Rouder et al. 2009) numerical integration for key t-tests
C. Wakefield approximate Bayes factor (ABF) for MR / GWAS associations

Run:  python compute_bayes_factors.py
"""
import math
import json
import pandas as pd
import numpy as np
from scipy.integrate import quad

OUT = r"C:\Users\ThinkPad\WorkBuddy\2026-10-20-11-06-57\outputs\data\bayes"
import os
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- A. BIC BF01

def bic_bf01_corr(r, n):
    """BIC approximation for testing rho=0. df = n-2. Fisher z gives t."""
    df = n - 2
    if df <= 2:
        return float("nan")
    # t equivalent of correlation
    t = r * math.sqrt(df / (1.0 - r * r)) if abs(r) < 1 else 1e6
    # two-model BIC with RSS0/RSS1 ratio (1 - r^2) vs 1
    # BIC1-BIC0 = n*ln(1-r^2) + ln(n)   [model 1 adds the predictor]
    bic1_minus_bic0 = n * math.log(max(1e-300, 1.0 - r * r)) + math.log(n)
    return math.exp(0.5 * bic1_minus_bic0)


def bic_bf01_t(t, n1, n2):
    """BIC approximation for two-sample t-test (Welch ~ pooled)."""
    N = n1 + n2
    df = N - 2
    bic1_minus_bic0 = N * math.log(max(1e-300, (df) / (df + t * t))) + math.log(N)
    return math.exp(0.5 * bic1_minus_bic0)


# ---------------------------------------------------------------- B. JZS BF10

def _inv_chi2_pdf(g, nu, scale):
    if g <= 0:
        return 0.0
    logp = (0.5 * nu) * math.log(scale) - (0.5 * nu) * math.log(2.0) \
           - math.lgamma(0.5 * nu) - (0.5 * nu + 1.0) * math.log(g) \
           - scale / (2.0 * g)
    return math.exp(logp)


def jzs_bf10_t(t, n, r=0.707):
    """Rouder et al. 2009 JZS BF10 for one-sample / paired t (df = n-1).
    For two-sample with n1,n2 use n = n1*n2/(n1+n2)*2 equivalent N_eff.
    """
    df = n - 1
    if df <= 1:
        return float("nan")

    def integrand(g):
        return ((1.0 + n * g) ** -0.5
                * (1.0 + t * t / ((1.0 + n * g) * df)) ** (-(df + 1.0) / 2.0)
                * _inv_chi2_pdf(g, 1, r * r))

    try:
        num, _ = quad(integrand, 1e-10, 50.0, limit=200)
    except Exception:
        num = np.nan
    den = (1.0 + t * t / df) ** (-(df + 1.0) / 2.0)
    return num / den


def jzs_bf10_twosample(t, n1, n2, r=0.707):
    """JZS for independent two-sample t: use N = n1+n2 as in BayesFactor pkg default."""
    return jzs_bf10_t(t, n1 + n2, r=r)


# ---------------------------------------------------------------- C. Wakefield ABF

def wakefield_bf01(beta, se, W):
    """Wakefield 2009 approximate Bayes factor. Returns BF01 (favouring null).
    log10(ABF) = 0.5*log10(V/(V+W)) + 0.5*z^2*W/(V+W)/ln(10)
    """
    V = se * se
    z = beta / se
    log_abf = 0.5 * math.log(V / (V + W)) + 0.5 * z * z * W / (V + W)
    return math.exp(-log_abf)   # BF01


# ---------------------------------------------------------------- run
results = {}

# --- 1. spatial correlations (n = 339 regions)
n_spatial = 339
spatial = {
    "raw_IL6R_rho_0.315": {"r": 0.315, "n": n_spatial},
    "partial_PC1_5_rho_0.005": {"r": 0.005, "n": n_spatial},
    "partial_PC1_only_rho_0.179": {"r": 0.179, "n": n_spatial},
}
sp_out = {}
for k, v in spatial.items():
    bf01 = bic_bf01_corr(v["r"], v["n"])
    t = v["r"] * math.sqrt((v["n"] - 2) / max(1e-12, 1 - v["r"] ** 2))
    bf10_jzs = jzs_bf10_t(abs(t), v["n"])
    sp_out[k] = {"r": v["r"], "n": v["n"], "t_equiv": round(t, 3),
                 "BF01_BIC": round(bf01, 2), "BF10_JZS": round(bf10_jzs, 2),
                 "BF01_JZS": round(1.0 / bf10_jzs, 2) if bf10_jzs == bf10_jzs else None}
results["spatial"] = sp_out

# --- 2. MR (Wald ratio on global mean CT) + PGC rs2228145
mr = {}
# outcome scale (SD per allele); z = -1.4931
beta_out, se_out = -0.0119313, 0.00799077
for W in (0.01 ** 2, 0.02 ** 2, 0.04 ** 2, 0.08 ** 2):
    mr[f"MR_globalCT_W_{math.sqrt(W):.3f}"] = {
        "beta": beta_out, "se": se_out,
        "W": W, "BF01_Wakefield": round(wakefield_bf01(beta_out, se_out, W), 3)
    }
# PGC MDD 2018 rs2228145: OR=0.989 per A, p=0.230 -> derive se
or_pgc, p_pgc = 0.989, 0.230
beta_pgc = math.log(or_pgc)
# two-sided p -> |z|
from math import erf, sqrt as _sqrt
# invert normal cdf
def p_to_z(p):
    # bisection on two-sided
    lo, hi = 0.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        tail = 2 * (1 - 0.5 * (1 + erf(mid / _sqrt(2))))
        if tail > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2
z_pgc = p_to_z(p_pgc)
se_pgc = abs(beta_pgc) / z_pgc
for W in (0.02 ** 2, 0.05 ** 2, 0.10 ** 2):
    mr[f"PGC_MDD_rs2228145_W_{math.sqrt(W):.2f}"] = {
        "beta_logOR": round(beta_pgc, 5), "se": round(se_pgc, 5),
        "W": W, "BF01_Wakefield": round(wakefield_bf01(beta_pgc, se_pgc, W), 3)
    }
results["mr_gwas"] = mr

# --- 3. postmortem DE (GSE53987): gene-level representative probes
tsv = r"C:\Users\ThinkPad\WorkBuddy\2026-10-20-11-06-57\outputs\data\gse53987\gse53987_de_results.tsv"
df = pd.read_csv(tsv, sep="\t")
best = df[df["is_best_probe"] == True].copy()
de_out = []
for _, row in best.iterrows():
    n1, n2 = int(row["n_case"]), int(row["n_ctrl"])
    t = float(row["t"])
    bf01_bic = bic_bf01_t(t, n1, n2)
    bf10_jzs = jzs_bf10_twosample(abs(t), n1, n2)
    de_out.append({
        "region": row["region"], "contrast": row["contrast"], "gene": row["gene"],
        "n_case": n1, "n_ctrl": n2, "t": round(t, 3),
        "p_ttest": float(row["p_ttest"]),
        "BF01_BIC": round(bf01_bic, 3),
        "BF10_JZS": round(bf10_jzs, 3),
        "BF01_JZS": round(1.0 / bf10_jzs, 3) if bf10_jzs == bf10_jzs else None,
    })
de = pd.DataFrame(de_out)
de.to_csv(os.path.join(OUT, "de_bf_all.csv"), index=False)

mdd = de[de["contrast"] == "MDD_vs_HC"]
scz = de[de["contrast"] == "SCZ_vs_HC"]
summary_de = {
    "MDD_vs_HC": {
        "n_tests": int(len(mdd)),
        "BF01_JZS_median": round(float(mdd["BF01_JZS"].median()), 2),
        "BF01_JZS_range": [round(float(mdd["BF01_JZS"].min()), 2),
                           round(float(mdd["BF01_JZS"].max()), 2)],
        "n_BF01_ge_3": int((mdd["BF01_JZS"] >= 3).sum()),
        "n_BF01_ge_10": int((mdd["BF01_JZS"] >= 10).sum()),
        "max_BF10_JZS": round(float(mdd["BF10_JZS"].max()), 2),
        "gene_max_BF10": mdd.loc[mdd["BF10_JZS"].idxmax(), ["region", "gene", "t"]].to_dict(),
    },
    "SCZ_vs_HC": {
        "n_tests": int(len(scz)),
        "n_BF10_JZS_ge_3": int((scz["BF10_JZS"] >= 3).sum()),
        "max_BF10_JZS": round(float(scz["BF10_JZS"].max()), 2),
        "top3": scz.nlargest(3, "BF10_JZS")[["region", "gene", "t", "BF10_JZS"]].to_dict("records"),
    },
}
results["de_summary"] = summary_de

# --- 4. MR panel direction consistency (12 regions) — report BF for strongest region
# ventral frontal: beta = -0.0203, se from CI ( -0.0354, -0.0051 ) -> se ~ (0.0354-0.0051)/(2*1.96)
se_vf = (0.0354 - 0.0051) / (2 * 1.96)
mr_panel = {"ventralFrontal": {"beta": -0.0203, "se": round(se_vf, 5)}}
for W in (0.01 ** 2, 0.02 ** 2, 0.04 ** 2):
    mr_panel["ventralFrontal"][f"BF01_W_{math.sqrt(W):.3f}"] = round(
        wakefield_bf01(-0.0203, se_vf, W), 3)
results["mr_panel"] = mr_panel

with open(os.path.join(OUT, "bayes_results.json"), "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)

print(json.dumps(results, indent=2, ensure_ascii=False))
print("\n[written]", os.path.join(OUT, "bayes_results.json"))
print("[written]", os.path.join(OUT, "de_bf_all.csv"))
