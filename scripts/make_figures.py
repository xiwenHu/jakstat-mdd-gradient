# -*- coding: utf-8 -*-
"""
Figures 1-3 for Path-1 manuscript (English labels only, 300 dpi, PNG+PDF+TIFF).
Run:  python make_figures.py
"""
import warnings, os, math, json
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import spearmanr, rankdata

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

ROOT = r"C:\Users\ThinkPad\WorkBuddy\2026-10-20-11-06-57\outputs"
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
BASE_F = r"F:\炎症因子文章"

C_MAIN = "#1F4E79"      # deep blue
C_ACC = "#C0504D"       # muted red
C_NULL = "#9E9E9E"      # gray
C_HI = "#E8A33D"        # amber
C_SCZ = "#6A3D9A"       # purple
C_MDD = "#1F4E79"

GENES = ['IL6R','IL6ST','JAK1','JAK2','STAT3','STAT4','STAT1','TYK2','STAT5A','STAT5B','STAT6','JAK3','SOCS1']


def save(fig, name):
    for ext in ("png", "pdf", "tiff"):
        p = os.path.join(FIG, f"{name}.{ext}")
        fig.savefig(p, dpi=300)
    plt.close(fig)
    print("[fig]", name)


# =====================================================================
# Load / recompute core arrays
# =====================================================================
sign = pd.read_csv(os.path.join(ROOT, "sign_check", "sign_check_13genes.csv"))
part = pd.read_csv(os.path.join(ROOT, "sign_check", "partial_control_PC1.csv"))
perparcel = pd.read_csv(os.path.join(ROOT, "sign_check", "sign_check_per_parcel.csv"))

# recompute PC1 and PC1-5 partials from full expression matrix
d = pd.read_csv(os.path.join(BASE_F, "traits_schaefer400", "trait_mddadult_case-controls_CortThick.csv"))
col = [c for c in d.columns if c != "parcel"][0]
d["parcel"] = d["parcel"].astype(int)
expr = pd.read_csv(os.path.join(BASE_F, "expr_schaefer400_allgenes.csv"), index_col=0)
expr.index = expr.index.astype(int)
expr = expr.loc[expr.notna().all(axis=1)]
expr = expr.loc[sorted(expr.index)]
common = sorted(set(expr.index) & set(d["parcel"]))
expr = expr.loc[common]
dv = d.set_index("parcel").loc[common, col].to_numpy(float)

X = expr.to_numpy(float)
Xz = (X - X.mean(0)) / X.std(0)
U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
pcs = U[:, :5] * S[:5]
pc1 = pcs[:, 0].copy()
pc1 /= pc1.std()
var1 = S[0] ** 2 / (S ** 2).sum()
r_pc1_d, p_pc1_d = spearmanr(pc1, dv)


def partial_spearman(x, y, Z):
    """Rank-based partial correlation, matching gradient_confound_check.py:
    rank x,y and each column of Z; residualise ranks; Pearson r of residuals."""
    from scipy.stats import pearsonr
    rx, ry = rankdata(x), rankdata(y)
    rZ = np.apply_along_axis(rankdata, 0, Z.reshape(len(x), -1) if Z.ndim == 1 else Z)
    A = np.hstack([rZ, np.ones((len(ry), 1))])
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None)
    by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    ex, ey = rx - A @ bx, ry - A @ by
    return pearsonr(ex, ey)


# per-gene: rho(d), partial PC1, partial PC1-5
rows = []
for g in GENES:
    gv = expr[g].to_numpy(float)
    r_d, p_d = spearmanr(gv, dv)
    r_abs, _ = spearmanr(gv, np.abs(dv))
    r_p1, p_p1 = partial_spearman(gv, dv, pcs[:, :1])
    r_p5, p_p5 = partial_spearman(gv, dv, pcs[:, :5])
    rows.append({"gene": g, "rho_d": r_d, "p_d": p_d, "rho_absd": r_abs,
                 "rho_pc1": r_p1, "p_pc1": p_p1,
                 "rho_pc5": r_p5, "p_pc5": p_p5,
                 "drop1": 100 * (1 - abs(r_p1) / abs(r_d)) if abs(r_d) > 1e-12 else np.nan,
                 "drop5": 100 * (1 - abs(r_p5) / abs(r_d)) if abs(r_d) > 1e-12 else np.nan})
gdf = pd.DataFrame(rows)
gdf.to_csv(os.path.join(ROOT, "figures", "fig1C_gradient_partial_data.csv"), index=False)
print(gdf[["gene", "rho_d", "rho_absd", "rho_pc1", "rho_pc5", "drop5"]].to_string(index=False))

# null distributions (test4)
t4 = os.path.join(ROOT, "test4")
count_null = np.load(os.path.join(t4, "partA_count_null_matched.npy"))
meanZ_null = np.load(os.path.join(t4, "partA_meanZ_null_matched.npy"))
spatial_null = pd.read_csv(os.path.join(t4, "test4b_spatial_geneset_null.csv"))

# MR panel
mr_panel = json.load(open(os.path.join(ROOT, "data", "mr", "il6r_mr_panel_results.json")))["panel"]
mr_main = json.load(open(os.path.join(ROOT, "data", "mr", "il6r_mr_results.json")))

# DE
de = pd.read_csv(os.path.join(ROOT, "data", "gse53987", "gse53987_de_results.tsv"), sep="\t")
best = de[de["is_best_probe"] == True].copy()

# Bayes
bayes = json.load(open(os.path.join(ROOT, "data", "bayes", "bayes_results.json")))

MDE_MR = 0.0224
BONF = 0.05 / 12

# =====================================================================
# FIGURE 1 — spatial signal and its dissolution
# =====================================================================
fig = plt.figure(figsize=(7.2, 6.0))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.32)

# ---- A: rho(d) vs rho(|d|)
axA = fig.add_subplot(gs[0, 0])
axA.axvline(0, color="k", lw=0.6)
axA.axhline(0, color="k", lw=0.6)
axA.scatter(gdf["rho_d"], gdf["rho_absd"], s=36, color=C_MAIN, edgecolor="white", linewidth=0.6, zorder=3)
lab_off = {"IL6R": (6, -2), "STAT4": (6, 6), "STAT1": (6, 4), "TYK2": (6, 4)}
for g, off in lab_off.items():
    row = gdf[gdf["gene"] == g].iloc[0]
    axA.annotate(g, (row["rho_d"], row["rho_absd"]),
                 textcoords="offset points", xytext=off, fontsize=7)
axA.set_xlabel(r"Spearman $\rho$ with signed effect sizes $\rho(d)$")
axA.set_ylabel(r"Spearman $\rho$ with $|d|$  $\rho(|d|)$")
axA.set_title("A  Sign structure of 13 genes", loc="left", fontweight="bold", fontsize=9)
axA.text(0.02, 0.04, "associated genes: $\\rho(d)>0$ and $\\rho(|d|)<0$",
         transform=axA.transAxes, fontsize=6.5, color="#444444")

# ---- B: PC1 vs d
axB = fig.add_subplot(gs[0, 1])
sc = axB.scatter(pc1, dv, s=14, c=dv, cmap="coolwarm", alpha=0.85, linewidths=0)
# binned trend line (deciles of PC1)
bins = np.quantile(pc1, np.linspace(0, 1, 11))
idx = np.digitize(pc1, bins[1:-1])
bm = [dv[idx == k].mean() for k in range(idx.max() + 1) if (idx == k).sum() > 2]
bx = [pc1[idx == k].mean() for k in range(idx.max() + 1) if (idx == k).sum() > 2]
axB.plot(bx, bm, color="k", lw=1.6, marker="o", ms=3, zorder=4)
axB.set_xlabel("PC1 (sensory-to-association gradient)")
axB.set_ylabel("MDD cortical thickness Cohen's $d$")
axB.set_title(f"B  Gradient prior: $\\rho$ = {r_pc1_d:.3f}, P = {p_pc1_d:.1e}", loc="left", fontweight="bold", fontsize=9)
cb = fig.colorbar(sc, ax=axB, fraction=0.05, pad=0.02)
cb.set_label("Cohen's $d$", fontsize=7)

# ---- C: attenuation
axC = fig.add_subplot(gs[1, 0])
order = gdf.reindex(gdf["rho_d"].abs().sort_values(ascending=False).index)
xpos = [0, 1, 2]
for _, r in order.iterrows():
    if abs(r["rho_d"]) < 0.05 and r["gene"] != "IL6R":
        continue  # skip near-zero genes whose partials wander
    ys = [r["rho_d"], r["rho_pc1"], r["rho_pc5"]]
    col = C_ACC if r["gene"] == "IL6R" else "#BBBBBB"
    lw = 2.0 if r["gene"] == "IL6R" else 0.9
    axC.plot(xpos, ys, color=col, lw=lw, marker="o", ms=3 if r["gene"] == "IL6R" else 2, zorder=3)
    if r["gene"] == "IL6R":
        for x, y in zip(xpos, ys):
            dy = -12 if x == 2 else 8
            axC.annotate(f"{y:+.3f}", (x, y), textcoords="offset points",
                         xytext=(0, dy), fontsize=7, ha="center", color=C_ACC, fontweight="bold")
        il6r_drop = r["drop5"]
axC.axhline(0, color="k", lw=0.6)
axC.set_xticks(xpos)
axC.set_xticklabels([r"raw $\rho(d)$", "partial | PC1", "partial | PC1–PC5"])
axC.set_ylabel("Spearman $\\rho$ with $d$")
axC.set_title("C  Gradient control dissolves the association", loc="left", fontweight="bold", fontsize=9)
axC.text(2, 0.30, f"IL6R: {il6r_drop:.1f}% attenuation\n$BF_{{01}}$ = {bayes['spatial']['partial_PC1_5_rho_0.005']['BF01_JZS']:.1f}",
         fontsize=7, ha="right", color=C_ACC)

# ---- D: matched nulls
axD = fig.add_subplot(gs[1, 1])
axD.hist(count_null, bins=np.arange(-0.5, 7.5, 1), color=C_NULL, edgecolor="white",
         linewidth=0.6, density=True, alpha=0.9)
axD.axvline(3, color=C_ACC, lw=2.0)
axD.set_ylim(0, 0.42)
axD.annotate("observed = 3", (2.85, 0.37), color=C_ACC,
             fontsize=7, ha="right", fontweight="bold")
axD.set_xlabel("genes at nominal $P<0.05$ per matched set")
axD.set_ylabel("density (10,000 matched sets)")
axD.set_title("D  Matched gene-set nulls", loc="left", fontweight="bold", fontsize=9)
txt = ("empirical P:\n"
       " nominal excess  0.116\n"
       " gene-set mean|$\\rho$|  0.985\n"
       " grid top-18 $\\Sigma|\\rho|$  0.978\n"
       " MAGMA competitive  0.744")
axD.text(0.97, 0.78, txt, transform=axD.transAxes, fontsize=6.5, va="top", ha="right",
         bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#CCCCCC"))

save(fig, "Fig1_spatial_dissolution")

# =====================================================================
# FIGURE 2 — genetic evidence
# =====================================================================
fig = plt.figure(figsize=(7.2, 6.0))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.45)

# ---- A: MR forest (Wald theta on SD/allele scale)
axA = fig.add_subplot(gs[0, :])
items = []
# global
w = mr_main["wald_main"]
items.append(("Global mean CT (primary)", w["theta"], w["ci_lo"], w["ci_hi"], True, "thickness"))
for r in mr_panel:
    if r["label"] == "Global mean CT":
        continue  # already plotted as the primary row
    theta = r["theta"]
    se = r["se"] / abs(mr_main["exposure"]["beta_sd"])  # panel se is on outcome scale
    # panel entries store outcome beta/se; theta = beta/beta_exp
    bexp = mr_main["exposure"]["beta_sd"]
    th = r["beta"] / bexp
    se_th = r["se"] / abs(bexp)
    items.append((r["label"], th, th - 1.96 * se_th, th + 1.96 * se_th, False, r["cat"]))
labels = [i[0] for i in items][::-1]
ths = [i[1] for i in items][::-1]
los = [i[2] for i in items][::-1]
his = [i[3] for i in items][::-1]
prim = [i[4] for i in items][::-1]
cats = [i[5] for i in items][::-1]
y = np.arange(len(labels))
for i in range(len(labels)):
    col = C_ACC if prim[i] else ("#4B7BA8" if cats[i] == "thickness" else "#7F9DB9")
    axA.plot([los[i], his[i]], [y[i], y[i]], color=col, lw=2.0 if prim[i] else 1.3)
    axA.plot(ths[i], y[i], "s" if prim[i] else "o", color=col, ms=5 if prim[i] else 4, zorder=3)
axA.axvline(0, color="k", lw=0.8)
axA.axvline(-MDE_MR, color=C_HI, lw=1.2, ls="--")
axA.axvline(MDE_MR, color=C_HI, lw=1.2, ls="--")
axA.text(-MDE_MR, len(labels) - 0.2, "MDE", color=C_HI, fontsize=6.5, ha="center")
axA.set_yticks(y)
axA.set_yticklabels(labels, fontsize=7)
axA.set_xlabel("Wald ratio $\\theta$ (SD per C allele of rs2228145)")
axA.set_title("A  cis-MR of IL6R signalling on cortical thickness (all estimates negative; global primary endpoint $P$ = 0.135)",
              loc="left", fontweight="bold", fontsize=8.5)

# ---- B: genetic triangulation numbers
axB = fig.add_subplot(gs[1, 0])
axB.axis("off")
tbl = [
    ["test", "estimate", "P (emp.)"],
    ["MAGMA competitive", "Z = 0.43", "0.744"],
    ["Nominal 3/13 (O/E 4.6)", "1.23 ± 1.04", "0.116"],
    ["PGC rs2228145", "OR = 0.989", "0.230"],
    ["MR global CT", "θ = −0.012", "0.135"],
    ["MR surface area", "θ = +0.011", "0.166"],
]
colx = [0.02, 0.58, 0.86]
for i, row in enumerate(tbl):
    for j, cell in enumerate(row):
        axB.text(colx[j], 0.94 - i * 0.155, cell, fontsize=6.8,
                 fontweight="bold" if (i == 0 or j == 0) else "normal",
                 transform=axB.transAxes, family="DejaVu Sans")
axB.plot([0.0, 1.0], [0.825, 0.825], color="#CCCCCC", lw=0.8, transform=axB.transAxes)
axB.set_title("B  Genetic triangulation: all null", loc="left", fontweight="bold", fontsize=9)

# ---- C: MDE vs observed
axC = fig.add_subplot(gs[1, 1])
names = ["MR global\nCT", "MR ventral\nfrontal", "DE MDD\nmax", "DE SCZ\nmax"]
mde_plot = [MDE_MR, MDE_MR, 0.12, 0.12]
obs_plot = [abs(mr_main["wald_main"]["theta"]), 0.02027, 2 ** 0.139 - 1, 2 ** 0.313 - 1]
x = np.arange(len(names))
axC.bar(x - 0.18, mde_plot, 0.34, color=C_NULL, label="MDE (80% power)")
axC.bar(x + 0.18, obs_plot, 0.34, color=[C_ACC, C_ACC, C_MDD, C_SCZ], label="observed")
for xi, o, m in zip(x, obs_plot, mde_plot):
    axC.text(xi + 0.18, o + 0.006, f"{100*o/m:.0f}%", fontsize=6.5, ha="center")
axC.set_xticks(x)
axC.set_xticklabels(names, fontsize=7)
axC.set_ylabel("effect size (SD / fold-change $-$ 1)")
axC.set_title("C  Observed effects vs minimum detectable effect", loc="left", fontweight="bold", fontsize=9)
axC.legend(fontsize=6.5, frameon=False, loc="upper left")

save(fig, "Fig2_genetic_evidence")

# =====================================================================
# FIGURE 3 — postmortem brain with positive control
# =====================================================================
fig = plt.figure(figsize=(7.2, 6.0))
gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.35)


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def volcano(ax, df, title, color, hl_fdr=False, hl_top=()):
    df = df.copy()
    df["fdr_plot"] = bh_fdr(df["p_ttest"].to_numpy(float))
    x = df["logFC"].to_numpy(float)
    y = -np.log10(df["p_ttest"].to_numpy(float))
    sig = df["fdr_plot"].to_numpy() < 0.05
    ax.scatter(x[~sig], y[~sig], s=22, color=color, alpha=0.75,
               edgecolor="white", linewidth=0.5)
    if sig.any():
        ax.scatter(x[sig], y[sig], s=42, color=C_ACC,
                   edgecolor="k", linewidth=0.6, zorder=4)
    ax.axhline(-np.log10(0.05), color="k", lw=0.7, ls=":")
    ax.set_xlabel(r"log$_2$ fold-change (case $-$ control)")
    ax.set_ylabel(r"$-\log_{10}P$")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=9)
    # label FDR-significant points (one row each = gene x region)
    if hl_fdr:
        for _, r in df[sig].iterrows():
            ax.annotate(f"{r['gene']} {r['region']}", (r["logFC"], -np.log10(r["p_ttest"])),
                        textcoords="offset points", xytext=(5, 2), fontsize=6.5,
                        fontweight="bold")
    # label the single top point of named genes
    for g in hl_top:
        sub = df[df["gene"] == g]
        if len(sub) == 0:
            continue
        r = sub.loc[sub["p_ttest"].idxmin()]
        ax.annotate(f"{g} {r['region']}", (r["logFC"], -np.log10(r["p_ttest"])),
                    textcoords="offset points", xytext=(5, 2), fontsize=6.5,
                    fontweight="bold")


axA = fig.add_subplot(gs[0, 0])
mdd = best[best["contrast"] == "MDD_vs_HC"]
volcano(axA, mdd, "A  MDD vs HC: 0/39 FDR", C_MDD, hl_fdr=False, hl_top=("STAT1",))
axA.text(0.98, 0.95, "max |log$_2$FC| = 0.139 (1.10-fold)\nnominal 2/39 (expected 1.95)",
         transform=axA.transAxes, fontsize=6.5, ha="right", va="top",
         bbox=dict(boxstyle="round,pad=0.3", fc="#F7F7F7", ec="#CCCCCC"))

axB = fig.add_subplot(gs[0, 1])
scz = best[best["contrast"] == "SCZ_vs_HC"]
volcano(axB, scz, "B  SCZ vs HC (positive control): 3/39 FDR", C_SCZ, hl_fdr=True)

axC = fig.add_subplot(gs[1, 0])
# MDE panel per region vs max observed
regions = ["PFC", "HIP", "STR"]
mde_r = [1.12, 1.16, 1.12]
maxobs = []
for rg in regions:
    sub = mdd[mdd["region"] == rg]
    maxobs.append(1 + sub["logFC"].abs().max())  # fold scale
sczmax = []
for rg in regions:
    sub = scz[scz["region"] == rg]
    sczmax.append(1 + sub["logFC"].abs().max())
x = np.arange(3)
axC.bar(x - 0.2, [v - 1 for v in mde_r], 0.38, color=C_NULL, label="MDE")
axC.bar(x + 0.2, [v - 1 for v in maxobs], 0.38, color=C_MDD, label="MDD max |FC|")
axC.plot(x, [v - 1 for v in sczmax], "D", color=C_SCZ, ms=6, label="SCZ max |FC|", zorder=4)
axC.axhline(0, color="k", lw=0.6)
axC.set_xticks(x)
axC.set_xticklabels(["PFC", "Hippocampus", "Striatum"])
axC.set_ylabel("fold-change $-$ 1")
axC.set_title("C  MDE vs observed", loc="left", fontweight="bold", fontsize=9)
axC.text(0.02, 0.62, "MDD max < MDE", transform=axC.transAxes,
         fontsize=6.5, va="top", color="#444444")
axC.legend(fontsize=6.5, frameon=False, loc="center right")

axD = fig.add_subplot(gs[1, 1])
up = [6, 4, 6]
down = [7, 9, 7]
axD.bar(x, up, 0.55, color="#C0504D", label="up in MDD")
axD.bar(x, down, 0.55, bottom=up, color="#4B7BA8", label="down in MDD")
for xi, u, dn in zip(x, up, down):
    axD.text(xi, u / 2, str(u), ha="center", va="center", fontsize=7, color="white")
    axD.text(xi, u + dn / 2, str(dn), ha="center", va="center", fontsize=7, color="white")
axD.axhline(6.5, color="k", lw=0.8, ls="--")
axD.text(2.35, 6.7, "chance level (6.5)", fontsize=6.5, ha="right")
axD.set_xticks(x)
axD.set_xticklabels(["PFC\n(sign $P$ = 1.00)", "HIP\n(0.27)", "STR\n(1.00)"], fontsize=7)
axD.set_ylabel("genes (of 13)")
axD.set_ylim(0, 16.5)
axD.set_title("D  No consistent direction of change", loc="left", fontweight="bold", fontsize=9)
axD.legend(fontsize=6.5, frameon=False, loc="upper center", ncol=2)

save(fig, "Fig3_postmortem_positive_control")

print("[done] all figures ->", FIG)
