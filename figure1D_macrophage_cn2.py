#!/usr/bin/env python3
"""Figure 1D: Macrophage abundance (left) and CN2 proportion (right)."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel

INPUT_OBS = "phenotyped-obs.csv"
INPUT_CN = "pheno_8clus-sample_pct.csv"
INPUT_CLIN = "nivo_clinical.csv"
OUTPUT = "Figure1D_macrophage_cn2"
COLOR_R, COLOR_NR = "#2980b9", "#e74c3c"

obs = pd.read_csv(INPUT_OBS)
clin = pd.read_csv(INPUT_CLIN)
cn_pct = pd.read_csv(INPUT_CN)

if "group" not in obs.columns:
    obs = obs.merge(clin[["pt", "group"]].drop_duplicates(), on="pt", how="left")
obs["pt_tp"] = obs["pt"].astype(str) + "_" + obs["tp"].astype(str)

# Macrophage %
counts = obs.groupby(["pt_tp", "pheno"]).size().reset_index(name="n")
totals = obs.groupby("pt_tp").size().reset_index(name="total")
counts = counts.merge(totals, on="pt_tp")
counts["pct"] = counts["n"] / counts["total"] * 100
meta = obs[["pt_tp", "pt", "tp", "group"]].drop_duplicates()
counts = counts.merge(meta, on="pt_tp")
mac = counts[counts.pheno == "Mac"].copy()
mac["x"] = mac["group"] + "_" + mac["tp"]

# CN2 %
cn_long = cn_pct.melt(id_vars=[c for c in cn_pct.columns if not c.startswith("CN")],
                       var_name="cn", value_name="cn_pct") if "CN0" in cn_pct.columns else cn_pct
if "cn" not in cn_long.columns:
    cn_long = cn_pct.melt(id_vars=["pt_tp"], var_name="cn", value_name="cn_pct")
if "group" not in cn_long.columns:
    pt_meta = obs[["pt_tp", "pt", "tp", "group"]].drop_duplicates()
    cn_long = cn_long.merge(pt_meta, on="pt_tp", how="left")
cn_agg = cn_long.groupby(["pt_tp", "cn"])["cn_pct"].mean().reset_index()
cn_agg = cn_agg.merge(meta, on="pt_tp")
cn2 = cn_agg[cn_agg.cn == "CN2"].copy()
cn2["x"] = cn2["group"] + "_" + cn2["tp"]

order = ["good_Baseline", "good_Post", "poor_Baseline", "poor_Post"]
labels = ["R\nBL", "R\nPost", "NR\nBL", "NR\nPost"]
colors = [COLOR_R, COLOR_R, COLOR_NR, COLOR_NR]

def plot_panel(ax, dat, ycol, ylabel, title):
    for i, (grp, lab, col) in enumerate(zip(order, labels, colors)):
        vals = dat.loc[dat.x == grp, ycol]
        bp = ax.boxplot([vals.values], positions=[i], widths=0.5, patch_artist=True, showfliers=False)
        bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha(0.6)
        bp["medians"][0].set_color("black")
        ax.scatter(np.full(len(vals), i) + np.random.default_rng(42).uniform(-0.1, 0.1, len(vals)),
                   vals, c=col, s=40, zorder=3, edgecolors="white", linewidths=0.5)
    for pt in dat.pt.unique():
        sub = dat[dat.pt == pt].set_index("x")
        for pre, post in [("good_Baseline", "good_Post"), ("poor_Baseline", "poor_Post")]:
            if pre in sub.index and post in sub.index:
                xi = order.index(pre); xf = order.index(post)
                ax.plot([xi, xf], [sub.loc[pre, ycol], sub.loc[post, ycol]],
                        color="#999999", lw=0.8, alpha=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
plot_panel(ax1, mac, "pct", "Macrophages (%)", "Macrophage Abundance")
plot_panel(ax2, cn2, "cn_pct", "CN2 Proportion (%)", "CN2 (Macrophage-enriched)")

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
