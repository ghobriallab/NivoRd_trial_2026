#!/usr/bin/env python3
"""Figure 2D: CD4/CD8 ratio vs CN2 and CN4 proportions (two subpanels)."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

INPUT_OBS = "phenotyped-obs.csv"
INPUT_CN = "pheno_8clus-sample_pct.csv"
INPUT_CLIN = "nivo_clinical.csv"
OUTPUT = "Figure2D_cd4cd8_vs_cn"
COLOR_R, COLOR_NR = "#2980b9", "#e74c3c"

obs = pd.read_csv(INPUT_OBS)
clin = pd.read_csv(INPUT_CLIN)
cn_pct = pd.read_csv(INPUT_CN)

if "group" not in obs.columns:
    obs = obs.merge(clin[["pt", "group"]].drop_duplicates(), on="pt", how="left")
obs["pt_tp"] = obs["pt"].astype(str) + "_" + obs["tp"].astype(str)

# CD4/CD8 ratio
tcells = obs[obs.pheno.isin(["T4", "T8"])]
tc = tcells.groupby(["pt_tp", "pheno"]).size().unstack(fill_value=0)
tc["ratio"] = tc["T4"] / tc["T8"].replace(0, np.nan)
tc = tc.dropna(subset=["ratio"]).reset_index()

# CN proportions
cn_long = cn_pct.melt(id_vars=[c for c in cn_pct.columns if not c.startswith("CN")],
                       var_name="cn", value_name="cn_pct") if "CN0" in cn_pct.columns else cn_pct
if "cn" not in cn_long.columns:
    cn_long = cn_pct.melt(id_vars=["pt_tp"], var_name="cn", value_name="cn_pct")
cn_agg = cn_long.groupby(["pt_tp", "cn"])["cn_pct"].mean().reset_index()

meta = obs[["pt_tp", "pt", "tp", "group"]].drop_duplicates()

def plot_corr(ax, cn_label, title):
    cn_sub = cn_agg[cn_agg.cn == cn_label].rename(columns={"cn_pct": "y"})
    merged = tc[["pt_tp", "ratio"]].merge(cn_sub[["pt_tp", "y"]], on="pt_tp")
    merged = merged.merge(meta, on="pt_tp")

    for grp, col, lab in [("good", COLOR_R, "R"), ("poor", COLOR_NR, "NR")]:
        m = merged.group == grp
        ax.scatter(merged.loc[m, "ratio"], merged.loc[m, "y"], c=col, s=50,
                   edgecolors="white", linewidths=0.5, label=lab, zorder=3)

    r, p = pearsonr(merged["ratio"], merged["y"])
    z = np.polyfit(merged["ratio"], merged["y"], 1)
    xr = np.linspace(merged["ratio"].min(), merged["ratio"].max(), 100)
    ax.plot(xr, np.polyval(z, xr), color="#555555", lw=1.5, ls="--", alpha=0.7)
    ax.text(0.05, 0.95, f"R={r:.2f}, p={p:.3f}", transform=ax.transAxes,
            fontsize=10, va="top", fontweight="bold")

    ax.set_xlabel("CD4/CD8 Ratio", fontsize=11, fontweight="bold")
    ax.set_ylabel(f"{cn_label} (%)", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(fontsize=10, frameon=False)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
plot_corr(ax1, "CN2", "CN2 (Macrophage-enriched)")
plot_corr(ax2, "CN4", "CN4 (Plasma cell-enriched)")

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
