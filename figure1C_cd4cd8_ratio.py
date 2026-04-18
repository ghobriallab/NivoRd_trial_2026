#!/usr/bin/env python3
"""Figure 1C: CD4/CD8 ratio across response groups and timepoints."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_OBS = "phenotyped-obs.csv"
INPUT_CLIN = "nivo_clinical.csv"
OUTPUT = "Figure1C_cd4cd8_ratio"
COLOR_R, COLOR_NR = "#2980b9", "#e74c3c"

obs = pd.read_csv(INPUT_OBS)
clin = pd.read_csv(INPUT_CLIN)

if "group" not in obs.columns:
    obs = obs.merge(clin[["pt", "group"]].drop_duplicates(), on="pt", how="left")

obs["pt_tp"] = obs["pt"].astype(str) + "_" + obs["tp"].astype(str)

tcells = obs[obs.pheno.isin(["T4", "T8"])].copy()
counts = tcells.groupby(["pt_tp", "pheno"]).size().unstack(fill_value=0)
counts["ratio"] = counts["T4"] / counts["T8"].replace(0, np.nan)
counts = counts.dropna(subset=["ratio"]).reset_index()

meta = obs[["pt_tp", "pt", "tp", "group"]].drop_duplicates()
counts = counts.merge(meta, on="pt_tp")

counts["x"] = counts["group"] + "_" + counts["tp"]
order = ["good_Baseline", "good_Post", "poor_Baseline", "poor_Post"]
labels = ["R\nBaseline", "R\nPost", "NR\nBaseline", "NR\nPost"]
colors = [COLOR_R, COLOR_R, COLOR_NR, COLOR_NR]

fig, ax = plt.subplots(figsize=(4, 4.5))
for i, (grp, lab, col) in enumerate(zip(order, labels, colors)):
    vals = counts.loc[counts.x == grp, "ratio"]
    bp = ax.boxplot([vals.values], positions=[i], widths=0.5, patch_artist=True, showfliers=False)
    bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha(0.6)
    bp["medians"][0].set_color("black")
    ax.scatter(np.full(len(vals), i) + np.random.default_rng(42).uniform(-0.1, 0.1, len(vals)),
               vals, c=col, s=40, zorder=3, edgecolors="white", linewidths=0.5)

for pt in counts.pt.unique():
    sub = counts[counts.pt == pt].set_index("x")
    for pre, post in [("good_Baseline", "good_Post"), ("poor_Baseline", "poor_Post")]:
        if pre in sub.index and post in sub.index:
            xi = order.index(pre); xf = order.index(post)
            ax.plot([xi, xf], [sub.loc[pre, "ratio"], sub.loc[post, "ratio"]],
                    color="#999999", lw=0.8, alpha=0.5)

ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel("CD4/CD8 Ratio", fontsize=12, fontweight="bold")
ax.set_title("CD4/CD8 Ratio", fontsize=12, fontweight="bold")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
