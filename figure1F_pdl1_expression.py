#!/usr/bin/env python3
"""Figure 1F: PD-L1 expression by compartment (plasma cells + macrophages)."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_OBS = "phenotyped-obs.csv"
INPUT_MARKERS = "combat_all_markers.csv"
INPUT_CLIN = "nivo_clinical.csv"
OUTPUT = "Figure1F_pdl1_expression"
COLOR_R, COLOR_NR = "#2980b9", "#e74c3c"

obs = pd.read_csv(INPUT_OBS)
markers = pd.read_csv(INPUT_MARKERS)
clin = pd.read_csv(INPUT_CLIN)

if "group" not in obs.columns:
    obs = obs.merge(clin[["pt", "group"]].drop_duplicates(), on="pt", how="left")
obs["pt_tp"] = obs["pt"].astype(str) + "_" + obs["tp"].astype(str)

merge_key = obs.columns.intersection(markers.columns).difference(["PDL1"]).tolist()
if not merge_key:
    obs = pd.concat([obs.reset_index(drop=True), markers[["PDL1"]].reset_index(drop=True)], axis=1)
else:
    obs = obs.merge(markers[merge_key + ["PDL1"]], on=merge_key, how="left")

order = ["good_Baseline", "good_Post", "poor_Baseline", "poor_Post"]
labels = ["R\nBL", "R\nPost", "NR\nBL", "NR\nPost"]
colors = [COLOR_R, COLOR_R, COLOR_NR, COLOR_NR]

def plot_pdl1(ax, pheno_label, pheno_code, title):
    sub = obs[obs.pheno == pheno_code].copy()
    threshold = sub.PDL1.quantile(0.9)
    agg = sub.groupby("pt_tp").apply(lambda g: (g.PDL1 >= threshold).mean() * 100).reset_index(name="pct_hi")
    meta = obs[["pt_tp", "pt", "tp", "group"]].drop_duplicates()
    agg = agg.merge(meta, on="pt_tp")
    agg["x"] = agg["group"] + "_" + agg["tp"]

    for i, (grp, lab, col) in enumerate(zip(order, labels, colors)):
        vals = agg.loc[agg.x == grp, "pct_hi"]
        bp = ax.boxplot([vals.values], positions=[i], widths=0.5, patch_artist=True, showfliers=False)
        bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha(0.6)
        bp["medians"][0].set_color("black")
        ax.scatter(np.full(len(vals), i) + np.random.default_rng(42).uniform(-0.1, 0.1, len(vals)),
                   vals, c=col, s=40, zorder=3, edgecolors="white", linewidths=0.5)
    for pt in agg.pt.unique():
        s = agg[agg.pt == pt].set_index("x")
        for pre, post in [("good_Baseline", "good_Post"), ("poor_Baseline", "poor_Post")]:
            if pre in s.index and post in s.index:
                xi = order.index(pre); xf = order.index(post)
                ax.plot([xi, xf], [s.loc[pre, "pct_hi"], s.loc[post, "pct_hi"]],
                        color="#999999", lw=0.8, alpha=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("% PD-L1 high", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
plot_pdl1(ax1, "Plasma Cells", "PC", "PD-L1: Plasma Cells")
plot_pdl1(ax2, "Macrophages", "Mac", "PD-L1: Macrophages")

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
