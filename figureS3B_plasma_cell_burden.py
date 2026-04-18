#!/usr/bin/env python3
"""Supplementary Figure S3B: Plasma cell burden (PC% and M-spike)."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_CLIN = "nivo_clinical.csv"
INPUT_MSPIKE = "nivo_MSPike.csv"
OUTPUT = "FigureS3B_plasma_cell_burden"
COLOR_R, COLOR_NR = "#2980b9", "#e74c3c"

clin = pd.read_csv(INPUT_CLIN)
mspike = pd.read_csv(INPUT_MSPIKE)

order = ["good_Baseline", "good_Post", "poor_Baseline", "poor_Post"]
labels = ["R\nBL", "R\nPost", "NR\nBL", "NR\nPost"]
colors = [COLOR_R, COLOR_R, COLOR_NR, COLOR_NR]

def prepare(df, val_col):
    df = df.copy()
    df[val_col] = pd.to_numeric(df[val_col], errors="coerce")
    df = df.dropna(subset=[val_col])
    if "tp" not in df.columns and "timepoint" in df.columns:
        df["tp"] = df["timepoint"]
    df["x"] = df["group"] + "_" + df["tp"]
    return df

def plot_panel(ax, df, val_col, ylabel, title):
    for i, (grp, lab, col) in enumerate(zip(order, labels, colors)):
        vals = df.loc[df.x == grp, val_col]
        if len(vals) == 0: continue
        bp = ax.boxplot([vals.values], positions=[i], widths=0.5, patch_artist=True, showfliers=False)
        bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha(0.6)
        bp["medians"][0].set_color("black")
        ax.scatter(np.full(len(vals), i) + np.random.default_rng(42).uniform(-0.1, 0.1, len(vals)),
                   vals, c=col, s=40, zorder=3, edgecolors="white", linewidths=0.5)
    for pt in df.pt.unique():
        s = df[df.pt == pt].set_index("x")
        for pre, post in [("good_Baseline", "good_Post"), ("poor_Baseline", "poor_Post")]:
            if pre in s.index and post in s.index:
                xi = order.index(pre); xf = order.index(post)
                ax.plot([xi, xf], [s.loc[pre, val_col], s.loc[post, val_col]],
                        color="#999999", lw=0.8, alpha=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

pc_col = "PC_pct" if "PC_pct" in clin.columns else "pc_pct"
ms_col = "MSpike" if "MSpike" in mspike.columns else "mspike"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
plot_panel(ax1, prepare(clin, pc_col), pc_col, "CD138+ % (trephine)", "Plasma Cell %")
plot_panel(ax2, prepare(mspike, ms_col), ms_col, "M-spike (g/dL)", "M-spike")

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
