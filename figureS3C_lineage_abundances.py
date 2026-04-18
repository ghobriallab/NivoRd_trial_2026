#!/usr/bin/env python3
"""Supplementary Figure S3C: Lineage abundances (B cells, CD11b+ myeloid, megakaryocytes)."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_OBS = "phenotyped-obs.csv"
INPUT_CLIN = "nivo_clinical.csv"
OUTPUT = "FigureS3C_lineage_abundances"
COLOR_R, COLOR_NR = "#2980b9", "#e74c3c"

obs = pd.read_csv(INPUT_OBS)
clin = pd.read_csv(INPUT_CLIN)

if "group" not in obs.columns:
    obs = obs.merge(clin[["pt", "group"]].drop_duplicates(), on="pt", how="left")
obs["pt_tp"] = obs["pt"].astype(str) + "_" + obs["tp"].astype(str)

counts = obs.groupby(["pt_tp", "pheno"]).size().reset_index(name="n")
totals = obs.groupby("pt_tp").size().reset_index(name="total")
counts = counts.merge(totals, on="pt_tp")
counts["pct"] = counts["n"] / counts["total"] * 100
meta = obs[["pt_tp", "pt", "tp", "group"]].drop_duplicates()
counts = counts.merge(meta, on="pt_tp")

order = ["good_Baseline", "good_Post", "poor_Baseline", "poor_Post"]
labels = ["R\nBL", "R\nPost", "NR\nBL", "NR\nPost"]
colors = [COLOR_R, COLOR_R, COLOR_NR, COLOR_NR]
phenotypes = ["B", "Myeloid", "MK"]
titles = ["B Cells", "CD11b+ Myeloid", "Megakaryocytes"]

fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
for ax, pheno, title in zip(axes, phenotypes, titles):
    sub = counts[counts.pheno == pheno].copy()
    sub["x"] = sub["group"] + "_" + sub["tp"]
    for i, (grp, lab, col) in enumerate(zip(order, labels, colors)):
        vals = sub.loc[sub.x == grp, "pct"]
        bp = ax.boxplot([vals.values], positions=[i], widths=0.5, patch_artist=True, showfliers=False)
        bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha(0.6)
        bp["medians"][0].set_color("black")
        ax.scatter(np.full(len(vals), i) + np.random.default_rng(42).uniform(-0.1, 0.1, len(vals)),
                   vals, c=col, s=40, zorder=3, edgecolors="white", linewidths=0.5)
    for pt in sub.pt.unique():
        s = sub[sub.pt == pt].set_index("x")
        for pre, post in [("good_Baseline", "good_Post"), ("poor_Baseline", "poor_Post")]:
            if pre in s.index and post in s.index:
                xi = order.index(pre); xf = order.index(post)
                ax.plot([xi, xf], [s.loc[pre, "pct"], s.loc[post, "pct"]],
                        color="#999999", lw=0.8, alpha=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Cells (%)", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
