#!/usr/bin/env python3
"""Supplementary Figure S3A: CN proportions (CN3, CN4, CN7) across groups."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_CN = "pheno_8clus-sample_pct.csv"
INPUT_CLIN = "nivo_clinical.csv"
OUTPUT = "FigureS3A_cn_proportions"
COLOR_R, COLOR_NR = "#2980b9", "#e74c3c"

cn_pct = pd.read_csv(INPUT_CN)
clin = pd.read_csv(INPUT_CLIN)

cn_long = cn_pct.melt(id_vars=[c for c in cn_pct.columns if not c.startswith("CN")],
                       var_name="cn", value_name="cn_pct") if "CN0" in cn_pct.columns else cn_pct
if "cn" not in cn_long.columns:
    cn_long = cn_pct.melt(id_vars=["pt_tp"], var_name="cn", value_name="cn_pct")

if "pt" not in cn_long.columns:
    cn_long["pt"] = cn_long["pt_tp"].str.extract(r"(Patient\d+|P\d+)", expand=False)
if "tp" not in cn_long.columns:
    cn_long["tp"] = cn_long["pt_tp"].apply(lambda x: "Baseline" if "Baseline" in x else "Post")
if "group" not in cn_long.columns:
    cn_long = cn_long.merge(clin[["pt", "group"]].drop_duplicates(), on="pt", how="left")

cn_agg = cn_long.groupby(["pt_tp", "cn", "pt", "tp", "group"])["cn_pct"].mean().reset_index()

order = ["good_Baseline", "good_Post", "poor_Baseline", "poor_Post"]
labels = ["R\nBL", "R\nPost", "NR\nBL", "NR\nPost"]
colors = [COLOR_R, COLOR_R, COLOR_NR, COLOR_NR]
cn_panels = ["CN7", "CN3", "CN4"]
cn_titles = ["CN7 (CD8-enriched)", "CN3 (PC-enriched)", "CN4 (PC-enriched)"]

fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
for ax, cn_label, cn_title in zip(axes, cn_panels, cn_titles):
    sub = cn_agg[cn_agg.cn == cn_label].copy()
    sub["x"] = sub["group"] + "_" + sub["tp"]
    for i, (grp, lab, col) in enumerate(zip(order, labels, colors)):
        vals = sub.loc[sub.x == grp, "cn_pct"]
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
                ax.plot([xi, xf], [s.loc[pre, "cn_pct"], s.loc[post, "cn_pct"]],
                        color="#999999", lw=0.8, alpha=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("CN Proportion (%)", fontsize=11, fontweight="bold")
    ax.set_title(cn_title, fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
