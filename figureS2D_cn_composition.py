#!/usr/bin/env python3
"""Supplementary Figure S2D: Cellular neighborhood composition stacked bar."""

import matplotlib.pyplot as plt
import pandas as pd

INPUT_CN = "pheno_8clus-cells_data.csv"
OUTPUT = "FigureS2D_cn_composition"

PALETTE = {
    "T4": "#e41a1c", "T8": "#6a3d9a", "B": "#888888", "Mac": "#1f4e9f",
    "Myeloid": "#33a02c", "MK": "#ff7f00", "PC": "#f781bf",
}

cells = pd.read_csv(INPUT_CN)

comp = cells.groupby(["cn_labels", "pheno"]).size().unstack(fill_value=0)
comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100

cn_order = sorted(comp_pct.index, key=lambda x: int(str(x).replace("CN", "")) if "CN" in str(x) else x)
comp_pct = comp_pct.reindex(cn_order)

fig, ax = plt.subplots(figsize=(8, 5))
bottom = pd.Series(0.0, index=comp_pct.index)
phenotypes = [p for p in PALETTE if p in comp_pct.columns] + \
             [p for p in comp_pct.columns if p not in PALETTE]

for pheno in phenotypes:
    if pheno not in comp_pct.columns:
        continue
    color = PALETTE.get(pheno, "#cccccc")
    ax.bar(range(len(comp_pct)), comp_pct[pheno], bottom=bottom,
           color=color, edgecolor="white", linewidth=0.5, label=pheno)
    bottom += comp_pct[pheno]

ax.set_xticks(range(len(comp_pct)))
ax.set_xticklabels(comp_pct.index, fontsize=11, fontweight="bold")
ax.set_ylabel("Proportion (%)", fontsize=12, fontweight="bold")
ax.set_xlabel("Cellular Neighborhood", fontsize=12, fontweight="bold")
ax.set_title("CN Composition", fontsize=13, fontweight="bold")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=10, frameon=False)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
