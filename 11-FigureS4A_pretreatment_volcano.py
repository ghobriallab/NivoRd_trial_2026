#!/usr/bin/env python3
"""Supplementary Figure S4A: Pre-treatment macrophage DGE volcano plot."""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from adjustText import adjust_text

INPUT = "pretreatment_macrophage_DGE_results.csv"
OUTPUT = "FigureS4A_macrophage_volcano_pretreatment"
FDR, LFC = 0.05, 0.5
COLOR_NR, COLOR_R, COLOR_NS = "#e74c3c", "#2980b9", "#d5d8dc"

df = pd.read_csv(INPUT).dropna(subset=["padj", "log2FoldChange"]).copy()
df["-log10padj"] = -np.log10(df["padj"].clip(lower=1e-300))
df["cat"] = "ns"
df.loc[(df.padj < FDR) & (df.log2FoldChange > LFC), "cat"] = "up_NR"
df.loc[(df.padj < FDR) & (df.log2FoldChange < -LFC), "cat"] = "up_R"
n_p, n_np = (df.cat == "up_NR").sum(), (df.cat == "up_R").sum()

gene_col = "Gene" if "Gene" in df.columns else "gene"

fig, ax = plt.subplots(figsize=(5.5, 5))
colors = {"up_NR": COLOR_NR, "up_R": COLOR_R, "ns": COLOR_NS}
for cat, z in [("ns", 1), ("up_R", 3), ("up_NR", 3)]:
    m = df.cat == cat
    if not m.any(): continue
    ax.scatter(df.loc[m, "log2FoldChange"], df.loc[m, "-log10padj"],
               c=colors[cat], s=6 if cat == "ns" else 20,
               alpha=0.25 if cat == "ns" else 0.8,
               edgecolors="none", rasterized=True, zorder=z)

ax.axhline(-np.log10(FDR), color="gray", ls="--", lw=0.7, alpha=0.5)
ax.axvline(LFC, color="gray", ls="--", lw=0.7, alpha=0.5)
ax.axvline(-LFC, color="gray", ls="--", lw=0.7, alpha=0.5)

sig = df[(df.padj < FDR) & (df.log2FoldChange.abs() > LFC)].nsmallest(15, "padj")
texts = [ax.text(r.log2FoldChange, getattr(r, "-log10padj"), getattr(r, gene_col),
                 fontsize=10, ha="center", va="bottom", fontstyle="italic", fontweight="bold")
         for _, r in sig.iterrows()]
if texts:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6),
                expand=(1.4, 1.6), force_text=(0.8, 1.0), force_points=(0.5, 0.5))

ax.set_xlabel("log\u2082FC (NR / R)", fontsize=13, fontweight="bold")
ax.set_ylabel("\u2212log\u2081\u2080(FDR)", fontsize=13, fontweight="bold")
ax.tick_params(labelsize=11)
ax.set_title(f"Pre-treatment Macrophages\n{n_p + n_np} DEGs (FDR<0.05, |log\u2082FC|>0.5)",
             fontsize=12, fontweight="bold", linespacing=1.3)
ax.legend(handles=[
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_NR, markersize=9,
           label=f"Up in Non-Responders ({n_p})"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_R, markersize=9,
           label=f"Up in Responders ({n_np})"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_NS, markersize=7, label="n.s."),
], loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=10, frameon=False)

plt.tight_layout(pad=0.5)
fig.subplots_adjust(bottom=0.15)
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
