#!/usr/bin/env python3
"""Figure 2C: Post-treatment macrophage GSEA hallmark pathways."""

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

INPUT = "master_gsea_results.csv"
OUTPUT = "Figure2C_macrophage_gsea_hallmarks"
COLOR_NR, COLOR_R = "#e74c3c", "#2980b9"
FDR_DISPLAY, TOP_N = 0.25, 10

# Load and filter
gsea = pd.read_csv(INPUT)
mac = gsea[(gsea.CellType.str.contains("Macrophages", case=False)) &
           (gsea.Timepoint.str.contains("Post", case=False)) &
           (gsea.Library == "Hallmarks")].copy()

sig = mac[mac["FDR q-val"] <= FDR_DISPLAY]
if len(sig) < TOP_N:
    sig = mac.nsmallest(TOP_N, "NOM p-val")
sig = sig.sort_values("NES").copy()
sig["label"] = sig.Term.str.replace("_", " ").str.title().str[:40]

# Plot
fig, ax = plt.subplots(figsize=(6, 4.5))
colors = [COLOR_R if n < 0 else COLOR_NR for n in sig.NES]
ax.barh(range(len(sig)), sig.NES.values, color=colors, edgecolor="none", height=0.7)
ax.set_yticks(range(len(sig)))
ax.set_yticklabels(sig.label.values, fontsize=11)
ax.set_xlabel("Normalized Enrichment Score", fontsize=12, fontweight="bold")
ax.tick_params(axis="x", labelsize=11)
ax.axvline(0, color="black", linewidth=0.8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for i, (nes, fdr) in enumerate(zip(sig.NES, sig["FDR q-val"])):
    q = f"q={fdr:.1e}" if fdr < 0.001 else f"q={fdr:.3f}" if fdr < 0.01 else f"q={fdr:.2f}"
    if abs(nes) > 0.8:
        ax.text(nes - 0.05 if nes > 0 else nes + 0.05, i, q, va="center",
                ha="right" if nes > 0 else "left", fontsize=9, color="white", fontweight="bold")
    else:
        ax.text(nes + 0.05 if nes > 0 else nes - 0.05, i, q, va="center",
                ha="left" if nes > 0 else "right", fontsize=9, color="#555555")

ax.set_title("Post-treatment Macrophage GSEA (Hallmark)", fontsize=12, fontweight="bold", pad=30)
ax.legend(handles=[Patch(facecolor=COLOR_R, label="Enriched in R"),
                   Patch(facecolor=COLOR_NR, label="Enriched in NR")],
          loc="lower center", bbox_to_anchor=(0.5, 1.005), fontsize=10, frameon=False,
          ncol=2, handletextpad=0.4, columnspacing=1.5)

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
