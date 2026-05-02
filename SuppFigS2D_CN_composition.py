"""Supp Fig S2D — CN composition stacked bar.

Regenerated from pheno_8clus-cn_centers.csv so the labeling matches the
main text and Supp S3 panels (which describe CN4 = 77% PC,
CN7 = 71% Mac, CN2 = PC-Mac contact zone).

The published version of this panel was built from an earlier k-means
run that produced different CN numbering (in that run CN2 was Mac-rich
and CN3 was PC-rich). Manuscript v4 uses the rerun CN labeling
throughout, so this panel is regenerated to match.

Output: figures/SuppFigS2D.{pdf,svg}
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import _common as cm

mpl.rcParams.update({"font.family": "Arial", "font.size": 10,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

CN_CENTERS = (cm.ROOT.parent / "_archive" / "pipeline" / "data"
              / "canonical_pipeline" / "pheno_8clus-cn_centers.csv")

# Stack order (bottom -> top) and colors matching the published S2D:
# B = light blue, Mac = dark blue, MK = light green, Myeloid = dark green,
# PC = pink, T4 = red, T8 = orange. Stack T8 at bottom, B at top.
STACK_ORDER = ["T8", "T4", "PC", "Myeloid", "MK", "Mac", "B"]
COLORS = {
    "B":       "#A6CEE3",  # light blue
    "Mac":     "#1F78B4",  # dark blue
    "MK":      "#B2DF8A",  # light green
    "Myeloid": "#33A02C",  # dark green
    "PC":      "#FB9A99",  # pink
    "T4":      "#E31A1C",  # red
    "T8":      "#FDBF6F",  # orange
}

df = pd.read_csv(CN_CENTERS, index_col=0)[STACK_ORDER]
# Cluster centers in the source CSV sum to 1.0 for 7 of 8 CNs; CN2 sums
# to 0.926 due to a numerical artifact in the upstream pipeline. Renormalize
# row-wise to 100% so each bar fills the axis cleanly. Compositional
# rankings are unchanged.
df = df.div(df.sum(axis=1), axis=0) * 100
df.index.name = "CN"

fig, ax = plt.subplots(figsize=(5.4, 3.6))
bottom = np.zeros(len(df))
xs = np.arange(len(df))
for ct in STACK_ORDER:
    vals = df[ct].values
    ax.bar(xs, vals, bottom=bottom, color=COLORS[ct], width=0.78,
           edgecolor="white", linewidth=0.5, label=ct)
    bottom += vals

ax.set_xticks(xs)
ax.set_xticklabels(df.index, fontsize=10)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylim(0, 101)
ax.set_ylabel("% cells", fontsize=11, fontweight="bold")
ax.tick_params(labelsize=10)
ax.spines[["top", "right"]].set_visible(False)

# Legend in the same order the publisher used (B, Mac, MK, Myeloid, PC, T4, T8)
legend_order = ["B", "Mac", "MK", "Myeloid", "PC", "T4", "T8"]
handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS[c]) for c in legend_order]
ax.legend(handles, legend_order, loc="center left",
          bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)

fig.tight_layout()
cm.save_fig(fig, "SuppFigS2D")
plt.close(fig)
print(f"wrote SuppFigS2D ({len(df)} CNs from {CN_CENTERS.name})")
