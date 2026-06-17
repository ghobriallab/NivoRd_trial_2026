#!/usr/bin/env python3
"""Supplementary Figure S5A-B: Ligand-receptor interaction analysis (exploratory).
(A) Plasma cell-macrophage axis  (B) CD4 T cell-macrophage axis

LIANA differential interactions (Cohen's d, DR vs P/NR). Reported as
exploratory in the manuscript Discussion.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

import sys; sys.path.insert(0, "."); from _00_common import FIG_DIR, DATA_LOCAL
INPUT_PM = str(DATA_LOCAL / "differential_plasma_macrophage_pretreatment.csv")
INPUT_CD4 = str(DATA_LOCAL / "differential_cd4_macrophage_pretreatment.csv")
OUTPUT = str(FIG_DIR / "SuppFigS5")

COLOR_DR, COLOR_PNR = "#2980b9", "#e74c3c"

HIGHLIGHT = ["LGALS9", "ANXA2", "TGFB1 - CXCR4"]

def clean_label(s):
    return (s.replace("Plasma Cell -> Macrophages: ", "PC \u2192 Mac: ")
             .replace("Macrophages -> Plasma Cell: ", "Mac \u2192 PC: ")
             .replace("CD4 T -> Macrophages: ", "CD4 \u2192 Mac: ")
             .replace("Macrophages -> CD4 T: ", "Mac \u2192 CD4: ")
             .replace("_", "/"))

def plot_panel(df, title, ax, n_top=15, xlim_cap=None):
    top = df.nlargest(n_top, "abs_d").sort_values("cohens_d")
    labels = [clean_label(x) for x in top.interaction]
    vals = top.cohens_d.values
    disp = np.clip(vals, -xlim_cap, xlim_cap) if xlim_cap else vals
    capped = np.abs(vals) > xlim_cap if xlim_cap else np.zeros(len(vals), dtype=bool)

    ax.barh(range(len(labels)), disp, color=[COLOR_DR if v > 0 else COLOR_PNR for v in vals],
            edgecolor="none", height=0.7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Cohen's d (DR - P/NR)", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for i, (v, c) in enumerate(zip(vals, capped)):
        if c:
            ax.text(disp[i] + (0.3 if v > 0 else -0.3), i, f"d={v:.1f}",
                    va="center", ha="left" if v > 0 else "right",
                    fontsize=9, fontstyle="italic", color="#555555")

    for i, interaction in enumerate(top.interaction):
        if any(k in interaction for k in HIGHLIGHT):
            ax.get_yticklabels()[i].set_fontweight("bold")

# Load
pm = pd.read_csv(INPUT_PM)
cd4m = pd.read_csv(INPUT_CD4)

# Combined figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
plot_panel(pm, "Plasma Cell - Macrophage Axis", ax1, xlim_cap=15)
plot_panel(cd4m, "CD4 T Cell - Macrophage Axis", ax2)
fig.legend(handles=[Patch(facecolor=COLOR_DR, label="Enriched in deep responders"),
                    Patch(facecolor=COLOR_PNR, label="Enriched in partial/non-responders")],
           loc="lower center", ncol=2, fontsize=11, frameon=False, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.04, 1, 1])
for fmt in ["pdf", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
