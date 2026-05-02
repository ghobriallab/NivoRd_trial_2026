"""Shared helper for single-panel volcano plots (Fig 2B, Fig 2C)."""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from adjustText import adjust_text

import _00_common as cm

warnings.filterwarnings("ignore")
mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

FDR, LFC = 0.05, 0.5
COLOR_NR, COLOR_R, COLOR_NS = "#e74c3c", "#2980b9", "#d5d8dc"


def make_volcano(csv_path, title, save_tag):
    df = pd.read_csv(csv_path).rename(columns={"Gene": "gene"})
    df = df.dropna(subset=["padj", "log2FoldChange"]).copy()
    # Drop LMM fits that hit boundary / singular RE covariance: their
    # SE/z/p are unreliable. Effect estimates can still feed GSEA, but
    # significance counts and labels are restricted to converged fits.
    if "converged" in df.columns:
        df = df[df["converged"].astype(bool)].copy()
    df["-log10padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["cat"] = "ns"
    df.loc[(df.padj < FDR) & (df.log2FoldChange > LFC), "cat"] = "up_NR"
    df.loc[(df.padj < FDR) & (df.log2FoldChange < -LFC), "cat"] = "up_R"
    n_NR = int((df.cat == "up_NR").sum())
    n_R  = int((df.cat == "up_R").sum())

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    for cat, z in [("ns", 1), ("up_R", 3), ("up_NR", 3)]:
        m = df.cat == cat
        if not m.any():
            continue
        ax.scatter(df.loc[m, "log2FoldChange"], df.loc[m, "-log10padj"],
                    c={"up_NR": COLOR_NR, "up_R": COLOR_R, "ns": COLOR_NS}[cat],
                    s=6 if cat == "ns" else 18,
                    alpha=0.25 if cat == "ns" else 0.85,
                    edgecolors="none", rasterized=True, zorder=z)
    ax.axhline(-np.log10(FDR), color="gray", ls="--", lw=0.7, alpha=0.5)
    ax.axvline(LFC, color="gray", ls="--", lw=0.7, alpha=0.5)
    ax.axvline(-LFC, color="gray", ls="--", lw=0.7, alpha=0.5)

    sig = df[(df.padj < FDR) & (df.log2FoldChange.abs() > LFC)].nsmallest(15, "padj")
    texts = [ax.text(r.log2FoldChange, getattr(r, "-log10padj"), r.gene,
                      fontsize=8, ha="center", va="bottom",
                      fontstyle="italic", fontweight="bold")
             for _, r in sig.iterrows()]
    if texts:
        adjust_text(texts, ax=ax,
                      arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6),
                      expand=(1.4, 1.6))

    ax.set_xlabel("log₂FC (NR / R)", fontsize=11, fontweight="bold")
    ax.set_ylabel("−log₁₀(FDR)", fontsize=11, fontweight="bold")
    ax.tick_params(labelsize=10)
    ax.set_title(f"{title}\n{n_R + n_NR} DEGs ({n_R} up R, {n_NR} up NR)",
                  fontsize=10.5, fontweight="bold", linespacing=1.3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_NR,
                markersize=8, label="Up in NR"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_R,
                markersize=8, label="Up in R"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_NS,
                markersize=6, label="n.s."),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3,
       fontsize=9, frameon=False)

    plt.tight_layout(pad=0.5)
    fig.subplots_adjust(bottom=0.18)
    cm.save_fig(fig, save_tag)
    plt.close(fig)
    return n_R, n_NR
