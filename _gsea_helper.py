"""Shared helper for single-panel GSEA bar plots (Fig 2D, Fig 2E)."""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
import gseapy as gp

import _common as cm

warnings.filterwarnings("ignore")
mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

GMT = str(cm.GMT_PATH)
COLOR_NR, COLOR_R = "#e74c3c", "#2980b9"


def run_gsea(csv_path):
    df = pd.read_csv(csv_path)
    rnk = df.dropna(subset=["coef", "Gene"])[["Gene", "coef"]]
    rnk.columns = ["Name", "Score"]
    rnk = (rnk.sort_values("Score", ascending=False)
                .drop_duplicates("Name", keep="first"))
    pre = gp.prerank(rnk=rnk, gene_sets=GMT, threads=4, seed=0,
                      permutation_num=5000, min_size=15, max_size=500,
                      no_plot=True, verbose=False, outdir=None)
    return pre.res2d.copy()


def make_gsea_panel(csv_path, save_csv_path, title, save_tag, n_top=10):
    res = run_gsea(csv_path)
    res.to_csv(save_csv_path, index=False)

    sig = res[res["FDR q-val"] <= 0.25]
    if len(sig) < n_top:
        sig = res.nsmallest(n_top, "NOM p-val")
    sig = sig.sort_values("NES").copy()
    sig["label"] = sig["Term"].str.replace("_", " ").str.title().str[:36]

    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    colors = [COLOR_R if n < 0 else COLOR_NR for n in sig["NES"]]
    ax.barh(range(len(sig)), sig["NES"].values, color=colors,
             edgecolor="none", height=0.7)
    ax.set_yticks(range(len(sig)))
    ax.set_yticklabels(sig["label"].values, fontsize=10)
    ax.set_xlabel("Normalized Enrichment Score", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", labelsize=10)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    for i, (nes, fdr) in enumerate(zip(sig["NES"], sig["FDR q-val"])):
        q = (f"q={fdr:.1e}" if fdr < 0.001 else
              f"q={fdr:.3f}" if fdr < 0.01 else f"q={fdr:.2f}")
        if abs(nes) > 0.8:
            ax.text(nes - 0.05 if nes > 0 else nes + 0.05, i, q,
                      va="center", ha="right" if nes > 0 else "left",
                      fontsize=9, color="white", fontweight="bold")
        else:
            ax.text(nes + 0.05 if nes > 0 else nes - 0.05, i, q,
                      va="center", ha="left" if nes > 0 else "right",
                      fontsize=9, color="#555555")
    ax.set_title(title, fontsize=11, fontweight="bold", pad=18)
    ax.legend(handles=[Patch(facecolor=COLOR_R, label="Enriched in R"),
                        Patch(facecolor=COLOR_NR, label="Enriched in NR")],
               loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=9,
               frameon=False, ncol=2, handletextpad=0.4, columnspacing=1.5)
    plt.tight_layout()
    cm.save_fig(fig, save_tag)
    plt.close(fig)
