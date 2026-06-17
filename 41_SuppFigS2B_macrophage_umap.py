#!/usr/bin/env python3
"""Supplementary Figure S2B: CD138- immune UMAP with macrophages highlighted.

Purpose:      Render Supp Fig S2B: UMAP of CD138- cells with macrophages highlighted.
Inputs:       CD138neg_immune_cells.h5ad.
Outputs:      figures/SuppFigS2B.pdf / .svg.
Dependencies: scanpy, matplotlib, _00_common, _01_data_loaders.
"""

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

import sys; sys.path.insert(0, "."); from _00_common import FIG_DIR, H5AD_IMMUNE
INPUT = str(H5AD_IMMUNE)
OUTPUT = str(FIG_DIR / "SuppFigS2B")
CELLTYPE_COL = "celltype_refined"

EXCLUDE = ["Doublets", "Proliferating", "Hemoglobin", "ILC-like"]
MAC_LABEL = "Macrophages"
COLOR_MAC = "#1f4e9f"
COLOR_BG = "#d5d8dc"

adata = ad.read_h5ad(INPUT)
obs = adata.obs.copy()
umap = adata.obsm["X_umap"].copy()
del adata

mask = ~obs[CELLTYPE_COL].isin(EXCLUDE)
obs = obs[mask].copy()
umap = umap[mask.values]
obs[CELLTYPE_COL] = obs[CELLTYPE_COL].replace({"B": "B Cells", "DC": "Dendritic Cells"})

is_mac = (obs[CELLTYPE_COL] == MAC_LABEL).values
n_mac = is_mac.sum()
n_total = len(obs)

fig, ax = plt.subplots(figsize=(5.5, 5))
ax.scatter(umap[~is_mac, 0], umap[~is_mac, 1], c=COLOR_BG,
           s=0.5, alpha=0.3, rasterized=True, linewidths=0, zorder=1)
ax.scatter(umap[is_mac, 0], umap[is_mac, 1], c=COLOR_MAC,
           s=2.0, alpha=0.8, rasterized=True, linewidths=0, zorder=2)

cx, cy = umap[is_mac].mean(axis=0)
ax.text(cx, cy, f"Macrophages\n(n={n_mac:,})", fontsize=11, fontweight="bold",
        ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                  edgecolor=COLOR_MAC, alpha=0.95, linewidth=2))

ax.set_title(f"CD138$^-$ Immune Cells (n={n_total:,})", fontsize=12, fontweight="bold")
ax.set_xlabel("UMAP 1", fontsize=13, fontweight="bold")
ax.set_ylabel("UMAP 2", fontsize=13, fontweight="bold")
ax.set_xticks([]); ax.set_yticks([])

plt.tight_layout(pad=0.5)
for fmt in ["pdf", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
