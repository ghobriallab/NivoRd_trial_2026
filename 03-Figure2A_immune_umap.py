#!/usr/bin/env python3
"""Figure 2A: CD138- bone marrow immune cell UMAP by cell type."""

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from adjustText import adjust_text

INPUT = "CD138neg_immune_cells.h5ad"
OUTPUT = "Figure2A_CD138neg_UMAP_by_celltype"
CELLTYPE_COL = "celltype_refined"

EXCLUDE = ["Doublets", "Proliferating", "Hemoglobin", "ILC-like"]
RENAME = {"B": "B Cells", "DC": "Dendritic Cells"}

PALETTE = {
    "CD4 Naive": "#e41a1c", "CD4 Effector": "#ff7f00", "Tregs": "#b8860b",
    "CD8 Naive": "#6a3d9a", "CD8 GZMK+ Effector": "#e31a8f",
    "CD8 GZMB+ Effector": "#f781bf", "NK/NKT": "#1b9e77",
    "Classical Monocytes": "#33a02c", "Non-classical Monocytes": "#00bfc4",
    "Macrophages": "#1f4e9f", "B Cells": "#888888", "Dendritic Cells": "#d95f02",
}

# Load and prepare
adata = ad.read_h5ad(INPUT)
obs = adata.obs.copy()
umap = adata.obsm["X_umap"].copy()
del adata

mask = ~obs[CELLTYPE_COL].isin(EXCLUDE)
obs = obs[mask].copy()
umap = umap[mask.values]
obs[CELLTYPE_COL] = obs[CELLTYPE_COL].replace(RENAME)

counts = obs[CELLTYPE_COL].value_counts()
celltypes = counts[counts > 0].index.tolist()
n_total = counts.sum()

# Density centroid for spread-out clusters
def density_centroid(coords, percentile=70):
    if len(coords) < 50:
        return coords.mean(axis=0)
    idx = np.random.default_rng(42).choice(len(coords), min(2000, len(coords)), replace=False)
    kde = gaussian_kde(coords[idx].T)
    d = kde(coords.T)
    return coords[d >= np.percentile(d, percentile)].mean(axis=0)

# Plot
fig, ax = plt.subplots(figsize=(5.5, 5))
for ct in reversed(celltypes):
    m = obs[CELLTYPE_COL] == ct
    ax.scatter(umap[m.values, 0], umap[m.values, 1], c=PALETTE.get(ct, "#888"),
               s=0.8, alpha=0.7, rasterized=True, linewidths=0)

ax.set_title(f"CD138$^-$ Bone Marrow Immune Cells (n={n_total:,})", fontsize=12, fontweight="bold")
ax.set_xlabel("UMAP 1", fontsize=13, fontweight="bold")
ax.set_ylabel("UMAP 2", fontsize=13, fontweight="bold")
ax.set_xticks([]); ax.set_yticks([])

texts = []
for ct in celltypes:
    m = obs[CELLTYPE_COL] == ct
    ct_umap = umap[m.values]
    cx, cy = density_centroid(ct_umap) if max(ct_umap.std(0)) > 3.0 else ct_umap.mean(0)
    t = ax.text(cx, cy, ct, fontsize=11, fontweight="bold", ha="center", va="center",
                color="black", bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                edgecolor=PALETTE.get(ct, "#888"), alpha=0.95, linewidth=2))
    texts.append(t)

adjust_text(texts, ax=ax, expand=(1.05, 1.05), force_text=(0.3, 0.3),
            force_points=(0.8, 0.8), ensure_inside_axes=True,
            arrowprops=dict(arrowstyle="-", lw=0, alpha=0))

plt.tight_layout(pad=0.5)
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
