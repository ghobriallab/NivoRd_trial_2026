#!/usr/bin/env python3
"""Supplementary Figure S2A: Immune cell type marker expression heatmap."""

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import issparse

import sys; sys.path.insert(0, "."); from _00_common import FIG_DIR, H5AD_IMMUNE
INPUT = str(H5AD_IMMUNE)
OUTPUT = str(FIG_DIR / "SuppFigS2A")
CELLTYPE_COL = "celltype_refined"

RENAME = {"B": "B Cells", "DC": "Dendritic Cells"}

MARKERS = {
    "CD4 Naive": ["IL7R", "LEF1", "CCR7", "SELL"],
    "CD4 Effector": ["TNF", "IFNG", "ISG15"],
    "Tregs": ["FOXP3", "IL2RA", "CTLA4", "TIGIT"],
    "CD8 Naive": ["CD8A", "CD8B", "TCF7"],
    "CD8 GZMK+ Effector": ["GZMK", "CCL5", "NKG7"],
    "CD8 GZMB+ Effector": ["GZMB", "PRF1", "GNLY", "FGFBP2"],
    "NK/NKT": ["NCAM1", "TYROBP", "KLRB1", "KLRD1"],
    "Classical Monocytes": ["CD14", "S100A8", "S100A9", "VCAN"],
    "Non-classical Monocytes": ["FCGR3A", "MS4A7", "CDKN1C", "CSF1R"],
    "Macrophages": ["CD68", "C1QA", "C1QB"],
    "B Cells": ["CD79A", "MS4A1", "CD79B", "BANK1"],
    "Dendritic Cells": ["FCER1A", "CD1C", "CLEC10A", "HLA-DRA"],
}

# Load
adata = ad.read_h5ad(INPUT)
obs = adata.obs.copy()
obs[CELLTYPE_COL] = obs[CELLTYPE_COL].replace(RENAME)
var_names = list(adata.var_names)

celltypes = list(MARKERS.keys())
all_markers = [g for ct in celltypes for g in MARKERS[ct]]
unique = list(dict.fromkeys(all_markers))
available = [g for g in unique if g in var_names]

idx = [var_names.index(g) for g in available]
X = adata.X[:, idx].toarray() if issparse(adata.X) else adata.X[:, idx]

mean_expr = pd.DataFrame(index=available, columns=celltypes, dtype=float)
for ct in celltypes:
    m = (obs[CELLTYPE_COL] == ct).values
    for i, gene in enumerate(available):
        mean_expr.at[gene, ct] = float(np.mean(X[m, i]))
del adata, X

zscore = mean_expr.astype(float).apply(
    lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0, axis=1)
display_genes = [g for g in all_markers if g in available]

# Plot
fig, ax = plt.subplots(figsize=(8, 12))
im = ax.imshow(zscore.loc[display_genes, celltypes].values, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
ax.set_xticks(range(len(celltypes)))
ax.set_xticklabels(celltypes, fontsize=10, fontweight="bold", rotation=45, ha="right")
ax.set_yticks(range(len(display_genes)))
ax.set_yticklabels(display_genes, fontsize=10, fontweight="bold")

gi = 0
for ct in celltypes:
    n = len([g for g in MARKERS[ct] if g in available])
    if gi > 0: ax.axhline(y=gi - 0.5, color="black", linewidth=0.8)
    gi += n

ax.set_title("Cell Type Marker Expression (Z-score)", fontsize=14, fontweight="bold", pad=10)
cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.04)
cbar.set_label("Z-score", fontsize=11, fontweight="bold")

plt.tight_layout()
for fmt in ["pdf", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
