#!/usr/bin/env python3
"""Supplementary Figure S2E: CD138+ bone marrow plasma cell UMAP by patient."""

import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt

import sys; sys.path.insert(0, "."); from _00_common import FIG_DIR, H5AD_TUMOR
INPUT = str(H5AD_TUMOR)
OUTPUT = str(FIG_DIR / "SuppFigS2E")
PATIENT_COL = "patient"

RESPONSE_MAP = {
    "Patient1": ("Case #1", "P"), "Patient4": ("Case #4", "P"),
    "Patient5": ("Case #5", "NP"), "Patient6": ("Case #6", "NP"),
    "Patient7": ("Case #7", "NP"),
}

COLORS = {
    "Case #1": "#e41a1c", "Case #4": "#ff7f00", "Case #5": "#4daf4a",
    "Case #6": "#377eb8", "Case #7": "#984ea3",
}

CASE_ORDER = ["Case #1", "Case #4", "Case #5", "Case #6", "Case #7"]

# Load and filter to treatment-arm tumor cells
adata = ad.read_h5ad(INPUT)
if "Source" in adata.obs.columns:
    adata = adata[adata.obs.Source == "Nivo"].copy()
if "Malignant_PC_Immune" in adata.obs.columns:
    adata = adata[adata.obs.Malignant_PC_Immune == "Malignant_PC"].copy()

adata.obs["Case"] = adata.obs[PATIENT_COL].map(
    {k: v[0] for k, v in RESPONSE_MAP.items() if k in str(adata.obs[PATIENT_COL].values)})
# Substring matching for pre/post sample IDs
for idx in adata.obs.index:
    pid = str(adata.obs.loc[idx, PATIENT_COL])
    for key, (case, _) in RESPONSE_MAP.items():
        if key in pid:
            adata.obs.loc[idx, "Case"] = case

adata = adata[adata.obs.Case.isin(CASE_ORDER)].copy()
if "X_umap" not in adata.obsm:
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)

# Plot
fig, ax = plt.subplots(figsize=(5.5, 5))
sc.pl.umap(adata, color="Case", palette=COLORS, ax=ax, show=False,
           frameon=True, title="", size=30, legend_loc="none")

n_total = adata.shape[0]
for case in CASE_ORDER:
    m = adata.obs.Case == case
    if not m.any(): continue
    coords = adata.obsm["X_umap"][m]
    ax.annotate(f"{case}\n(n={m.sum()})", xy=(coords[:, 0].mean(), coords[:, 1].mean()),
                fontsize=10, fontweight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor=COLORS[case], alpha=0.9, linewidth=2))

ax.set_title(f"CD138$^+$ Bone Marrow Plasma Cells (n={n_total:,})", fontsize=12, fontweight="bold")
ax.set_xlabel("UMAP 1", fontsize=13, fontweight="bold")
ax.set_ylabel("UMAP 2", fontsize=13, fontweight="bold")

plt.tight_layout(pad=0.5)
for fmt in ["pdf", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
