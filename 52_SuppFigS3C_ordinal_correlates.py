"""Supp Fig S3C — Ordinal-DOR correlates (mega-panel, 3 sub-panels).

Purpose:      Render Supp Fig S3C: DOR_ordinal correlations (CN4, Mac PD-L1, baseline PC PD-L1).
Inputs:       imc_extras_persample.csv + pheno_8clus-sample_pct.csv.
Outputs:      figures/SuppFigS3C.pdf / .svg.
Dependencies: matplotlib, scipy.stats, _00_common, _01_data_loaders.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from adjustText import adjust_text

import _00_common as cm
from _01_data_loaders import load_imc_sample_pct, load_imc_extras

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})


def panel_cn4_post():
    cn = load_imc_sample_pct()
    a = cn.groupby(["Pt", "Tp"], as_index=False)["CN4"].mean()
    a = a[a["Tp"] == "PostNivo"].copy()
    a["DOR_ordinal"] = a["Pt"].map(cm.DOR_ORDINAL)
    a["DOR_group"] = a["Pt"].map(cm.DOR)
    return a.dropna(subset=["DOR_ordinal", "CN4"]).rename(columns={"CN4": "x"}), \
           "CN4 post-treatment (% cells)"


def panel_dmac_pdl1():
    df = load_imc_extras().copy()
    piv = df.pivot_table(index="Pt", columns="Tp",
                          values="pct_pdl1_in_macro", aggfunc="first")
    piv["delta"] = piv["PostNivo"] - piv["Baseline"]
    piv["DOR_ordinal"] = piv.index.map(cm.DOR_ORDINAL)
    piv["DOR_group"] = piv.index.map(cm.DOR)
    out = piv.dropna(subset=["delta", "DOR_ordinal"]).reset_index()[
        ["Pt", "delta", "DOR_ordinal", "DOR_group"]]
    return out.rename(columns={"delta": "x"}), \
           r"$\Delta$ Mac PD-L1+ % (post $-$ pre)"


def panel_basepc_pdl1():
    df = load_imc_extras().copy()
    sub = df[df["Tp"] == "Baseline"][["Pt", "pct_pdl1_in_PC"]].dropna()
    sub["DOR_ordinal"] = sub["Pt"].map(cm.DOR_ORDINAL)
    sub["DOR_group"] = sub["Pt"].map(cm.DOR)
    return sub.rename(columns={"pct_pdl1_in_PC": "x"}), \
           "Baseline PD-L1+ plasma cells (%)"


def draw_scatter(ax, df, xlabel):
    rho, p, _ = cm.spearman(df["x"], df["DOR_ordinal"])
    texts = []
    for _, r in df.iterrows():
        color = (cm.COLORS["DR_PostNivo"] if r["DOR_group"] == "DR"
                 else cm.COLORS["P-NR_PostNivo"])
        ax.scatter(r["x"], r["DOR_ordinal"],
                    s=100, alpha=0.9, edgecolor="black", linewidths=0.7,
                    facecolor=color, zorder=3)
        texts.append(ax.text(r["x"], r["DOR_ordinal"], r["Pt"], fontsize=8))

    slope, intercept = np.polyfit(df["x"], df["DOR_ordinal"], 1)
    pad = (df["x"].max() - df["x"].min()) * 0.08 + 0.5
    xs = np.linspace(df["x"].min() - pad, df["x"].max() + pad, 50)
    ax.plot(xs, slope * xs + intercept, c="black", lw=1.0, ls="--", alpha=0.7)

    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("Response (1=SD ... 4=CR)", fontsize=10)
    ax.set_yticks([1, 2, 3, 4])
    ax.set_yticklabels(["SD", "PR", "VGPR", "CR"], fontsize=9)
    ax.set_title(rf"Spearman $\rho={rho:+.2f}$, $p={p:.3f}$",
                  fontsize=10, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    adjust_text(texts, ax=ax, expand=(1.4, 1.6),
                  arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5))
    return rho, p


fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.6))

for ax, builder in [
    (axes[0], panel_cn4_post),
    (axes[1], panel_dmac_pdl1),
    (axes[2], panel_basepc_pdl1),
]:
    df, xlabel = builder()
    rho, p = draw_scatter(ax, df, xlabel)
    print(f"  rho={rho:+.3f}, p={p:.4f}, n={len(df)}")

fig.tight_layout()
cm.save_fig(fig, "SuppFigS3C")
plt.close(fig)
print("wrote SuppFigS3C (3 sub-panels: CN4 ordinal, dMacPDL1, baselinePCPDL1)")
