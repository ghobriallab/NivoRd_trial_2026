"""Box-plot helper used by Fig 1B/1D/1F/Supp boxplots.

Layout per panel (matches original IMC_figures.Rmd aesthetic):
  4 boxes left->right: R-pre, R-post, NR-pre, NR-post.
  Patient-paired connector lines for pre->post within each group.
  Brackets with stars for default comparisons (paired t for same-group
  pre vs post; Welch t for between-group).
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

from _00_common import (DOR, COLORS, GROUP_ORDER, LEGEND_LABELS, pformat)

warnings.filterwarnings("ignore")


def _bracket(ax, x1, x2, y, p, d=None, lw=0.8, fontsize=7.0):
    """Draw a bracket from x1 to x2 at height y with a 'd=…, p=…' label.
    p < 0.05 → bold; p < 1e-4 → scientific notation; NaN → 'n.s.'.
    d is Cohen's d (between-group: pooled-SD; within-group paired: dz =
    mean(diff)/SD(diff)). Sign matches the direction of the t-statistic.
    """
    h = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.012
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], c="black", lw=lw)
    ptext, weight = pformat(p)
    if d is not None and np.isfinite(d):
        text = f"d={d:+.2f}, {ptext}"
    else:
        text = ptext
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom",
            fontsize=fontsize, fontweight=weight)


def boxplot_panel(ax, df, value_col, ylabel, title=None,
                  classifier=None, comparisons=None,
                  ymin=None, ymax=None):
    classifier = classifier if classifier is not None else DOR
    comparisons = comparisons or [
        ("R_Baseline", "R_PostNivo"),
        ("R_PostNivo", "NR_PostNivo"),
        ("NR_Baseline", "NR_PostNivo"),
        ("R_Baseline", "NR_Baseline"),
    ]
    df = df.copy()
    df["Group"] = df["Pt"].map(classifier)
    df = df.dropna(subset=["Group", value_col])
    df["x"] = df["Group"] + "_" + df["Tp"].astype(str)

    positions = {g: i for i, g in enumerate(GROUP_ORDER)}
    box_data = [df.loc[df["x"] == g, value_col].values for g in GROUP_ORDER]

    bp = ax.boxplot(box_data, positions=range(len(GROUP_ORDER)),
                    widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=1.0),
                    whiskerprops=dict(color="black", lw=0.8),
                    capprops=dict(color="black", lw=0.8),
                    boxprops=dict(lw=0.8))
    for patch, g in zip(bp["boxes"], GROUP_ORDER):
        patch.set_facecolor(COLORS[g])
        patch.set_alpha(0.6)
        patch.set_edgecolor("black")

    pts = sorted(df["Pt"].unique())
    for pt in pts:
        sub = df[df["Pt"] == pt]
        for _, row in sub.iterrows():
            ax.scatter(positions[row["x"]], row[value_col],
                       c="black", s=18, zorder=3,
                       edgecolors="white", linewidths=0.4)
        for grp in ("R", "NR"):
            pair = sub[sub["Group"] == grp]
            if len(pair) == 2:
                pre = pair[pair["Tp"] == "Baseline"][value_col]
                post = pair[pair["Tp"] == "PostNivo"][value_col]
                if len(pre) and len(post):
                    ax.plot([positions[f"{grp}_Baseline"],
                             positions[f"{grp}_PostNivo"]],
                            [pre.iloc[0], post.iloc[0]],
                            c="grey", lw=0.7, alpha=0.7, zorder=2)

    all_vals = (np.concatenate([d for d in box_data if len(d)])
                if any(len(d) for d in box_data) else np.array([0.0]))
    if ymin is None:
        ymin = float(all_vals.min() - (all_vals.max() - all_vals.min()) * 0.05)
    if ymax is None:
        ymax = float(all_vals.max() + (all_vals.max() - all_vals.min()) * 0.55)
    ax.set_ylim(ymin, ymax)

    bracket_h = (ymax - ymin) * 0.06
    base_y = all_vals.max() + (ymax - all_vals.max()) * 0.25
    for i, (g1, g2) in enumerate(comparisons):
        v1 = df.loc[df["x"] == g1, value_col].values
        v2 = df.loc[df["x"] == g2, value_col].values
        if len(v1) < 2 or len(v2) < 2:
            continue
        grp1, grp2 = g1.split("_")[0], g2.split("_")[0]
        if grp1 == grp2:
            # paired pre/post within group
            sub = df[(df["Group"] == grp1) &
                     (df["Tp"].isin([g1.split("_")[1], g2.split("_")[1]]))]
            wide = sub.pivot_table(index="Pt", columns="Tp",
                                    values=value_col).dropna()
            if len(wide) >= 2:
                tp1 = g1.split("_")[1]
                tp2 = g2.split("_")[1]
                p = stats.ttest_rel(wide[tp1], wide[tp2])[1]
                # Cohen's dz = mean(diff) / SD(diff)
                diff = (wide[tp1] - wide[tp2]).values
                sd = np.std(diff, ddof=1)
                d_eff = float(np.mean(diff) / sd) if sd > 1e-12 else np.nan
            else:
                continue
        else:
            # unpaired Welch between groups
            p = stats.ttest_ind(v1, v2, equal_var=False)[1]
            n1, n2 = len(v1), len(v2)
            sp = np.sqrt(((n1-1) * np.var(v1, ddof=1) +
                          (n2-1) * np.var(v2, ddof=1)) / (n1 + n2 - 2))
            d_eff = float((np.mean(v1) - np.mean(v2)) / sp) if sp > 1e-12 else np.nan
        x1, x2 = positions[g1], positions[g2]
        _bracket(ax, x1, x2, base_y + i * bracket_h, p, d=d_eff)

    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels(["", "", "", ""])
    ax.tick_params(axis="x", which="both", bottom=False)
    ax.set_xlim(-0.6, len(GROUP_ORDER) - 0.4)
    ax.set_ylabel(ylabel, fontsize=10)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_legend(out_path_stem):
    fig, ax = plt.subplots(figsize=(2.2, 1.6))
    ax.axis("off")
    handles = [Patch(facecolor=COLORS[g], edgecolor="black", alpha=0.6,
                      label=LEGEND_LABELS[g]) for g in GROUP_ORDER]
    ax.legend(handles=handles, loc="center", frameon=False, fontsize=10,
              title="DOR group / timepoint", title_fontsize=10)
    for fmt in ("pdf", "svg"):
        fig.savefig(f"{out_path_stem}.{fmt}", dpi=300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
