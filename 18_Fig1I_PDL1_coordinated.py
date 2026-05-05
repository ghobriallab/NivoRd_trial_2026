"""Fig 1I — Mac PD-L1+ vs PC PD-L1+ scatter (coordinated PD-L1 program).

PC PD-L1 on x-axis (the upstream environmental signal), Mac PD-L1 on
y-axis (the downstream amplification). DR-only and P-NR-only OLS lines
plus an overall reference line. All significance tests sit in a side
panel to the right of the axes (no overlap with regression lines).

Reports:
  - within-group Spearman ρ and OLS slope for DR and P-NR
  - LMM ANCOVA test of slope difference between groups
    (Mac PD-L1 ~ PC PD-L1 × Group + (1|Patient); LRT vs additive)
"""
from __future__ import annotations
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from statsmodels.regression.mixed_linear_model import MixedLM

import _00_common as cm
from _01_data_loaders import load_imc_extras

warnings.filterwarnings("ignore")
mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

df = load_imc_extras().copy()
df["DOR_group"] = df["Pt"].map(cm.DOR)
df["bucket"] = df["DOR_group"] + "_" + df["Tp"]
df = df.dropna(subset=["pct_pdl1_in_macro", "pct_pdl1_in_PC", "DOR_group"])

# ---- linear mixed-effects ANCOVA: Mac ~ PC * Group + (1|Patient) ----
sub = df.copy()
sub["x"]  = sub["pct_pdl1_in_PC"] - sub["pct_pdl1_in_PC"].mean()
sub["g"]  = (sub["DOR_group"] == "P-NR").astype(float)
sub["xg"] = sub["x"] * sub["g"]

def fit(exog_cols):
    X = sub[exog_cols].copy(); X.insert(0, "Intercept", 1.0)
    return MixedLM(sub["pct_pdl1_in_macro"].values, X,
                    groups=sub["Pt"].values).fit(reml=False, method="lbfgs")

m_full = fit(["x", "g", "xg"])
m_add  = fit(["x", "g"])
chi2_int = 2 * (m_full.llf - m_add.llf)
p_int = 1 - stats.chi2.cdf(chi2_int, df=1)

# ---- group-specific stats for the side panel ----
strat = {}
for grp in ("DR", "P-NR"):
    s = df[df["DOR_group"] == grp]
    xg = s["pct_pdl1_in_PC"].values
    yg = s["pct_pdl1_in_macro"].values
    rho_g, p_g = stats.spearmanr(xg, yg)
    slope_g, _ = np.polyfit(xg, yg, 1)
    strat[grp] = (rho_g, p_g, slope_g)

# ---- figure layout: axes on left, stat box on right ----
fig = plt.figure(figsize=(8.4, 3.4))
ax = fig.add_axes([0.08, 0.16, 0.42, 0.78])  # [left, bottom, w, h]

# scatter by DR/P-NR × pre/post
for bk in cm.GROUP_ORDER:
    s = df[df["bucket"] == bk]
    ax.scatter(s["pct_pdl1_in_PC"], s["pct_pdl1_in_macro"],
                s=70, alpha=0.85, edgecolor="black", linewidths=0.6,
                facecolor=cm.COLORS[bk], label=cm.LEGEND_LABELS[bk], zorder=3)

# overall reference line (faint)
slope_all, b_all = np.polyfit(df["pct_pdl1_in_PC"], df["pct_pdl1_in_macro"], 1)
xs_all = np.linspace(df["pct_pdl1_in_PC"].min(), df["pct_pdl1_in_PC"].max(), 50)
ax.plot(xs_all, slope_all * xs_all + b_all,
         c="black", lw=1.0, ls="--", alpha=0.45, zorder=2)

# group lines
for grp, color in (("DR",   cm.COLORS["DR_PostNivo"]),
                    ("P-NR", cm.COLORS["P-NR_PostNivo"])):
    s = df[df["DOR_group"] == grp]
    xg = s["pct_pdl1_in_PC"].values
    yg = s["pct_pdl1_in_macro"].values
    sl, ic = np.polyfit(xg, yg, 1)
    xline = np.linspace(xg.min(), xg.max(), 30)
    ax.plot(xline, sl * xline + ic, c=color, lw=1.6, alpha=0.85, zorder=2)

ax.set_xlabel("PD-L1+ plasma cells (%)", fontsize=10)
ax.set_ylabel("PD-L1+ macrophages (%)", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(fontsize=7, loc="lower right", frameon=False, title=None)

# ---- side panel with all stats (one line per category) ----
ax_stats = fig.add_axes([0.52, 0.16, 0.46, 0.78])
ax_stats.axis("off")

def fmt_p(p):
    if not np.isfinite(p): return ("n.s.", "normal")
    bold = p < 0.05
    if p < 1e-3:
        e = int(np.floor(np.log10(p))); m = p / 10**e
        body = rf"p={m:.1f}\times 10^{{{e}}}"
        text = rf"$\mathbf{{{body}}}$" if bold else f"${body}$"
    elif p < 0.1:
        text = f"p={p:.3f}"
    else:
        text = f"p={p:.2f}"
    return text, ("bold" if bold else "normal")

COL_DR  = cm.COLORS["DR_PostNivo"]
COL_PNR = cm.COLORS["P-NR_PostNivo"]

# DR line
rho_DR, p_DR, slope_DR = strat["DR"]
p_text_DR, w_DR = fmt_p(p_DR)
ax_stats.text(0, 0.85,
               rf"$\bf{{DR:}}$  $\rho={rho_DR:+.2f}$,  $\bf{{slope={slope_DR:+.2f}}}$,  {p_text_DR}",
               transform=ax_stats.transAxes, fontsize=10, color=COL_DR,
               va="center", ha="left", fontweight=w_DR)

# P-NR line
rho_N, p_N, slope_N = strat["P-NR"]
p_text_N, w_N = fmt_p(p_N)
ax_stats.text(0, 0.65,
               rf"$\bf{{P/NR:}}$  $\rho={rho_N:+.2f}$,  $\bf{{slope={slope_N:+.2f}}}$,  {p_text_N}",
               transform=ax_stats.transAxes, fontsize=10, color=COL_PNR,
               va="center", ha="left", fontweight=w_N)

# Interaction (slope-difference) line
p_text_int, w_int = fmt_p(p_int)
ax_stats.text(0, 0.40,
               rf"$\bf{{DR/P/NR\ interaction\ (slope\ difference):}}$",
               transform=ax_stats.transAxes, fontsize=10, color="black",
               va="center", ha="left")
ax_stats.text(0, 0.25,
               rf"LMM ANCOVA (LRT, df=1):  $\chi^2={chi2_int:.1f}$,  {p_text_int}",
               transform=ax_stats.transAxes, fontsize=10, color="black",
               va="center", ha="left", fontweight=w_int)
# slope unit footnote
ax_stats.text(0, 0.05, "(slope = Δ Mac PD-L1+ % per +1 unit PC PD-L1+ %)",
               transform=ax_stats.transAxes, fontsize=8, color="0.45",
               va="center", ha="left", fontstyle="italic")

cm.save_fig(fig, "Fig1I")
plt.close(fig)
print(f"wrote Fig1I  DR: ρ={strat['DR'][0]:+.3f}, slope={strat['DR'][2]:+.3f}, p={strat['DR'][1]:.4f}; "
      f"P-NR: ρ={strat['P-NR'][0]:+.3f}, slope={strat['P-NR'][2]:+.3f}, p={strat['P-NR'][1]:.4f}; "
      f"interaction χ²={chi2_int:.2f}, p={p_int:.2e}")
