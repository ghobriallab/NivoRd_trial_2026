"""ANCOVA / mixed-effects test for the Fig 1I PD-L1 coordination panel.

Q1 (interaction): does the Mac-PC PD-L1 slope differ between R and NR?
Q2 (intercept):   after accounting for that slope, is Mac PD-L1
                  systematically higher in NR at matched PC PD-L1?

Model: Mac_PDL1 ~ PC_PDL1 * Group + (1 | Patient)

Reports:
- the interaction p (test of slope difference between groups)
- the Group main effect p (intercept difference at matched PC_PDL1)
- the LRT chi-square comparing full vs. additive models
- both R-only and NR-only OLS slopes for sanity
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.regression.mixed_linear_model import MixedLM

import _common as cm
from _data_loaders import load_imc_extras

warnings.filterwarnings("ignore")

df = load_imc_extras()[
    ["Pt", "Tp", "pct_pdl1_in_macro", "pct_pdl1_in_PC"]
].dropna().copy()
df["Group"] = df["Pt"].map(cm.DOR)
df = df.dropna(subset=["Group"])
df["Group_NR"] = (df["Group"] == "NR").astype(float)
df["x_centered"] = df["pct_pdl1_in_PC"] - df["pct_pdl1_in_PC"].mean()
df["x_times_NR"] = df["x_centered"] * df["Group_NR"]

print(f"n samples: {len(df)}  ({(df.Group=='R').sum()} R, {(df.Group=='NR').sum()} NR)")
print(f"n patients: {df.Pt.nunique()}")
print()

# 1) OLS within each group for sanity
for g in ["R", "NR"]:
    sub = df[df.Group == g]
    rho, p_rho = stats.spearmanr(sub.pct_pdl1_in_PC, sub.pct_pdl1_in_macro)
    slope, intercept, r, p, se = stats.linregress(
        sub.pct_pdl1_in_PC, sub.pct_pdl1_in_macro)
    print(f"  {g}-only (n={len(sub)}): "
          f"Spearman rho={rho:+.3f}, p={p_rho:.3f}; "
          f"OLS slope={slope:+.3f} +- {se:.3f}, intercept={intercept:+.2f}")
print()

# 2) Mixed-effects: full model (with interaction)
endog = df["pct_pdl1_in_macro"].values
exog_full = df[["x_centered", "Group_NR", "x_times_NR"]].copy()
exog_full.insert(0, "Intercept", 1.0)
exog_add  = df[["x_centered", "Group_NR"]].copy()
exog_add.insert(0, "Intercept", 1.0)
exog_null = df[["x_centered"]].copy()
exog_null.insert(0, "Intercept", 1.0)

groups = df["Pt"].values

m_full = MixedLM(endog, exog_full, groups=groups).fit(reml=False, method="lbfgs")
m_add  = MixedLM(endog, exog_add,  groups=groups).fit(reml=False, method="lbfgs")
m_null = MixedLM(endog, exog_null, groups=groups).fit(reml=False, method="lbfgs")

print("=== Mixed-effects model: Mac_PDL1 ~ PC_PDL1*Group + (1|Patient) ===")
print(f"Fixed effects (full model):")
print(m_full.summary().tables[1])
print()

# Wald p-values from the full model
p_int_wald   = m_full.pvalues["x_times_NR"]
p_grp_wald   = m_full.pvalues["Group_NR"]
print(f"Wald p (interaction PC_PDL1 x Group):       {p_int_wald:.4f}")
print(f"Wald p (Group main effect at mean PC_PDL1): {p_grp_wald:.4f}")
print()

# LRTs
lr_int  = 2 * (m_full.llf - m_add.llf)
p_int_lrt = 1 - stats.chi2.cdf(lr_int, df=1)
lr_grp  = 2 * (m_add.llf - m_null.llf)
p_grp_lrt = 1 - stats.chi2.cdf(lr_grp, df=1)
print(f"LRT chi-square interaction:   {lr_int:.3f}, df=1, p={p_int_lrt:.4f}")
print(f"LRT chi-square group (after adjusting for PC_PDL1, no interaction): "
      f"{lr_grp:.3f}, df=1, p={p_grp_lrt:.4f}")
print()

# 3) Effect size on the Group main effect from the additive model
print("=== Additive (no-interaction) model: ===")
print(m_add.summary().tables[1])
beta_grp = m_add.fe_params["Group_NR"]
print(f"\nNR macrophages have {beta_grp:+.2f} %-points more PD-L1+ "
      f"at matched PC PD-L1 (additive model).")
