"""ANCOVA / mixed-effects test for the Fig 1I PD-L1 coordination panel.

Purpose:      LMM ANCOVA for Mac PD-L1 vs PC PD-L1 x Group on the IMC cohort.
Inputs:       imc_extras_persample.csv (NivoRd cohort).
Outputs:      Console summary + LRT chi-square / p-values used by Figure 1I.
Dependencies: statsmodels MixedLM, pandas, _00_common.
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.regression.mixed_linear_model import MixedLM

import _00_common as cm
from _01_data_loaders import load_imc_extras

warnings.filterwarnings("ignore")

df = load_imc_extras()[
    ["Pt", "Tp", "pct_pdl1_in_macro", "pct_pdl1_in_PC"]
].dropna().copy()
df["Group"] = df["Pt"].map(cm.DOR)
df = df.dropna(subset=["Group"])
df["Group_PNR"] = (df["Group"] == "P-NR").astype(float)
df["x_centered"] = df["pct_pdl1_in_PC"] - df["pct_pdl1_in_PC"].mean()
df["x_times_PNR"] = df["x_centered"] * df["Group_PNR"]

print(f"n samples: {len(df)}  ({(df.Group=='DR').sum()} DR, {(df.Group=='P-NR').sum()} P-NR)")
print(f"n patients: {df.Pt.nunique()}")
print()

# 1) OLS within each group for sanity
for g in ["DR", "P-NR"]:
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
exog_full = df[["x_centered", "Group_PNR", "x_times_PNR"]].copy()
exog_full.insert(0, "Intercept", 1.0)
exog_add  = df[["x_centered", "Group_PNR"]].copy()
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
p_int_wald   = m_full.pvalues["x_times_PNR"]
p_grp_wald   = m_full.pvalues["Group_PNR"]
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
beta_grp = m_add.fe_params["Group_PNR"]
print(f"\nP-NR macrophages have {beta_grp:+.2f} %-points more PD-L1+ "
      f"at matched PC PD-L1 (additive model).")
