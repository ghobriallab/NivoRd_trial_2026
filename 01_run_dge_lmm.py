"""Per-cell linear mixed model DGE for macrophages, post-treatment and

Purpose:      Per-cell linear mixed-effects DGE (DR vs P-NR) on pre- and post-treatment macrophages.
Inputs:       CD138neg_immune_cells.h5ad (NivoRd cohort macrophages).
Outputs:      {macrophage_pretreatment,macrophage_posttreatment}_DGE_LMM_DOR.csv in OUT_DIR.
Dependencies: anndata, statsmodels MixedLM, scipy.stats, pandas, _00_common.
"""
from __future__ import annotations
import warnings
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed

import _00_common as cm
from _01_data_loaders import load_scrna_immune

warnings.filterwarnings("ignore")

CELLTYPE_COL = "celltype_refined"
CELLTYPE = "Macrophages"
MIN_PCT = 0.05
MIN_CELLS_PER_PATIENT = 5
N_JOBS = 4

SEX_GENES = {"XIST", "RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "UTY",
             "ZFY", "USP9Y", "TTTY15"}


def fit_gene(gene, expr, group, patient):
    try:
        exog = pd.DataFrame({"Intercept": 1.0, "Group": group})
        res = MixedLM(endog=expr, exog=exog, groups=patient).fit(
            reml=True, method="lbfgs", maxiter=200)
        return {"Gene": gene, "coef": res.fe_params["Group"],
                "se": res.bse_fe["Group"], "zstat": res.tvalues["Group"],
                "pvalue": res.pvalues["Group"], "converged": res.converged}
    except Exception:
        return None


def run_one(adata, tp_label, out_csv, summary_lines):
    a = adata[(adata.obs[CELLTYPE_COL].astype(str) == CELLTYPE)
              & (adata.obs["Tp"] == tp_label)].copy()
    a = a[a.obs["DOR_group"].isin(["DR", "P-NR"])].copy()
    counts = a.obs.groupby("Pt").size()
    keep = counts[counts >= MIN_CELLS_PER_PATIENT].index
    a = a[a.obs["Pt"].isin(keep)].copy()
    pts = sorted(a.obs["Pt"].unique().tolist())
    grp_counts = a.obs.groupby("DOR_group").size().to_dict()
    line = (f"  [{tp_label}] cells={a.n_obs}, patients={pts}, "
            f"DR={grp_counts.get('DR',0)}, P-NR={grp_counts.get('P-NR',0)}")
    print(line, flush=True)
    summary_lines.append(line)

    if a.obs["DOR_group"].nunique() < 2 or len(pts) < 2:
        msg = "  insufficient groups/patients; skipping"
        print(msg, flush=True); summary_lines.append(msg)
        return None

    X = a.X.toarray() if issparse(a.X) else np.asarray(a.X)
    X = np.asarray(X, dtype=np.float64)
    # Group_PNR = 1 for P-NR, 0 for DR; positive coef => higher in P-NR
    group = (a.obs["DOR_group"] == "P-NR").astype(float).values
    patient = a.obs["Pt"].values
    var = np.array(a.var_names)

    keep_genes = np.zeros(X.shape[1], dtype=bool)
    for g_label in ["DR", "P-NR"]:
        m = a.obs["DOR_group"].values == g_label
        keep_genes |= ((X[m] > 0).mean(axis=0) >= MIN_PCT)
    genes = [g for g in var[keep_genes]
             if g not in SEX_GENES and not str(g).startswith("MT-")]
    idx = {g: i for i, g in enumerate(var)}
    print(f"  testing {len(genes)} genes", flush=True)

    results = Parallel(n_jobs=N_JOBS, batch_size=200)(
        delayed(fit_gene)(g, X[:, idx[g]], group, patient) for g in genes)
    results = [r for r in results if r is not None]

    df = pd.DataFrame(results).set_index("Gene")
    df["log2FoldChange"] = df["coef"] / np.log(2)
    pvals = df["pvalue"].values
    valid = ~np.isnan(pvals)
    padj = np.full(len(pvals), np.nan)
    if valid.sum() > 0:
        _, padj[valid], _, _ = multipletests(pvals[valid], method="fdr_bh")
    df["padj"] = padj
    df = df.sort_values("padj")
    df.to_csv(out_csv)

    sig_all = df[(df.padj < 0.05) & (df.log2FoldChange.abs() > 0.5)]
    sig = sig_all[sig_all["converged"].astype(bool)]
    n_up_PNR = int((sig.log2FoldChange > 0).sum())
    n_up_DR = int((sig.log2FoldChange < 0).sum())
    n_excluded = len(sig_all) - len(sig)
    line = (f"  -> {len(sig)} sig (FDR<0.05, |LFC|>0.5, converged): "
            f"up_DR={n_up_DR}, up_PNR={n_up_PNR} "
            f"({n_excluded} non-converged dropped); CSV: {out_csv.name}")
    print(line, flush=True); summary_lines.append(line)
    return df


def main():
    print("Loading scRNA immune h5ad…", flush=True)
    a = load_scrna_immune(nivord_only=True)
    print(f"  NivoRd cells: {a.n_obs}; cohort: "
          f"{sorted(a.obs['Pt'].dropna().unique().tolist())}",
          flush=True)
    summary_lines = [f"Macrophage def: {CELLTYPE_COL}=='{CELLTYPE}'",
                      f"Min cells/pt: {MIN_CELLS_PER_PATIENT}; "
                      f"min %expr: {MIN_PCT*100:.0f}%; "
                      f"BH-FDR; thresholds: FDR<0.05, |log2FC|>0.5",
                      f"Cohort under DOR: DR={[p for p,v in cm.DOR.items() if v=='DR' and p in cm.SCRNA_COHORT]}, "
                      f"P-NR={[p for p,v in cm.DOR.items() if v=='P-NR' and p in cm.SCRNA_COHORT]}",
                      ""]

    print()
    print("=== POST-TREATMENT ===", flush=True)
    summary_lines.append("=== POST-TREATMENT ===")
    out_post = cm.OUT_DIR / "macrophage_posttreatment_DGE_LMM_DOR.csv"
    run_one(a, "PostNivo", out_post, summary_lines)

    print()
    print("=== PRE-TREATMENT ===", flush=True)
    summary_lines.append("\n=== PRE-TREATMENT ===")
    out_pre = cm.OUT_DIR / "macrophage_pretreatment_DGE_LMM_DOR.csv"
    run_one(a, "Baseline", out_pre, summary_lines)

    (cm.OUT_DIR / "dge_summary.txt").write_text("\n".join(summary_lines))
    print()
    print("DONE")


if __name__ == "__main__":
    main()
