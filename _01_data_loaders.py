"""Data loaders for the manuscript_v4_DOR pipeline.

Purpose:      Canonical loaders for IMC h5ad, scRNA-seq immune/tumor h5ads, IMC extras CSV.
Inputs:       h5ads + CSVs under DATA_LOCAL (resolved via NIVO_DATA_DIR or _00_common defaults).
Outputs:      None (returns AnnData / pandas objects).
Dependencies: anndata, pandas, _00_common.
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
import anndata as ad

from _00_common import (
    H5AD_IMC, H5AD_IMMUNE, H5AD_TUMOR,
    EXTRAS_CSV, SAMPLE_PCT_CSV, PHENO_CSV,
    IMC_PT_TO_CANON, SCRNA_PT_TO_CANON, SCRNA_POST_TO_CANON,
    DOR, COT, DOR_ORDINAL, ORDINAL_RESPONSE, SCRNA_COHORT,
)

warnings.filterwarnings("ignore")


def _add_canonical_imc(obs):
    obs = obs.copy()
    obs["Pt"] = obs["pt"].map(IMC_PT_TO_CANON)
    obs["Tp"] = obs["tp"].astype(str)
    obs["DOR_group"] = obs["Pt"].map(DOR)
    obs["COT_group"] = obs["Pt"].map(COT)
    obs["DOR_ordinal"] = obs["Pt"].map(DOR_ORDINAL)
    return obs


def _add_canonical_scrna(obs):
    obs = obs.copy()
    pt = obs["patient"].astype(str)
    canon = pt.map(SCRNA_PT_TO_CANON).fillna(pt.map(SCRNA_POST_TO_CANON))
    obs["Pt"] = canon
    src = obs["Source"].astype(str)
    obs["Tp"] = np.where(src == "Nivo_post", "PostNivo",
                np.where(src == "Nivo", "Baseline", src))
    obs["DOR_group"] = obs["Pt"].map(DOR)
    obs["COT_group"] = obs["Pt"].map(COT)
    obs["DOR_ordinal"] = obs["Pt"].map(DOR_ORDINAL)
    return obs


def load_imc(attach_pheno=True):
    a = ad.read_h5ad(H5AD_IMC)
    a.obs_names_make_unique()
    a.obs = _add_canonical_imc(a.obs)
    if attach_pheno and PHENO_CSV.exists():
        ph = pd.read_csv(PHENO_CSV)
        if len(ph) == a.n_obs:
            a.obs["pheno"] = ph["pheno"].values
            if "leiden" in ph.columns:
                a.obs["leiden"] = ph["leiden"].astype(str).values
    return a


def load_scrna_immune(nivord_only=True, exclude_doublets=True):
    a = ad.read_h5ad(H5AD_IMMUNE)
    if hasattr(a.obs.index, "categories"):
        a.obs.index = a.obs.index.astype(str)
    a.obs_names_make_unique()
    a.obs = _add_canonical_scrna(a.obs)
    keep = np.ones(a.n_obs, dtype=bool)
    if nivord_only:
        keep &= a.obs["Source"].isin(["Nivo", "Nivo_post"]).values
    if exclude_doublets and "Immune_cell_2" in a.obs.columns:
        keep &= a.obs["Immune_cell_2"].astype(str).values != "Doublets"
    return a[keep].copy()


def load_scrna_tumor(nivord_only=True, malignant_only=True):
    a = ad.read_h5ad(H5AD_TUMOR)
    if hasattr(a.obs.index, "categories"):
        a.obs.index = a.obs.index.astype(str)
    a.obs_names_make_unique()
    a.obs = _add_canonical_scrna(a.obs)
    keep = np.ones(a.n_obs, dtype=bool)
    if nivord_only:
        keep &= (a.obs["Source"].astype(str).values == "Nivo")
    if malignant_only and "Malignant_PC_Immune" in a.obs.columns:
        is_mal = (
            (a.obs["Malignant_PC_Immune"].astype(str).values == "Malignant_PC")
            | (a.obs.get("InferCNV", pd.Series(["?"] * a.n_obs)).astype(str).values == "Malignant"))
        keep &= is_mal
    return a[keep].copy()


def patient_response_meta():
    """One row per patient with Pt, response (CR/VGPR/PR/SD), DOR_ordinal,
    DOR_group, COT_group. Source: IMC h5ad obs.reponse (note misspelling)."""
    a = ad.read_h5ad(H5AD_IMC)
    obs = _add_canonical_imc(a.obs)
    cols = ["Pt", "DOR_ordinal", "DOR_group", "COT_group", "reponse"]
    if "exceptional" in obs.columns:
        cols.append("exceptional")
    meta = obs[cols].drop_duplicates("Pt").sort_values("Pt").reset_index(drop=True)
    return meta


def load_imc_extras():
    return pd.read_csv(EXTRAS_CSV)


def load_imc_sample_pct():
    return pd.read_csv(SAMPLE_PCT_CSV)
