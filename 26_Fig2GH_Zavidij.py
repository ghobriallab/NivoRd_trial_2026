"""Fig 2G + Fig 2H — Zavidij (GSE124310) external validation.

Two outputs: Fig2G.{pdf,svg} (IFN-γ) and Fig2H.{pdf,svg} (TNFα-NFkB).
Adapted from 2-Scripts/06-Figure2EF_external_validation.py.
"""
from __future__ import annotations
import os, ftplib, tarfile, shutil, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.io import mmread
from scipy.stats import mannwhitneyu
import gseapy

import _00_common as cm

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

DATA = cm.OUT_DIR.parent / "data" / "zavidij"
EXTRACT = DATA / "extracted"
PER_SAMPLE = DATA / "per_sample"
CONCAT = DATA / "zavidij_concat.h5ad"
ANNOT = DATA / "zavidij_annotated.h5ad"

GEO_HOST = "ftp.ncbi.nlm.nih.gov"
GEO_PATH = "/geo/series/GSE124nnn/GSE124310/suppl/"
STAGE_ORDER = ["HD", "MGUS", "SMMl", "SMMh", "MM"]
STAGE_COLORS = {"HD": "#95a5a6", "MGUS": "#f1c40f", "SMMl": "#3498db",
                 "SMMh": "#2980b9", "MM": "#e74c3c"}
MARKERS = {
    "T_cell":       ["CD3D", "CD3E", "CD3G"],
    "CD4_T":        ["CD4", "IL7R", "LEF1"],
    "CD8_T":        ["CD8A", "CD8B"],
    "NK":           ["NCAM1", "NKG7", "KLRD1", "KLRF1"],
    "B_cell":       ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "Plasma":       ["SDC1", "MZB1", "XBP1"],
    "Macrophage":   ["CD68", "C1QA", "C1QB", "CD163"],
    "Monocyte_cl":  ["CD14", "VCAN", "S100A8", "S100A9", "LYZ"],
    "Monocyte_ncl": ["FCGR3A", "MS4A7", "CSF1R"],
    "DC":           ["FCER1A", "CLEC10A", "CD1C", "HLA-DRA"],
    "Erythroid":    ["HBB", "HBA1", "HBA2"],
    "Progenitor":   ["CD34", "MPO", "ELANE"],
}


def download_geo():
    EXTRACT.mkdir(parents=True, exist_ok=True)
    if [f for f in EXTRACT.iterdir() if f.suffix == ".gz"]:
        return
    ftp = ftplib.FTP(GEO_HOST); ftp.login(); ftp.cwd(GEO_PATH)
    files = [f for f in ftp.nlst() if f.endswith(".tar.gz")]
    for f in files:
        with open(EXTRACT / f, "wb") as fh:
            ftp.retrbinary(f"RETR {f}", fh.write)
    ftp.quit()


def parse_name(fname):
    base = fname.replace(".filtered_gene_bc_matrices.tar.gz", "")
    gsm, rest = base.split("_", 1)
    sample, fraction = rest.split(".")
    stage = ("HD" if sample.startswith("NBM") else
             "MGUS" if sample.startswith("MGUS") else
             "SMMh" if sample.startswith("SMMh") else
             "SMMl" if sample.startswith("SMMl") else
             "MM" if sample.startswith("MM") else "Other")
    return gsm, sample, fraction, stage


def load_sample(fname):
    gsm, sample, fraction, stage = parse_name(fname)
    sdir = PER_SAMPLE / f"{sample}_{fraction}"
    if sdir.exists():
        shutil.rmtree(sdir)
    sdir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(EXTRACT / fname, "r:gz") as t:
        t.extractall(sdir)
    mdir = glob.glob(str(sdir / "filtered_gene_bc_matrices/GRCh38"))[0]
    X = mmread(os.path.join(mdir, "matrix.mtx")).T.tocsr().astype(np.float32)
    bc = pd.read_csv(os.path.join(mdir, "barcodes.tsv"),
                       header=None, sep="\t")[0].values
    g = pd.read_csv(os.path.join(mdir, "genes.tsv"),
                       header=None, sep="\t",
                       names=["gene_id", "gene_symbol"])
    a = ad.AnnData(X=X,
                    obs=pd.DataFrame({"barcode": bc}, index=bc),
                    var=pd.DataFrame({"gene_id": g.gene_id.values,
                                        "gene_symbol": g.gene_symbol.values},
                                      index=g.gene_symbol.values))
    a.var_names_make_unique()
    a.obs[["gsm", "sample", "fraction", "disease_stage"]] = gsm, sample, fraction, stage
    a.obs_names = [f"{sample}_{b}" for b in a.obs_names]
    return a


def load_and_qc():
    if CONCAT.exists():
        return sc.read_h5ad(CONCAT)
    PER_SAMPLE.mkdir(parents=True, exist_ok=True)
    files = sorted(f.name for f in EXTRACT.iterdir() if f.suffix == ".gz")
    parts = [load_sample(f) for f in files]
    a = ad.concat(parts, join="inner", label="sample_key", index_unique=None)
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None,
                                  log1p=False, inplace=True)
    a = a[(a.obs.n_genes_by_counts >= 500)
           & (a.obs.n_genes_by_counts <= 5000)
           & (a.obs.pct_counts_mt < 15)].copy()
    sc.pp.filter_genes(a, min_cells=3)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.write_h5ad(CONCAT, compression="gzip")
    return a


def annotate(a):
    if ANNOT.exists():
        return sc.read_h5ad(ANNOT)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat",
                                  batch_key="sample")
    a.raw = a
    hv = a[:, a.var.highly_variable].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=30, svd_solver="arpack")
    a.obsm["X_pca"] = hv.obsm["X_pca"]
    try:
        import harmonypy as hm
        Z = np.asarray(hm.run_harmony(a.obsm["X_pca"], a.obs, "sample",
                                         max_iter_harmony=20).Z_corr)
        a.obsm["X_pca_harmony"] = Z if Z.shape[0] == a.n_obs else Z.T
        rep = "X_pca_harmony"
    except ImportError:
        rep = "X_pca"
    sc.pp.neighbors(a, n_neighbors=30, use_rep=rep)
    sc.tl.umap(a)
    sc.tl.leiden(a, resolution=0.6)
    for ct, genes in MARKERS.items():
        present = [g for g in genes if g in a.var_names]
        if present:
            sc.tl.score_genes(a, gene_list=present,
                                 score_name=f"score_{ct}", random_state=0)
        else:
            a.obs[f"score_{ct}"] = 0.0
    score_cols = [c for c in a.obs.columns if c.startswith("score_")]
    by_clu = a.obs.groupby("leiden")[score_cols].mean()
    a.obs["celltype"] = a.obs["leiden"].map(
        {c: by_clu.loc[c].idxmax().replace("score_", "")
         for c in by_clu.index})
    mac_m = [g for g in ["CD68", "C1QA", "C1QB"] if g in a.var_names]
    mono_m = [g for g in ["S100A8", "S100A9", "VCAN"] if g in a.var_names]
    if mac_m and mono_m:
        sc.tl.score_genes(a, mac_m, score_name="score_Mac_tight", random_state=0)
        sc.tl.score_genes(a, mono_m, score_name="score_Mono_tight", random_state=0)
        a.obs["is_macrophage"] = (
            a.obs.celltype.isin(["Macrophage", "Monocyte_cl", "Monocyte_ncl"])
            & (a.obs.score_Mac_tight > a.obs.score_Mono_tight))
    else:
        a.obs["is_macrophage"] = a.obs.celltype == "Macrophage"
    a.write_h5ad(ANNOT, compression="gzip")
    return a


def score_macs(a):
    mac = a[a.obs.is_macrophage].copy()
    lib = gseapy.get_library(name="MSigDB_Hallmark_2020", organism="human")
    sets = {"HALLMARK_INTERFERON_GAMMA_RESPONSE": "Interferon Gamma Response",
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF-alpha Signaling via NF-kB"}
    for name, key in sets.items():
        found = key if key in lib else next(
            k for k in lib if ("Interferon Gamma" if "INTERFERON" in name else "TNF") in k)
        genes = [g for g in lib[found] if g in mac.var_names]
        sc.tl.score_genes(mac, gene_list=genes, score_name=name,
                             random_state=0)
    return mac.obs[["sample", "disease_stage",
                     "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                     "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / (np.arange(n) + 1)
    for i in range(n - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    out = np.empty_like(q); out[o] = np.clip(q, 0, 1)
    return out


def fmt_q(q):
    if q == 0: return "q<1e-300"
    e = int(np.floor(np.log10(q))); m = q / (10 ** e)
    return f"q={m:.1f}×10$^{{{e}}}$"


def draw(ax, df, metric, title):
    data = [df.loc[df.disease_stage == s, metric].values for s in STAGE_ORDER]
    pos = list(range(len(STAGE_ORDER)))
    bp = ax.boxplot(data, positions=pos, widths=0.6, patch_artist=True,
                     showfliers=False)
    for patch, s in zip(bp["boxes"], STAGE_ORDER):
        patch.set_facecolor(STAGE_COLORS[s]); patch.set_alpha(0.8)
        patch.set_edgecolor("#222222"); patch.set_linewidth(1.0)
    for m in bp["medians"]: m.set_color("#111111"); m.set_linewidth(1.8)
    for w in bp["whiskers"]: w.set_color("#444444"); w.set_linewidth(1.0)
    for c in bp["caps"]: c.set_color("#444444"); c.set_linewidth(1.0)
    ax.set_xticks(pos); ax.set_xticklabels(STAGE_ORDER, fontsize=11)
    ax.tick_params(axis="y", labelsize=10)
    ax.set_title(f"{title}\nZavidij et al. · GSE124310", fontsize=12, pad=8)
    ax.set_ylabel("Module score", fontsize=11)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    hd = df.loc[df.disease_stage == "HD", metric].values
    ps = [mannwhitneyu(df.loc[df.disease_stage == s, metric].values, hd,
                          alternative="two-sided").pvalue
          for s in STAGE_ORDER[1:]]
    qs = bh(ps)
    allv = np.concatenate(data)
    lo, hi = float(allv.min()), float(allv.max()); span = hi - lo
    base = hi + 0.07 * span; step = 0.13 * span
    for i, q in enumerate(qs, 1):
        y = base + (i - 1) * step
        ax.plot([0, 0, i, i], [y, y + 0.022 * span, y + 0.022 * span, y],
                 color="#222222", lw=1.1)
        ax.text(i / 2, y + 0.028 * span, fmt_q(q), ha="center",
                  va="bottom", fontsize=10)
    ax.set_ylim(lo - 0.05 * span, base + len(STAGE_ORDER) * step)


def main():
    download_geo()
    adata = load_and_qc()
    adata = annotate(adata)
    df = score_macs(adata)
    panels = [("HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "Hallmark IFN-γ response", "Fig2G"),
              ("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
               "Hallmark TNF-α via NF-κB", "Fig2H")]
    for metric, title, tag in panels:
        fig, ax = plt.subplots(1, 1, figsize=(4.4, 4.4))
        draw(ax, df, metric, title)
        plt.tight_layout()
        cm.save_fig(fig, tag)
        plt.close(fig)
    print("wrote Fig2G and Fig2H")


if __name__ == "__main__":
    main()
