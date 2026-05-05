"""Shared helpers for the manuscript_v4_DOR figure pipeline.

DOR = IMWG depth of response. Primary classifier for this manuscript:
  DR    = deep responders, VGPR or CR best response (n=4)
  P-NR  = partial / non-responders, PR or SD best response (n=4)

Per IMWG criteria, PR is a recognized response; the binary contrast in
this manuscript reflects depth of response (>=VGPR vs <VGPR), not
presence vs absence of response. Code-side identifiers use the
hyphenated `P-NR`; prose / figure legends render this as `P/NR`.

Ordinal scale (SD=1, PR=2, VGPR=3, CR=4) is the additional analytic frame
used for one pre-planned Spearman test per axis.

For sensitivity / robustness reporting only, COT (manuscript's published
change-of-treatment-within-24-mo rule) and FAST (biochemical PD within
12 mo) are also defined; their group labels follow the same DR / P-NR
convention to keep downstream code uniform.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# All data files are downloaded from Zenodo into ``data/`` next to this
# module by ``download_data.py``. Outputs and figures go into local dirs
# (created on import). The ``NIVO_DATA_DIR`` environment variable can
# override the data directory location.
import os

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA_LOCAL = Path(os.environ.get("NIVO_DATA_DIR", str(HERE / "data")))
OUT_DIR = HERE / "outputs"
LOG_DIR = HERE / "logs"
FIG_DIR = ROOT / "figures" if (ROOT / "figures").parent.exists() else HERE / "figures"
for d in (OUT_DIR, LOG_DIR, FIG_DIR):
    d.mkdir(parents=True, exist_ok=True)

# All data inputs live under DATA_LOCAL after running download_data.py
H5AD_IMC    = DATA_LOCAL / "IMC-import.h5ad"
H5AD_IMMUNE = DATA_LOCAL / "CD138neg_immune_cells.h5ad"
H5AD_TUMOR  = DATA_LOCAL / "CD138pos_tumor_cells.h5ad"

PHENO_CSV      = DATA_LOCAL / "phenotyped_canonical.csv"
EXTRAS_CSV     = DATA_LOCAL / "imc_extras_persample.csv"
SAMPLE_PCT_CSV = DATA_LOCAL / "pheno_8clus-sample_pct.csv"
CN_CENTERS_CSV = DATA_LOCAL / "pheno_8clus-cn_centers.csv"
GMT_PATH       = DATA_LOCAL / "MSigDB_Hallmark_2020.gmt"

# ---------------------------------------------------------------------------
# Classifiers
# ---------------------------------------------------------------------------
DOR = {  # primary classifier for this manuscript: VGPR/CR -> DR; SD/PR -> P-NR
    "Pt05": "DR", "Pt06": "DR", "Pt07": "DR", "Pt09": "DR",
    "Pt01": "P-NR", "Pt03": "P-NR", "Pt04": "P-NR", "Pt08": "P-NR",
}
COT = {  # change-of-treatment within 24 mo (published, kept for robustness only)
    "Pt05": "DR", "Pt06": "DR", "Pt07": "DR", "Pt08": "DR", "Pt09": "DR",
    "Pt01": "P-NR", "Pt03": "P-NR", "Pt04": "P-NR",
}
FAST = {  # biochemical PD within 12 mo (alternative classifier)
    "Pt06": "DR", "Pt07": "DR", "Pt08": "DR", "Pt09": "DR",
    "Pt01": "P-NR", "Pt03": "P-NR", "Pt04": "P-NR", "Pt05": "P-NR",
}

ORDINAL_RESPONSE = {"SD": 1, "PR": 2, "VGPR": 3, "CR": 4}
ORDINAL_LABEL    = {1: "SD", 2: "PR", 3: "VGPR", 4: "CR"}

# Per-Pt ordinal score (uses IMWG best response from CLASSIFIER_AUDIT.xlsx)
DOR_ORDINAL = {
    "Pt01": 1,  # SD
    "Pt03": 2,  # PR
    "Pt04": 2,  # PR
    "Pt05": 3,  # VGPR
    "Pt06": 4,  # CR
    "Pt07": 3,  # VGPR
    "Pt08": 2,  # PR
    "Pt09": 4,  # CR
}

# Patient mappings.
#
# IMC h5ad obs.pt = "PatientN" with N ∈ {1,3,4,5,6,7,8,9} → direct map.
IMC_PT_TO_CANON = {f"Patient{i}": f"Pt{i:02d}" for i in (1, 3, 4, 5, 6, 7, 8, 9)}
#
# scRNA h5ad obs.patient uses *sequential* renumbering. Three trial
# patients (Pt03, Pt08, Pt09) were excluded for scRNA QC, so the h5ad's
# patient field re-numbers the remaining ones. Verified against the
# archive pipeline's run_macrophage_dge.py and against the paired
# pre/post structure: the 5 paired ({patient1, patient3, patient4,
# patient5, patient6}) are the 5 sequenced trial patients
# (Pt01, Pt04, Pt05, Pt06, Pt07).
SCRNA_PT_TO_CANON = {
    "patient1": "Pt01",
    "patient2": "Pt03",   # excluded for QC (still appears in h5ad)
    "patient3": "Pt04",
    "patient4": "Pt05",
    "patient5": "Pt06",
    "patient6": "Pt07",
    "patient7": "Pt08",   # excluded for QC
}
SCRNA_POST_TO_CANON = {f"{k}_post": v for k, v in SCRNA_PT_TO_CANON.items()}

# Patients in scRNA-seq cohort (paired pre + post). Pt03/Pt08/Pt09 excluded for QC.
SCRNA_COHORT = ["Pt01", "Pt04", "Pt05", "Pt06", "Pt07"]

# ---------------------------------------------------------------------------
# Visual styling — DOR-keyed palette (matches published 4-color scheme)
# ---------------------------------------------------------------------------
COLORS = {
    "DR_Baseline":   "#ADD8E6",   # lightblue
    "DR_PostNivo":   "#1E90FF",   # dodgerblue
    "P-NR_Baseline": "#FFB6C1",   # light coral
    "P-NR_PostNivo": "#B22222",   # firebrick
}
GROUP_ORDER = ["DR_Baseline", "DR_PostNivo", "P-NR_Baseline", "P-NR_PostNivo"]
LEGEND_LABELS = {"DR_Baseline":   "DR, pre",   "DR_PostNivo":   "DR, post",
                 "P-NR_Baseline": "P/NR, pre", "P-NR_PostNivo": "P/NR, post"}

DOR_COLOR_BY_RANK = {1: "#7f7f7f", 2: "#f1c40f", 3: "#3498db", 4: "#27ae60"}

# ---------------------------------------------------------------------------
# Stats helpers
# ---------------------------------------------------------------------------
def cohens_d(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    n1, n2 = len(a), len(b)
    s1, s2 = np.std(a, ddof=1), np.std(b, ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    return np.nan if sp < 1e-9 else (np.mean(a) - np.mean(b)) / sp


def welch(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    return stats.ttest_ind(a, b, equal_var=False)


def paired_t(pre_by_pt, post_by_pt):
    pts = sorted(set(pre_by_pt) & set(post_by_pt))
    if len(pts) < 2:
        return np.nan, np.nan, len(pts)
    res = stats.ttest_rel([pre_by_pt[p] for p in pts],
                            [post_by_pt[p] for p in pts])
    return res.statistic, res.pvalue, len(pts)


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    r = stats.spearmanr(x[m], y[m])
    return float(r.statistic), float(r.pvalue), int(m.sum())


def bh(p):
    p = np.asarray(p, float)
    valid = ~np.isnan(p)
    n = valid.sum()
    out = np.full(len(p), np.nan)
    if n == 0:
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    q_sorted = np.minimum.accumulate(
        (ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(q_sorted, 0, 1)
    out[np.where(valid)[0]] = q
    return out


def pformat(p, alpha=0.05):
    """Return (text, weight) for a p-value.
    text: 'p=0.012', 'p=0.51', or LaTeX-rendered scientific notation
          (e.g. 'p=1.2×10⁻⁵') for p < 0.001 (so we never display 'p=0.000').
    weight: 'bold' if p < alpha else 'normal'.
    When weight is bold and the text contains math, the math content is
    wrapped in \\mathbf{} so it renders bold (matplotlib mathtext does
    not honor fontweight inside $...$ blocks).
    Handles NaN as ('n.s.', 'normal') and p=0 as a lower-bound annotation.
    """
    if not np.isfinite(p):
        return ("n.s.", "normal")
    bold = p < alpha
    if p <= 0:
        body = r"p<10^{-300}"
        text = rf"$\mathbf{{{body}}}$" if bold else f"${body}$"
    elif p < 1e-3:
        e = int(np.floor(np.log10(p)))
        m = p / (10 ** e)
        body = rf"p={m:.1f}\times 10^{{{e}}}"
        text = rf"$\mathbf{{{body}}}$" if bold else f"${body}$"
    elif p < 0.1:
        text = f"p={p:.3f}"
    else:
        text = f"p={p:.2f}"
    weight = "bold" if bold else "normal"
    return (text, weight)


def stars(p):
    if not np.isfinite(p):
        return "n.s."
    if p < 0.0001: return "****"
    if p < 0.001:  return "***"
    if p < 0.01:   return "**"
    if p < 0.05:   return "*"
    return f"p = {p:.2g}"


def save_fig(fig, name, formats=("pdf", "svg")):
    """Save fig under FIG_DIR/{name}.{ext} for each format."""
    for fmt in formats:
        fig.savefig(FIG_DIR / f"{name}.{fmt}", dpi=300,
                    bbox_inches="tight", facecolor="white")
