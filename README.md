# NivoRd Figure Scripts

Scripts to reproduce figures from the NivoRd manuscript (Letter to Blood).

## Layout

```
2-Scripts/
├── IMC/                                     # Imaging mass cytometry pipeline
│   ├── IMC_analysis.ipynb                   # Preprocessing, clustering, CN, S2C dotplots
│   ├── IMC_figures.Rmd                      # Downstream plots (Fig 1B/1C/1D/1F/2D, S2D, S3A/B/C)
│   └── figure1E_imc_composites.py           # 4-channel TIFF composites (Fig 1E)
├── scRNA-seq/                               # Single-cell RNA-seq scripts (Python)
│   ├── figure2A_immune_umap.py              # CD138⁻ immune UMAP
│   ├── figure2B_macrophage_volcano.py       # Post-treatment macrophage DGE volcano
│   ├── figure2C_macrophage_gsea.py          # Hallmark GSEA bars
│   ├── figure2EF_zavidij_macrophage_scoring.py  # External Zavidij cohort (Fig 2E-F)
│   ├── figureS2A_marker_heatmap.py          # Immune cell-type marker heatmap
│   ├── figureS2B_macrophage_umap.py         # Macrophage-highlighted UMAP
│   ├── figureS2E_tumor_umap.py              # CD138⁺ tumor UMAP by patient
│   ├── figureS2F_cytogenetic_heatmap.py     # Cytogenetic validation heatmap
│   ├── figureS4A_pretreatment_volcano.py    # Pre-treatment macrophage volcano
│   └── figureS4BC_cell_interactions.py      # LIANA ligand-receptor axes
└── clinical/
    └── figureS1F_consort.py                 # CONSORT flow diagram
```

## Figure → source mapping

### Main Figure 1

| Panel | Content | Script |
|-------|---------|--------|
| 1A | Trial schematic | n/a (design diagram) |
| 1B | CD8⁺ T cell abundance | `IMC/IMC_figures.Rmd` — "IMC cell %" (T8) |
| 1C | CD4/CD8 ratio | `IMC/IMC_figures.Rmd` — "IMC CD4/8 ratio" |
| 1D | Macrophage abundance + CN2 | `IMC/IMC_figures.Rmd` — "IMC cell %" (Mac) + "IMC CN %" (CN2) |
| 1E | Representative IMC composites | `IMC/figure1E_imc_composites.py` |
| 1F | PD-L1 (plasma cells + macrophages) | `IMC/IMC_figures.Rmd` — "PD-L1 expression" |

### Main Figure 2

| Panel | Content | Script |
|-------|---------|--------|
| 2A | CD138⁻ immune UMAP | `scRNA-seq/figure2A_immune_umap.py` |
| 2B | Post-treatment macrophage volcano | `scRNA-seq/figure2B_macrophage_volcano.py` |
| 2C | GSEA Hallmark pathways | `scRNA-seq/figure2C_macrophage_gsea.py` |
| 2D | CD4/CD8 vs CN2 and CN4 correlations | `IMC/IMC_figures.Rmd` — "IMC CD4/8 ratio and CN %" |
| 2E-F | Zavidij external cohort | `scRNA-seq/figure2EF_zavidij_macrophage_scoring.py` |

### Supplementary Figures

| Panel | Content | Script |
|-------|---------|--------|
| S1A-E | Swimmer + Kaplan-Meier curves | n/a (clinical team; Prism/R) |
| S1F | CONSORT flow diagram | `clinical/figureS1F_consort.py` |
| S2A | Immune cell marker heatmap | `scRNA-seq/figureS2A_marker_heatmap.py` |
| S2B | Macrophage-highlighted UMAP | `scRNA-seq/figureS2B_macrophage_umap.py` |
| S2C | IMC lineage marker expression | `IMC/IMC_analysis.ipynb` — "Publication plots" |
| S2D | CN composition stacked bar | `IMC/IMC_figures.Rmd` — "IMC CN overview" |
| S2E | CD138⁺ tumor UMAP | `scRNA-seq/figureS2E_tumor_umap.py` |
| S2F | Cytogenetic validation | `scRNA-seq/figureS2F_cytogenetic_heatmap.py` |
| S3A | CN proportions (CN3/CN4/CN7) | `IMC/IMC_figures.Rmd` — "IMC CN %" |
| S3B | Plasma cell burden + M-spike | `IMC/IMC_figures.Rmd` — "Clinical tumour burden metrics" |
| S3C | Lineage abundances (B, myeloid, MK) | `IMC/IMC_figures.Rmd` — "IMC cell %" (loop) |
| S4A | Pre-treatment macrophage volcano | `scRNA-seq/figureS4A_pretreatment_volcano.py` |
| S4B-C | LIANA ligand-receptor axes | `scRNA-seq/figureS4BC_cell_interactions.py` |

## Pipeline overview

```
IMC raw 34-channel TIFFs ─► (image import) ─► IMC-import.h5ad
IMC-import.h5ad ─► IMC_analysis.ipynb ─► phenotyped-obs.csv + CN CSVs + S2C dotplots
phenotyped-obs.csv + clinical CSVs ─► IMC_figures.Rmd ─► Fig 1B/1C/1D/1F/2D + Supp boxplots

scRNA-seq h5ad ─► scRNA-seq/*.py ─► Fig 2A/2B/2C + Supp Figs
Zavidij GEO (auto-downloaded) ─► figure2EF_zavidij_macrophage_scoring.py ─► Fig 2E/2F
```

## Data inputs

All processed data are deposited at **Zenodo [DOI to be assigned]**. Raw sequencing reads are at dbGaP [accession to be assigned]. Raw IMC images (34-channel TIFFs) are available from the corresponding authors subject to informed consent.

### Data files on Zenodo

| File | Size | Used by |
|------|------|---------|
| `CD138neg_immune_cells.h5ad` | 2.6 GB | scRNA-seq immune scripts |
| `CD138pos_tumor_cells.h5ad` | 652 MB | scRNA-seq tumor scripts |
| `IMC-import.h5ad` | 11 MB | `IMC_analysis.ipynb` (starting AnnData after image import) |
| `posttreatment_macrophage_DGE_results.csv` | — | `figure2B_macrophage_volcano.py` |
| `pretreatment_macrophage_DGE_results.csv` | — | `figureS4A_pretreatment_volcano.py` |
| `master_gsea_results.csv` | — | `figure2C_macrophage_gsea.py` |
| `gene_expression_summary.csv` | — | `figureS2F_cytogenetic_heatmap.py` |
| `cytogenetic_validation.csv` | — | `figureS2F_cytogenetic_heatmap.py` |
| `differential_plasma_macrophage_pretreatment.csv` | — | `figureS4BC_cell_interactions.py` |
| `differential_cd4_macrophage_pretreatment.csv` | — | `figureS4BC_cell_interactions.py` |
| `phenotyped-obs.csv` | — | `IMC_figures.Rmd` (most IMC plots) |
| `pheno_8clus-cells_data.csv` | — | `IMC_figures.Rmd` + `IMC_analysis.ipynb` |
| `pheno_8clus-sample_pct.csv` | — | `IMC_figures.Rmd` |
| `combat_all_markers.csv` | — | `IMC_figures.Rmd` (PD-L1 section) |
| `nivo_clinical.csv` | — | `IMC_figures.Rmd` (clinical burden section) |
| `nivo_md.csv` | — | `IMC_figures.Rmd` |
| `nivo_MSPike.csv` | — | `IMC_figures.Rmd` |

## Configuring data paths

The IMC notebook and R markdown expect processed data under a single root directory, `DATA_DIR`, with the layout:

```
$DATA_DIR/
├── nivo_clinical.csv
├── nivo_md.csv
├── nivo_MSPike.csv
├── outs/
│   ├── phenotyped-obs.csv
│   ├── combat_all_markers.csv
│   └── import.h5ad          # output of IMC_analysis.ipynb (equivalent to Zenodo IMC-import.h5ad upstream of phenotyping)
├── CN_data/
│   ├── pheno_8clus-cells_data.csv
│   └── pheno_8clus-sample_pct.csv
└── figures/                 # created automatically; plot outputs go here
```

Set the root via the environment variable `NIVO_IMC_DATA_DIR`, e.g.:

```bash
export NIVO_IMC_DATA_DIR=/path/to/your/nivo_imc_data/
```

Both `IMC_figures.Rmd` (R) and `IMC_analysis.ipynb` (Python) read this variable and fall back to `/mnt/disks/data/imc/Nivo/` (the original processing environment) if unset.

## Requirements

### Python (3.10+)

```
anndata
scanpy
matplotlib
numpy
pandas
scipy
adjustText
gseapy
tifffile
Pillow
inmoose            # IMC ComBat normalization
```

### R (4.2+) — for `IMC_figures.Rmd`

```
tidyverse
ggpubr
effsize
```

## Usage

```bash
# scRNA-seq figures (need h5ad files)
python scRNA-seq/figure2A_immune_umap.py

# IMC figures (after running IMC_analysis.ipynb to generate CSVs)
export NIVO_IMC_DATA_DIR=/path/to/nivo_imc_data/
Rscript -e "rmarkdown::render('IMC/IMC_figures.Rmd')"

# IMC composites (need raw 34-channel TIFFs)
python IMC/figure1E_imc_composites.py --input path/to/tiffs --output path/to/out

# External cohort (auto-downloads from GEO)
python scRNA-seq/figure2EF_zavidij_macrophage_scoring.py
```

Each script outputs PDF, PNG, and SVG at 300 DPI unless otherwise noted.
