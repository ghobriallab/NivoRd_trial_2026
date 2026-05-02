# NivoRd — figure code

**Manuscript:** *Macrophage reprogramming, not CD8 T cell dynamics, characterizes durable disease control following PD-1 blockade in smoldering myeloma*

**Status:** Letter to *Leukemia* (under review)

**Trial:** NivoRd, single-arm phase II (NCT02903381; Dana-Farber IRB 16-242).
n=8 patients with high-risk smoldering multiple myeloma treated with
nivolumab + lenalidomide + dexamethasone.

**Data:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19430235.svg)](https://doi.org/10.5281/zenodo.19430235)
The associated processed data deposit on Zenodo is referenced via concept
DOI `10.5281/zenodo.19430235`; the v4-DOR snapshot is at record
`19988651` (DOI `10.5281/zenodo.19988651`). See also the deposit's
own README and DATA_DICTIONARY.

## Setup

```bash
git clone <this-repo>
cd 2-Scripts
pip install -r requirements.txt
python download_data.py        # fetches data files from Zenodo into data/
python Fig1B_CD8_abundance.py  # any figure script
```

By default `download_data.py` writes into `./data/`. Override with
`NIVO_DATA_DIR=/some/other/path python download_data.py`. The same
variable is read by `_common.py` so figure scripts will pick up the
custom path.

The scRNA-seq h5ads (`CD138neg_immune_cells.h5ad`,
`CD138pos_tumor_cells.h5ad`) are part of the v4-DOR deposit
(`19988651`), filtered to NivoRd trial cells only (Source ∈
{`Nivo`, `Nivo_post`}; external NBM / ELO reference cohorts removed
prior to deposit and available from their original published sources).
`download_data.py` fetches them by default along with the IMC and CSV
files; the full deposit is ~230 MB.

## Figure-to-script mapping

### Main Figure 1 — IMC

| Panel  | Script                               | Notes                                   |
| ------ | ------------------------------------ | --------------------------------------- |
| 1A     | `Fig1A_placeholder.py`               | Placeholder — final schematic in BioRender |
| 1B     | `Fig1B_CD8_abundance.py`             | CD8 lineage % (IMC)                     |
| 1C     | `Fig1C_Mac_abundance.py`             | Macrophage lineage % (IMC)              |
| 1D     | `Fig1D_CN7.py`                       | CN7 (Mac-dominant neighborhood)         |
| 1E     | `Fig1E_dLogMacCD4.py`                | log Mac/CD4 dynamics                    |
| 1F     | `Fig1F_placeholder.py`               | Placeholder — IMC composite TIFFs       |
| 1G     | `Fig1G_PCPDL1.py`                    | PD-L1+ plasma cells                     |
| 1H     | `Fig1H_MacPDL1.py`                   | PD-L1+ macrophages                      |
| 1I     | `Fig1I_PDL1_coordinated.py`          | Coordinated PD-L1 increase              |

### Main Figure 2 — scRNA-seq + external validation

| Panel  | Script                                  | Notes                                   |
| ------ | --------------------------------------- | --------------------------------------- |
| 2A     | `Fig2A_immune_umap.py`                  | Needs scRNA h5ads                        |
| 2B     | `Fig2B_pretreatment_volcano.py`         | Macrophage DGE pretreatment              |
| 2C     | `Fig2C_posttreatment_volcano.py`        | Macrophage DGE posttreatment             |
| 2D     | `Fig2D_pretreatment_GSEA.py`            | Macrophage GSEA pretreatment             |
| 2E     | `Fig2E_posttreatment_GSEA.py`           | Macrophage GSEA posttreatment            |
| 2F     | `Fig2F_CD86_ordinal.py`                 | CD86 vs DOR_ordinal                      |
| 2G/2H  | `Fig2GH_Zavidij.py`                     | External validation (Zavidij et al.)     |

### Supplementary

| Panel    | Script                              | Notes                                     |
| -------- | ----------------------------------- | ----------------------------------------- |
| S1A–F    | `SuppFigS1_placeholder.py`          | Placeholder — clinical CONSORT/swimmer    |
| S2A      | `SuppFigS2A_marker_heatmap.py`      | Needs scRNA h5ads                         |
| S2B      | `SuppFigS2B_macrophage_umap.py`     | Needs scRNA h5ads                         |
| S2C      | `SuppFigS2C_imc_marker_heatmap.py`  | IMC marker × cluster heatmap               |
| S2D      | `SuppFigS2D_CN_composition.py`      | CN composition by lineage                 |
| S2E      | `SuppFigS2E_tumor_umap.py`          | Needs scRNA h5ads                         |
| S2F      | `SuppFigS2F_cytogenetic_heatmap.py` | Cytogenetic vs expression validation      |
| S3A      | `SuppFigS3A_lineage_abundances.py`  | Lineage abundance summary                 |
| S3B      | `SuppFigS3B_CN_abundances.py`       | CN abundance summary                      |
| S3C      | `SuppFigS3C_ordinal_correlates.py`  | DOR_ordinal correlations                   |
| S4A/B    | `SuppFigS4_cell_interactions.py`    | LIANA differential interactions           |

### Helper / analysis scripts

- `_common.py` — shared classifiers (DOR/COT/FAST), palettes, stats helpers
- `_data_loaders.py` — h5ad loaders with canonical patient/timepoint mapping
- `_boxplot.py`, `_volcano_helper.py`, `_gsea_helper.py` — plotting helpers
- `analyze_pdl1_interaction.py` — PD-L1 supplementary analysis
- `run_dge_lmm.py` — macrophage DGE driver (pre and post)

## Citation

> Perron N. et al. *Macrophage reprogramming, not CD8 T cell dynamics,
> characterizes durable disease control following PD-1 blockade in
> smoldering myeloma.* (Manuscript under review, 2026.)

## Contact

Noé Perron — noeperron01@gmail.com — Ghobrial Lab, Dana-Farber Cancer Institute
