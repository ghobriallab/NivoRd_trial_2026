# NivoRd Figure Scripts

Scripts to reproduce figures from the NivoRd manuscript (Letter to Blood).

## Requirements

```
python >= 3.10
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
```

## Input Data

Scripts expect input files in the working directory unless otherwise noted. Processed datasets are available at [Zenodo DOI to be assigned]. Raw sequencing reads are in dbGaP [accession to be assigned].

### scRNA-seq datasets (h5ad)

| File | Description |
|------|-------------|
| `CD138neg_immune_cells.h5ad` | CD138- bone marrow immune cells (90,849 cells) |
| `CD138pos_tumor_cells.h5ad` | CD138+ bone marrow plasma cells (23,362 cells) |

### scRNA-seq analysis results (CSV)

| File | Description |
|------|-------------|
| `posttreatment_macrophage_DGE_results.csv` | Post-treatment macrophage differential expression |
| `pretreatment_macrophage_DGE_results.csv` | Pre-treatment macrophage differential expression |
| `master_gsea_results.csv` | GSEA results across cell types and timepoints |
| `gene_expression_summary.csv` | Cytogenetic gene expression per case |
| `cytogenetic_validation.csv` | Cytogenetic validation statistics |
| `differential_plasma_macrophage_pretreatment.csv` | Plasma cell-macrophage LIANA results |
| `differential_cd4_macrophage_pretreatment.csv` | CD4-macrophage LIANA results |

### IMC datasets (CSV)

| File | Description |
|------|-------------|
| `phenotyped-obs.csv` | IMC cell phenotypes (52,043 cells, 7 lineages) |
| `pheno_8clus-cells_data.csv` | Cellular neighborhood clustering per cell |
| `pheno_8clus-sample_pct.csv` | CN proportions per sample |
| `combat_all_markers.csv` | Batch-corrected IMC marker expression |
| `nivo_clinical.csv` | Clinical metadata (patient, timepoint, response) |
| `nivo_MSPike.csv` | M-spike levels |

### IMC images (TIFF)

| File | Description |
|------|-------------|
| `Patient*_Baseline_*.tiff` | 34-channel IMC TIFFs (baseline) |
| `Patient*_PostNivo_*.tiff` | 34-channel IMC TIFFs (post-treatment) |

### External data (auto-downloaded)

| Source | Description |
|--------|-------------|
| GSE124310 (Zavidij et al.) | BM scRNA-seq across disease stages; fetched from GEO FTP |

## Scripts

### Main Figure 1

| Script | Panel | Description |
|--------|-------|-------------|
| `figure1B_cd8_abundance.py` | 1B | CD8+ T cell abundance by response group |
| `figure1C_cd4cd8_ratio.py` | 1C | CD4/CD8 ratio by response group |
| `figure1D_macrophage_cn2.py` | 1D | Macrophage abundance + CN2 proportion |
| `figure1E_imc_composites.py` | 1E | Representative 4-channel IMC composites |
| `figure1F_pdl1_expression.py` | 1F | PD-L1 expression (plasma cells + macrophages) |

### Main Figure 2

| Script | Panel | Description |
|--------|-------|-------------|
| `figure2A_immune_umap.py` | 2A | CD138- immune UMAP (12 cell types) |
| `figure2B_macrophage_volcano.py` | 2B | Post-treatment macrophage DGE volcano |
| `figure2C_macrophage_gsea.py` | 2C | GSEA hallmark pathways |
| `figure2D_cd4cd8_vs_cn.py` | 2D | CD4/CD8 ratio vs CN2 and CN4 correlation |
| `figure2EF_zavidij_macrophage_scoring.py` | 2E-F | Zavidij external cohort macrophage scoring |

### Supplementary Figures

| Script | Panel | Description |
|--------|-------|-------------|
| `figureS1F_consort.py` | S1F | CONSORT flow diagram |
| `figureS2A_marker_heatmap.py` | S2A | Immune cell type marker heatmap |
| `figureS2B_macrophage_umap.py` | S2B | Macrophage-highlighted UMAP |
| `figureS2C_imc_marker_dotplot.py` | S2C | IMC lineage marker expression |
| `figureS2D_cn_composition.py` | S2D | Cellular neighborhood composition |
| `figureS2E_tumor_umap.py` | S2E | CD138+ plasma cell UMAP by patient |
| `figureS2F_cytogenetic_heatmap.py` | S2F | Cytogenetic gene validation heatmap |
| `figureS3A_cn_proportions.py` | S3A | CN proportions (CN3, CN4, CN7) |
| `figureS3B_plasma_cell_burden.py` | S3B | Plasma cell burden + M-spike |
| `figureS3C_lineage_abundances.py` | S3C | B cell, myeloid, megakaryocyte abundances |
| `figureS4A_pretreatment_volcano.py` | S4A | Pre-treatment macrophage DGE volcano |
| `figureS4BC_cell_interactions.py` | S4B-C | LIANA ligand-receptor analysis |

### Figures without scripts

| Panel | Content | Reason |
|-------|---------|--------|
| 1A | Trial schematic | Design diagram |
| S1A-E | Swimmer plot, Kaplan-Meier curves | Clinical team (Prism/R) |

## Usage

```bash
# scRNA-seq figures (need h5ad files in working directory)
python figure2A_immune_umap.py

# IMC figures (need CSV files in working directory)
python figure1B_cd8_abundance.py

# IMC composites (need TIFF files)
python figure1E_imc_composites.py --input path/to/tiffs --output path/to/out

# External cohort (auto-downloads from GEO)
python figure2EF_zavidij_macrophage_scoring.py
```

Each script outputs PDF, PNG, and SVG at 300 DPI.
