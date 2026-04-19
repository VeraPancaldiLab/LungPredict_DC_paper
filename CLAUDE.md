# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bioinformatics research project analyzing dendritic cell (DC) subset spatial organization in Non-Small Cell Lung Cancer (NSCLC) Stage 4 and its impact on immunotherapy response. The pipeline covers 116 bulk RNAseq samples and 68 spatial transcriptomics samples, culminating in a 3-group patient stratification and validation in an independent cohort (TIPIT, 134 patients).

## Environment Setup

All scripts are in `scripts/` and use an isolated renv environment:

```r
install.packages('renv')
renv::activate()   # activate isolated environment
renv::restore()    # install locked package versions (R 4.3.1)
```

Open `scripts/LP_MOSAIC.Rproj` in RStudio to work within the project. R version 4.3.1, OS Ubuntu 22.04.

## Running Analyses

All analyses are R Markdown files — run sequentially in this order:

1. `scripts/Patient_names.Rmd` — harmonizes sample IDs across modalities
2. `scripts/Spatial_preprocessing.Rmd` — computes DC densities from raw spatial XLSX files → `input/Spatial_densities.RData`
3. `scripts/Cell_densities_analysis.Rmd` — DC density comparisons (R vs NR, biopsy site, survival)
4. `scripts/Megaclusters_analysis.Rmd` — megacluster spatial analysis, patient clustering → `output/Patients_groups.RData`
5. `scripts/RNAseq_analysis.Rmd` — TF activity (VIPER/ARACNE), pathway enrichment (PROGENy), immune deconvolution → `output/TF_matrix_*.csv`, `output/Pathway_matrix_*.csv`
6. `scripts/Signature_projection.Rmd` — validates signatures in independent cohort

## Code Architecture

### Key inputs (`input/`)
- `Clinical_data.csv` — patient metadata (treatment, response, survival)
- `Bulk_counts_raw.csv` — raw gene expression (116 samples)
- `Spatial_densities.RData` — DC subset densities per patient (68 samples)
- `MegaClusters_density.csv` / `MegaClusters_entropy.csv` — megacluster metrics
- `signatures/` — DC subset signatures from scRNA-seq for deconvolution

### Key outputs (`output/`)
- `Patients_groups.RData` — final patient stratification (3 groups, 66 patients); this is the central result used downstream
- `TF_matrix_*.csv` / `Pathway_matrix_*.csv` — activity scores per patient group
- `Databases/` — TF and pathway collections used for enrichment

### Shared utilities (`scripts/src/`)
- `utils.R` — `compute.anova.test()`, `extract_markers()`, `compute.survival.analysis()`, custom ggplot helpers
- `immunoscores.R` — immune scoring and cell type composition functions

### Data flow
Spatial XLSX → `Spatial_preprocessing` → densities → `Megaclusters_analysis` → patient groups (`Patients_groups.RData`) → `RNAseq_analysis` → TF/pathway matrices → `Signature_projection` (external validation)

## Coding Conventions

- Analyses produce figures saved to `Results/` and intermediate data to `output/`
- Patient group labels: `G1`, `G2`, `G3` (from hierarchical clustering on megacluster densities)
- Responder labels: `R` (responder), `NR` (non-responder)
- Tissue compartments: `Stroma`, `Tumor`, `Total`
- Biopsy sites: primary, node, metastasis
