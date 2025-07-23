# Single-cell RNA-seq Analysis Pipeline

This repository contains scripts for analyzing single-cell RNA-seq data from tumor and immune cells.

## Analysis Workflow

### 1. Initial Integration and Cell Type Separation
- Script: `01_integration.R`
- Performs batch correction on tumor and immune cells using Harmony
- Uses top 30 PCs for dimension reduction
- Separates tumor and immune cells based on clustering results

### 2. Immune Cell Analysis
- Script: `02_immune_analysis.R`
- Re-analyzes immune cells after separation
- Performs PCA and integration
- Constructs KNN graph using top 30 PCs
- Identifies marker genes and generates UMAP visualization

### 3. Biopsy Data Integration
- Script: `03_biopsy_integration.R`
- Integrates normalized data from biopsies (22MH0073 and 22MH0081)
- Uses fastMNN for integration
- Generates UMAP visualization using top 20 MNN dimensions

## Script Usage

1. First, ensure all required R packages are installed:
```R
install.packages(c("Seurat", "harmony", "SeuratWrappers"))
```

2. Run the scripts in order:
```bash
Rscript scripts/01_integration.R
Rscript scripts/02_immune_analysis.R
Rscript scripts/03_biopsy_integration.R
```

## Output Files

- Integrated Seurat objects
- UMAP visualizations
- Marker gene lists
- Clustering results

## Dependencies

- R >= 4.0.0
- Seurat >= 4.0.0
- harmony
- SeuratWrappers
- SingleCellExperiment 