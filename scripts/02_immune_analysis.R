#!/usr/bin/env Rscript

# Load required libraries
library(Seurat)
library(SeuratWrappers)
library(SingleCellExperiment)

# Function to process immune cells
process_immune_cells <- function(seu) {
    # Normalize data
    seu <- NormalizeData(seu)
    
    # Find variable features
    seu <- FindVariableFeatures(seu)
    
    # Scale data
    seu <- ScaleData(seu)
    
    # Run PCA
    seu <- RunPCA(seu, features = VariableFeatures(seu))
    
    # Run UMAP
    seu <- RunUMAP(seu, 
                  dims = 1:30,
                  seed.use = 42)
    
    return(seu)
}

# Function to find marker genes
find_marker_genes <- function(seu, output_file) {
    # Find all markers
    markers <- FindAllMarkers(seu,
                            only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = 0.25)
    
    # Save markers
    write.csv(markers, output_file)
    
    return(markers)
}

# Main analysis
main <- function() {
    # Load integrated data
    big <- readRDS("data/integrated_seurat.rds")
    
    # Extract immune cells
    immune_cells <- subset(big, subset = type == "immune")
    
    # Process immune cells
    immune_cells <- process_immune_cells(immune_cells)
    
    # Find clusters
    immune_cells <- FindNeighbors(immune_cells, 
                                dims = 1:30)
    immune_cells <- FindClusters(immune_cells, 
                               resolution = 2)
    
    # Find marker genes
    markers <- find_marker_genes(immune_cells, 
                               "data/immune_markers.csv")
    
    # Generate visualizations
    pdf("data/immune_umap.pdf")
    print(DimPlot(immune_cells, 
                 reduction = "umap",
                 label = TRUE,
                 pt.size = 0.5))
    dev.off()
    
    # Generate feature plots for key markers
    key_markers <- c("CD79A", "CD8A", "KLRG1", "CD4")
    pdf("data/immune_markers.pdf")
    for(marker in key_markers) {
        print(FeaturePlot(immune_cells,
                         features = marker,
                         order = TRUE))
    }
    dev.off()
    
    # Save processed immune cells
    saveRDS(immune_cells, "data/immune_cells_processed.rds")
    
    # Convert to SingleCellExperiment and save
    sce_immune <- as.SingleCellExperiment(immune_cells)
    saveRDS(sce_immune, "data/immune_cells_sce.rds")
}

# Run main function
main() 