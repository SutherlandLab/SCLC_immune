 #!/usr/bin/env Rscript

# Load required libraries
library(Seurat)
library(SeuratWrappers)
library(SingleCellExperiment)

# Function to process biopsy data
process_biopsy_data <- function(seu) {
    # Normalize data
    seu <- NormalizeData(seu)
    
    # Find variable features
    seu <- FindVariableFeatures(seu)
    
    # Scale data
    seu <- ScaleData(seu)
    
    return(seu)
}

# Main analysis
main <- function() {
    # Load biopsy data (replace with your data loading code)
    # biopsy_0073 <- readRDS("path_to_22MH0073_data.rds")
    # biopsy_0081 <- readRDS("path_to_22MH0081_data.rds")
    
    # Process each biopsy dataset
    biopsy_0073 <- process_biopsy_data(biopsy_0073)
    biopsy_0081 <- process_biopsy_data(biopsy_0081)
    
    # Create list of objects for integration
    objects <- list(biopsy_0073, biopsy_0081)
    names(objects) <- c("22MH0073", "22MH0081")
    
    # Run fastMNN integration
    integrated <- RunFastMNN(object.list = objects,
                           features = 3000)
    
    # Run UMAP
    integrated <- RunUMAP(integrated,
                        reduction = "mnn",
                        dims = 1:20,
                        seed.use = 42)
    
    # Generate visualizations
    pdf("data/biopsy_umap.pdf")
    print(DimPlot(integrated,
                 reduction = "umap",
                 group.by = "orig.ident",
                 label = TRUE,
                 pt.size = 0.5))
    dev.off()
    
    # Save integrated data
    saveRDS(integrated, "data/biopsy_integrated.rds")
    
    # Convert to SingleCellExperiment and save
    sce_biopsy <- as.SingleCellExperiment(integrated)
    saveRDS(sce_biopsy, "data/biopsy_sce.rds")
}

# Run main function
main() 