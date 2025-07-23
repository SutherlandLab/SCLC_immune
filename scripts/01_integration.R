#!/usr/bin/env Rscript

# Load required libraries
library(Seurat)
library(harmony)
library(SeuratWrappers)
library(SingleCellExperiment)

# Function to create and normalize Seurat object
create_seurat <- function(counts, project, min.cells = 3, min.features = 200) {
    seu <- CreateSeuratObject(counts = counts, 
                            project = project, 
                            min.cells = min.cells, 
                            min.features = min.features)
    seu <- NormalizeData(seu)
    return(seu)
}

# Function to remove unwanted genes
remove_unwanted_genes <- function(seu) {
    # Remove low count genes
    seu <- seu[rowSums(seu@assays$RNA@counts) > 10,]
    
    # Remove ribosomal and mitochondrial genes
    rps <- rownames(seu)[grep("^RPS", rownames(seu))]
    rpl <- rownames(seu)[grep("^RPL", rownames(seu))]
    mt <- rownames(seu)[grep("^MT-", rownames(seu))]
    
    # Remove lncRNAs
    lncrna_patterns <- c("^AC", "^AP", "^AL", "^AF")
    lncrna <- lapply(lncrna_patterns, function(pattern) {
        genes <- rownames(seu)[grep(pattern, rownames(seu))]
        genes[grepl("\\.", genes)]
    })
    lncrna <- unlist(lncrna)
    
    # Remove all unwanted genes
    seu <- seu[!rownames(seu) %in% c(rps, rpl, mt, lncrna),]
    
    return(seu)
}

# Main analysis
main <- function() {
    # Load data (replace with your data loading code)
    # seu_immune <- readRDS("path_to_immune_data.rds")
    # seu_tumor <- readRDS("path_to_tumor_data.rds")
    
    # Merge datasets
    big <- merge(seu_immune, y = seu_tumor, 
                add.cell.ids = c("immune", "tumor"), 
                project = "all")
    
    # Remove unwanted genes
    big <- remove_unwanted_genes(big)
    
    # Normalize data
    big <- NormalizeData(big)
    
    # Run PCA
    big <- FindVariableFeatures(big)
    big <- ScaleData(big)
    big <- RunPCA(big, features = VariableFeatures(big))
    
    # Run Harmony
    big <- RunHarmony(big, 
                     group.by.vars = "patient",
                     dims.use = 1:30,
                     plot_convergence = TRUE)
    
    # Run UMAP
    big <- RunUMAP(big, 
                  reduction = "harmony",
                  dims = 1:30,
                  seed.use = 42)
    
    # Find clusters
    big <- FindNeighbors(big, 
                        reduction = "harmony",
                        dims = 1:30)
    big <- FindClusters(big, 
                       resolution = 1.4)
    
    # Save results
    saveRDS(big, "data/integrated_seurat.rds")
    
    # Generate UMAP plot
    pdf("data/umap_integrated.pdf")
    print(DimPlot(big, 
                 reduction = "umap",
                 label = TRUE,
                 pt.size = 0.5))
    dev.off()
}

# Run main function
main() 