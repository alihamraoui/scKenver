library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

clustering <- function (matrix, min.features, nfeatures, dim) {
  ScObject = CreateSeuratObject(counts = matrix, min.cells = 3, min.features = min.features)
  ScObject = NormalizeData(ScObject, normalization.method = "LogNormalize", verbose = F)
  ScObject = FindVariableFeatures(ScObject, selection.method = "vst", nfeatures = nfeatures, verbose = F)
  all.genes <- rownames(ScObject)
  ScObject <- ScaleData(ScObject, features = all.genes, verbose = F)
  ScObject <- RunPCA(ScObject, features = VariableFeatures(object = ScObject), verbose = F)
  ScObject <- FindNeighbors(ScObject, dims = 1:dim, verbose = F)
  ScObject <- FindClusters(ScObject, resolution = 0.5, verbose = F)
  ScObject <- RunTSNE(ScObject, dims = 1:dim, verbose = F)
}

dim.plot <- function (Sobj, npipe) { DimPlot(Sobj, reduction = 'tsne') + 
    ggtitle(paste(npipe)) +
    theme(aspect.ratio = 1)}

compare_pseudobulk_normalized <- function(seurat_obj1, seurat_obj2, celltype_column = "CellType", method = "pearson", normalize = TRUE, plot = TRUE) {
  
  # -------- Helper function: generate pseudo-bulk --------
  generate_pseudobulk <- function(seurat_obj) {
    Idents(seurat_obj) <- seurat_obj[[celltype_column]][,1]
    counts <- GetAssayData(seurat_obj, slot = "counts")
    cell_types <- Idents(seurat_obj)
    
    valid_cells <- !is.na(cell_types)
    cell_types <- cell_types[valid_cells]
    
    pseudo_bulk_list <- list()
    for (ct in unique(cell_types)) {
      cells_in_type <- WhichCells(seurat_obj, idents = ct)
      counts_sum <- rowSums(counts[, cells_in_type, drop = FALSE])
      pseudo_bulk_list[[ct]] <- counts_sum
    }
    pseudo_bulk_matrix <- do.call(cbind, pseudo_bulk_list)
    return(pseudo_bulk_matrix)
  }
  
  # -------- Step 1: Generate pseudo-bulk for both objects --------
  message("Generating pseudo-bulk matrices...")
  pb1 <- generate_pseudobulk(seurat_obj1)
  pb2 <- generate_pseudobulk(seurat_obj2)
  # -------- Step 2: Filter common genes and cell types --------
  common_genes <- intersect(rownames(pb1), rownames(pb2))
  common_celltypes <- intersect(colnames(pb1), colnames(pb2))
  
  
  if (length(common_celltypes) == 0) {
    stop("No common cell types found between the two objects!")
  }
  
  pb1 <- pb1[common_genes, common_celltypes, drop = FALSE]
  pb2 <- pb2[common_genes, common_celltypes, drop = FALSE]
  
  # -------- Step 3: Normalization --------
  if (normalize) {
    message("Normalizing with CPM + log1p...")
    
    normalize_cpm_log <- function(mat) {
      lib_sizes <- colSums(mat)
      cpm <- sweep(mat, 2, lib_sizes, FUN = "/") * 1e6
      log_cpm <- log1p(cpm)
      return(log_cpm)
    }
    
    pb1 <- normalize_cpm_log(pb1)
    pb2 <- normalize_cpm_log(pb2)
  }
  
  # -------- Step 4: Compute correlations --------
  cor_results <- data.frame(
    CellType = common_celltypes,
    Correlation = NA
  )
  
  for (ct in common_celltypes) {
    counts1 <- pb1[, ct]
    counts2 <- pb2[, ct]
    
    cor_val <- cor(counts1, counts2, method = method)
    cor_results$Correlation[cor_results$CellType == ct] <- cor_val
  }
  
  # -------- Step 5: Optional plot --------
  if (plot) {
    ggplot(cor_results, aes(x = CellType, y = Correlation)) +
      geom_col(fill = "steelblue") +
      ylim(0, 1) +
      theme_minimal() +
      labs(title = paste("Correlation of Normalized Pseudo-bulk Counts (", method, ")", sep = ""), 
           y = "Correlation", x = "Cell Type") -> p
    print(p)
  }
  
  return(cor_results)
}

