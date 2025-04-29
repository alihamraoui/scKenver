mix <- function(tool.obj, illumina.obj, toolName){
  integ.list <- list()
  integ.list[[toolName]] <- tool.obj
  integ.list[["CellRanger"]] <- illumina.obj
  
  for (i in 1:length(integ.list)) {
    integ.list[[i]] <- NormalizeData(integ.list[[i]], verbose = F)
    integ.list[[i]] <- FindVariableFeatures(integ.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
  }
  anchors  <- FindIntegrationAnchors(object.list = integ.list, dims = 1:30, verbose = F)
  integrated_obj  <- IntegrateData(anchorset = anchors, dims = 1:30, verbose = F)
}


integrate <- function(tool.obj, illumina.obj, toolName){
  obj <- merge(tool.obj, illumina.obj)
  obj <- NormalizeData(obj, verbose = F)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = F)
  obj <- ScaleData(obj, verbose = F)
  obj <- RunPCA(obj, npcs =15, verbose = F)
  obj <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                         verbose = F)
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- RunUMAP(obj, dims = 1:15, reduction = "integrated.cca", verbose = F)
}


UMAP_projection <- function(integrated_obj){
  DefaultAssay(integrated_obj) <- "RNA"
  integrated_obj <- NormalizeData(integrated_obj, verbose = F)
  integrated_obj <- FindVariableFeatures(integrated_obj, selection.method = "vst", nfeatures = 2000, verbose = F)
  integrated_obj <- ScaleData(integrated_obj, verbose = F)
  integrated_obj <- RunPCA(integrated_obj, npcs = 5, verbose = F)
  integrated_obj <- RunUMAP(integrated_obj, reduction = "pca", dims = 1:15, verbose = F)
}


calc_meanLISI <- function(obj, reduction, ident="orig.ident"){
  if (ident == "cell_type"){
    meta_data <- data.frame(colnames(obj), ident=obj@meta.data$cell_type)
  }else{
    meta_data <- data.frame(colnames(obj), ident=obj@meta.data$orig.ident)
    colnames(meta_data) <- c("CB", "ident")
    meta_data <- find_common_barcodes(meta_data)
    obj <- subset(obj, cells = meta_data$CB)
  }
  cord.umap <- Embeddings(obj, reduction = reduction)[, 1:2]
  lisi <- lisi::compute_lisi(cord.umap, meta_data, c("ident"))
  medianLISI <- median(lisi$ident)
  
}


calc_LISI <- function(obj, reduction, ident="orig.ident"){
  if (ident == "cell_type"){
    meta_data <- data.frame(colnames(obj), ident=obj@meta.data$cell_type)
  }else{
    meta_data <- data.frame(colnames(obj), ident=obj@meta.data$orig.ident)
    colnames(meta_data) <- c("CB", "ident")
    meta_data <- find_common_barcodes(meta_data)
    obj <- subset(obj, cells = meta_data$CB)
  }
  cord.umap <- Embeddings(obj, reduction = reduction)[, 1:2]
  lisi <- lisi::compute_lisi(cord.umap, meta_data, c("ident"))
  meanLISI <- lisi$ident
}


find_common_barcodes <- function(df){
  
  df$barcode_seq <- sub("_.$", "", df$CB)
  
  tool1 <- names(table(df$ident))[1]
  tool2 <- names(table(df$ident))[2]
  
  tool1_df <- subset(df, ident == tool1)
  tool2_df <- subset(df, ident == tool2)
  
  common_barcodes <- intersect(tool1_df$barcode_seq, tool2_df$barcode_seq)
  
  common_rows <- df[df$barcode_seq %in% common_barcodes, ]
  
  df <- df[df$barcode_seq %in% common_barcodes, ][,-which(colnames(df) == "barcode_seq")]
  return(df)
}


clustering <- function (matrix, min.features, nfeatures, dim) {
  ScObject = CreateSeuratObject(counts = matrix, min.cells = 3, min.features = min.features)
  ScObject = NormalizeData(ScObject, normalization.method = "LogNormalize", verbose = F)
  ScObject = FindVariableFeatures(ScObject, selection.method = "vst", nfeatures = nfeatures, verbose = F)
  all.genes <- rownames(ScObject)
  ScObject <- ScaleData(ScObject, features = all.genes, verbose = F)
  ScObject <- RunPCA(ScObject, features = VariableFeatures(object = ScObject), verbose = F)
  ScObject <- FindNeighbors(ScObject, dims = 1:dim, verbose = F)
  ScObject <- FindClusters(ScObject, resolution = 0.5, verbose = F)
  ScObject <- RunUMAP(ScObject, dims = 1:dim, verbose = F)
}


annotate.aud <- function (ScObject) {
  if ("cell_type" %in% colnames(ScObject@meta.data)) {
    ScObject$cell_type = NULL
  }
  ScObject = cell_annot_custom(ScObject,
                               assay = "RNA",          
                               newname = "cell_type",  
                               markers = cell_markers) 
}


box_lisi <- function(vec){
  lisi_data <- data.frame(
    group = rep('iLISI', length(vec)),
    value = vec)
  
  bp <- ggplot(lisi_data, aes(x = group, y = value)) +
    geom_boxplot(width = 1.15, outlier.shape = NA) +
    labs(
      x = "",
      y = ""
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    ) + coord_flip()
  return(inset_element(bp, left = 0.005, bottom = 0.95, right = 1, top = 1.15, on_top = F))
}

DimPlot.ctype<- function (Sobj, npipe) { DimPlot(Sobj, group.by = "cell_type", cols = color_markers) + 
    ggtitle(paste(npipe, "| median cLISI", round(calc_meanLISI(Sobj, "umap", "cell_type"),2))) +
    nn_leg +
    theme(aspect.ratio = 1)}

DimPlot.integ <- function(Sobj) { DimPlot(Sobj, reduction = reduction,  group.by = c("orig.ident"), cols = c('gray57', '#d175b8')) + 
    ggtitle(paste('mean iLISI - 1 = ', round(calc_meanLISI(Sobj, reduction),2)-1)) + 
    theme(aspect.ratio = 1)}