
fx_1NN = function(i,location_in){
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  return(min(line_i))
}


fx_CHAOS = function(clusterlabel, location){
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }
  matched_location = scale(matched_location)
  dist_val = rep(0,length(unique(clusterlabel)))
  count = 0
  for(k in unique(clusterlabel)){
    count = count + 1
    location_cluster = matched_location[which(clusterlabel == k),]
    if(length(location_cluster)==2){next}
    results = mclapply(1:dim(location_cluster)[1], fx_1NN, location_in=location_cluster, mc.cores = 5)
    dist_val[count] = sum(unlist(results))
  }
  dist_val = na.omit(dist_val)
  return(sum(dist_val)/length(clusterlabel))

}


cluster_slots <- function(obj){
      for(assays in names(slot(obj, "assays"))){
        
          DefaultAssay(obj) <- assays
          
          obj <- NormalizeData(object = obj, assay = assays, verbose = FALSE)
          obj = Seurat::ScaleData(obj,features = rownames(obj), assay = assays, verbose = F)
          obj = Seurat::FindVariableFeatures(obj,
                                          assay = assays,
                                          normalization.method = "vst")
          
          obj <- RunPCA(obj, assay = assays, reduction.name= paste0(assays, '_PCA'),  verbose = FALSE)
          obj <- FindNeighbors(obj, assay = assays, 
                               reduction = paste0(assays, "_PCA"), 
                               dims = 1:15,  verbose = FALSE)
          
          obj <- FindClusters(obj, assay = assays, verbose = FALSE)
          obj <- RunUMAP(obj, assay = assays, 
                         reduction = paste0(assays, "_PCA"), 
                         reduction.name = paste0(assays, '_UMAP'), 
                         dims = 1:15,  
                         verbose = FALSE)
    }
    return(obj)
}


calc_LISI <- function(obj, reduction, ident="orig.ident"){
    meta_data <- data.frame(colnames(obj), ident=obj@meta.data$seurat_clusters)
    cord.umap <- Embeddings(obj, reduction = reduction)[, 1:2]
    lisi <- lisi::compute_lisi(cord.umap, meta_data, c("ident"))
    meanLISI <- lisi$ident
}


boxplot_lisi <- function(df, title){
    ggplot(df, aes(x = forcats::fct_reorder(workflow, LISI, .fun = median, .desc = T), y = LISI, fill = workflow)) +
              geom_boxplot(width = 0.5, outlier.shape = NA) +
              theme_classic()+
              labs(x = "", y = "LISI", title=title) +
              theme(axis.text.x =  element_text(colour = "black"),
                    axis.text = element_text(face="bold"),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank()) 
}

bar_metric <- function(df, title, desc=F){
    ggplot(df, aes(x = forcats::fct_reorder(workflow, metric, .fun = sum, .desc = desc), y = metric, fill = workflow)) +
              geom_bar(stat = "identity") +
              theme_classic()+
              labs(x = "", y = "", title=title) +
              theme(axis.text.x =  element_text(colour = "black"),
                   axis.text = element_text(face="bold"),
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank())
}