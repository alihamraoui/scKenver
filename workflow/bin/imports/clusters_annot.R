#' @title Annotation extension to clusters
#' @description This function is supposed to be used after cells.annot.custom and when there is a column in metadata corresponding to cluster ID for each cell. It calculates percentage of each cell type in each cluster and attributes to all cells of the cluster the cell type corresponding to the highest proportion.
#' @param sobj a Seurat object (no default)
#' @param annotation_column CHARACTER : column name in @meta.data corresponding to the cell annotation (default to "cell_type")
#' @param cluster_column CHARACTER : column name in @meta.data corresponding to the cluster ID (no default)
#' @param new_column CHARACTER : a new column name for the cluster annotation extension (default to "cluster_type")
#' @param force LOGICAL : whether to erase or not the column 'new_column' if it already exists (default to FALSE, no forcing)
#' @return the input Seurat object with a new column indicating cluste-corrected cell type
#' @export
clusters_annot = function(sobj,
                          annotation_column = "cell_type",
                          cluster_column = NULL,
                          new_column = "cluster_type",
                          force = FALSE) {
  ## Check parameters
  if (is.null(sobj)) stop("sobj must be a Seurat object")
  
  metadata = sobj@meta.data
  
  if (new_column %in% colnames(sobj@meta.data) && !force) {
    stop(new_column, " should be a column that do not exist in sobj")
  } 
  
  if (!(annotation_column %in% colnames(sobj@meta.data))) {
    stop(annotation_column, " must be a column in metadata")
  } else {
    annotations = as.factor(sobj@meta.data[, annotation_column])
  }
  if (!(cluster_column %in% colnames(sobj@meta.data))) {
    stop(cluster_column, " must be a column in metadata")
  } else {
    clusters = as.factor(sobj@meta.data[, cluster_column])
  }
  
  ## Find cluster annotation
  tab = table(annotations, clusters)
  annotations_by_clusters = rownames(tab)[apply(tab, 2, which.max)]
  names(annotations_by_clusters) = colnames(tab)
  annotations_by_clusters = as.data.frame(annotations_by_clusters)
  annotations_by_clusters$clusters = rownames(annotations_by_clusters)
  colnames(annotations_by_clusters) = c(new_column, cluster_column)
  
  ## Add in metadata
  metadata = dplyr::left_join(metadata,
                              annotations_by_clusters,
                              by = cluster_column)
  metadata[, new_column] = as.factor(as.character(metadata[, new_column]))
  sobj@meta.data[, new_column] = metadata[, new_column]
  
  ## Output
  return(sobj)
}

