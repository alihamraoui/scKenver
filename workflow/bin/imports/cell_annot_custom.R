#' @title Cell type annotation with Seurat::AddModuleScore
#' @description This function performs automatic annotation of clusters in a Seurat object, based on a provided list of markers and cellular populations names. This function uses Seurat::AddModuleScore to attribute a score to each cell. Then, it gives to the cell to cell type corresponding to the highest score.
#' @param sobj A Seurat object (no default)
#' @param assay CHARACTER : on which assay from sobj to perform cell type annotation ? (default to 'RNA')
#' @param newname CHARACTER : name for the new column in metadata (default to 'newgroup')
#' @param markers LIST : a named list containing markers for each population (no default)
#' @param add_score LOGICAL : whether to add the columns containing scores to metadata in sobj (default to FALSE)
#' @param prefix CHARACTER : if add_score is set to TRUE, the prefix for new column names in metadata. For example, if a population type is called "pop1", which names to give to the new column in metadata containing score ? If prefix is "score_", then new column will be "score_pop1" (default to "score_")
#' @param seed INTEGER : the seed to be used by Seurat::AddModuleScore (default to 1337L)
#' @return This function returns the input Seurat object with a new column in metadata, containing factor levels from names(markers)
#' @importFrom Seurat AddModuleScore
#' @export
cell_annot_custom = function(sobj,
                             assay = "RNA",
                             newname = 'newgroup',
                             markers = markers,
                             add_score = FALSE,
                             prefix = "score_",
                             seed = 1337L) {
  if (is.null(markers)) {stop("markers must be a named list")}
  if (is.null(sobj)) {stop("sobj must be a Seurat object")}
  if (newname %in% colnames(sobj@meta.data)) {warning(newname, " is already a column in metadata and will be over-written.")}

  ## Attribute scores according to markers
  scores_list = lapply(markers, FUN = function(population) {
    suppressWarnings(Seurat::AddModuleScore(sobj,
                                            features = list(population),
                                            search = FALSE,
                                            seed = seed,
                                            assay = assay,
                                            name = "SCORE"))@meta.data$SCORE1})

  ## Organize score list in dataframe
  scores_DF = base::do.call(base::cbind.data.frame, scores_list) # colnames are names(markers)
  scores_names = paste(prefix, colnames(scores_DF), sep = "")
  colnames(scores_DF) = scores_names

  if (add_score) {
    sobj@meta.data = cbind(sobj@meta.data, scores_DF)
  }

  rm(scores_list)

  ## Attribute class according to higher score
  attribute = function(score_list) {
    which_one = which(score_list == max(score_list))[1] # in case of egality, keep first
    score_name = scores_names[which_one]
    newgroup = base::substr(score_name, start = nchar(prefix) + 1, stop = nchar(score_name))
    return(newgroup)
  }

  sobj@meta.data[, newname] = as.character(NA)
  sobj@meta.data[, newname] = apply(scores_DF, 1, attribute)
  sobj@meta.data[, newname] = as.factor(sobj@meta.data[, newname])

  ## Output
  return(sobj)
}
