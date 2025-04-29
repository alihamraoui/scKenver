library(ggplot2)
library(tidyverse)
#library(ggstatsplot)
library(hrbrthemes)
library(Seurat)
library(ggpubr)
library(UpSetR)
library(RColorBrewer)
library(reshape2)
#library(ggvenn)
#library(VennDiagram)
library(viridis)
library(biomaRt)
#library(ggplotify)
#library(rtracklayer)
#library(ggVennDiagram)
library(plyr)
library(purrr)
library(patchwork)
library(data.table)
# Search gene name from GTF
search.genes.names <- function(flames.gene.mtx, gtf){
  genes.gtf = data.frame(mcols(gtf)[,c("gene_id","gene_name","transcript_id")])

  # match geneID gene nomes
  founded.geneNames <- data.frame(genes.gtf$gene_name[match(flames.gene.mtx$gene_id, genes.gtf$gene_id, nomatch = NA)])
  names(founded.geneNames) <- "gene_id"
  flames.gene.mtx$gene_id <- founded.geneNames$gene_id

  return(flames.gene.mtx)
}


## translate gene id -> gene name.
name.genes <- function(flames.gene.mtx, gtf){

  # remplace gene ID to Gene Symbol
  flames.gene.mtx <- search.genes.names(flames.gene.mtx, gtf)
  flames.gene.mtx <- aggregate(flames.gene.mtx[,-c(1,2)], #-length(flames.gene.mtx)],
                               by = list(flames.gene.mtx$gene_id),
                               FUN = sum)
  rownames(flames.gene.mtx) <- flames.gene.mtx$Group.1
  flames.gene.mtx <- flames.gene.mtx[,-1]
  return(flames.gene.mtx)
}

name.genes2 <- function(flames.gene.mtx, gtf){
  genes.gtf = data.frame(mcols(gtf)[,c("gene_id","gene_name","transcript_id")])

  # match geneID gene nomes
  founded.geneNames <- data.frame(genes.gtf$gene_name[match(rownames(flames.gene.mtx), genes.gtf$gene_id, nomatch = NA)])
  names(founded.geneNames) <- "gene_id"

  # remplace gene ID to Gene Symbol
  flames.gene.mtx$gene_id <- founded.geneNames$gene_id
  flames.gene.mtx <- aggregate(flames.gene.mtx[,-length(flames.gene.mtx)],
                               by = list(flames.gene.mtx$gene_id),
                               FUN = sum)
  rownames(flames.gene.mtx) <- flames.gene.mtx$Group.1
  flames.gene.mtx <- flames.gene.mtx[,-1]
  return(flames.gene.mtx)
}

# draw venn diagramm / return a grid object
draw.venn <- function (pipe.name, pipe.data) {
  venn_object = venn.diagram(
    x = list(pipe.name = pipe.data, 'Illumina'=data.bc_umi.illumina$V1),
    category.names = c(pipe.name, "Illumina "),
    fill = c('darkcyan', 'darkgoldenrod3'),
    print.mode=c("raw","percent"),
    lty = 'blank',
    #cex = .9,
    #cat.cex = 1,
    #cat.fontface = "bold",
    fontface = "italic",
    filename = NULL, #filename = 'venn-sicelore2.png',
    output=TRUE
  )
  return(venn_object)
}

# clustering / returns a seurat object
clustuerer <- function (matrix, min.features, nfeatures, dim) {
  ScObject = CreateSeuratObject(counts = matrix, min.cells = 3, min.features = min.features)
  ScObject = NormalizeData(ScObject, normalization.method = "LogNormalize", verbose = F) # , scale.factor = 10000)
  ScObject = FindVariableFeatures(ScObject, selection.method = "vst", nfeatures = nfeatures, verbose = F)
  all.genes <- rownames(ScObject)
  ScObject <- ScaleData(ScObject, features = all.genes, verbose = F)
  ScObject <- RunPCA(ScObject, features = VariableFeatures(object = ScObject), verbose = F)#, npcs = 10)
  ScObject <- FindNeighbors(ScObject, dims = 1:dim, verbose = F)
  ScObject <- FindClusters(ScObject, resolution = 0.5, verbose = F)
  ScObject <- RunUMAP(ScObject, dims = 1:dim, verbose = F)
}

# cluster annotation
annotate.aud <- function (ScObject) {
  if ("cell_type" %in% colnames(ScObject@meta.data)) {
    ScObject$cell_type = NULL
  }
  ScObject = cell_annot_custom(ScObject,
                               assay = "RNA",          # which assay to use to perform annotation ?
                               newname = "cell_type",  # new column name in meta.data
                               markers = cell_markers) # markers for cell type
  #summary(ScObject$cell_type)
  #DimPlot(ScObject, group.by = "cell_type", cols = color_markers)
  return(ScObject)
}

# scatters plot
sctp.gpc <- function  (pip){
  sp1 <- ggscatter(df2.gpc, x = "illumina", y = pip,
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE # Add confidence interval
  ) + stat_cor(method = "pearson")
  return(sp1)
}
