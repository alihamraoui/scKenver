library(dplyr)
library(tidyr)
library(ggpubr)
library(patchwork)

plot_scatter <- function(x, y, title, col){
  x= log(x)
  x[x == -Inf] <- 0
  y= log(y)
  y[y == -Inf] <- 0
  plot(x, y, main = title,
       xlab = "True Log2 TPM", 
       ylab = "Estimated Log2 TPM",
       ylim=c(-5,15),
       xlim=c(-5,15),
       pch = 19, cex = .2, frame = FALSE, col = col)
  abline(lm(y ~ x, data = mtcars), col = "black")
  text(8, 13, paste("r =", round(cor(x,y),2)))
}

filter_cols <- function(mtx){
  return(mtx[,which(!colnames(mtx) %in%
                      c("Geneid","Start","End","Chr","Strand","Length", "transcript_id", "gene_name"))])
}

is_gene_id <- function(genes) {
  return(grepl("^ENS[G|T|P|C|G][0-100]{9,}", genes))
}

mtx.to.gene <- function(mtx){
  mtx <- mtx %>% filter_cols()
  
  mtx <- data.table::as.data.table(mtx)
  
  numeric_cols <- which(sapply(mtx, is.numeric))
  
  mtx <- as.data.frame(mtx[, lapply(.SD, sum), by = gene_id, .SDcols = numeric_cols])
  
  mtx <- tibble::column_to_rownames(mtx, var = "gene_id")
}

is.cellbarcode <- function(mtx){
  sapply(colnames(mtx), function(x) all(grepl("^[ATGC]+$", x)))
}

rename.genes <- function(mtx, gtf){
  genes.gtf = data.frame(GenomicRanges::mcols(gtf)[,c("gene_id","gene_name")])
  
  founded.geneNames <- data.frame(gene_id=genes.gtf$gene_id[match(rownames(mtx), genes.gtf$gene_id, nomatch = NA)])
  
  if (all(is.na(founded.geneNames$gene_id))){
    founded.geneNames <- data.frame(gene_id=genes.gtf$gene_id[match(rownames(mtx), genes.gtf$gene_name, nomatch = NA)])
  }
  
  mtx$gene_id <- founded.geneNames$gene_id
  
  mtx <- aggregate(mtx[,-length(mtx)],
                   by = list(mtx$gene_id),
                   FUN = sum)
  rownames(mtx) <- mtx$Group.1
  mtx <- mtx[,-1]
  return(mtx)
}


filter.trns.mtx <- function(mtx){
  if ("transcript_id" %in% colnames(mtx)){
    mtx <- mtx %>%
      `rownames<-`(.$transcript_id)   %>%
      dplyr::select(-transcript_id)
  }else if ("transcriptId" %in% colnames(mtx)){
    mtx <- mtx %>% 
      mutate(transcriptId = ifelse(transcriptId == "undef", paste0(geneId, "_undef"), transcriptId)) %>%
      `rownames<-`(.$transcriptId)   %>%
      dplyr::select(-transcriptId)
  } else if ("X" %in% colnames(mtx)){
    mtx <- mtx %>% separate(X, c('gene_name', 'transcript_id'), sep="_") %>%
      `rownames<-`(.$transcript_id)   %>%
      dplyr::select(-c(transcript_id, gene_name))}
  mtx <- mtx[,is.cellbarcode(mtx)]
}

normalize <- function(x, mn, sd) {
  n <- ncol(x)
  for (i in 1:n) {
    std_dev <- sd
    if (std_dev == 0) std_dev <- 1
    x[, i] <- (x[, i] - mn) / sd
  }
  return(x)
}