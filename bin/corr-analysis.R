cor.per.cell <- function(pipe.mtx, illumina){
  whitelist_cells <- colnames(illumina)
  pipe.cells <- colnames(pipe.mtx)
  common.genes <- intersect(rownames(pipe.mtx), rownames(illumina))
  cor.list <- list() #vector(length=length(whitelist_cells))
  idx = 1
  for (cell in whitelist_cells){
    if (cell %in% pipe.cells){
      r.coef <- cor(pipe.mtx[common.genes, cell], illumina[common.genes, cell])
      cor.list[[cell]] = r.coef #cor.list[idx] = r.coef
    }else{
      cor.list[[cell]] = NaN #cor.list[idx] = NaN
    }
    idx = idx + 1
  }
  return(cor.list)
}

getVariableName <- function(var) {
  return(deparse(substitute(var)))
}