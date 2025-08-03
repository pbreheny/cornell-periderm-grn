#' Read in ATAC Priors and return the genes with positive ATAC peaks
#' 
#' @param prior_dat       ATAC-Seq prior matrix (each row is a gene, each column is a transcription factor)
#' @param TF              The transcription factor name of interest (corresponds to a column name in prior_dat)
#' 
#' @return A vector of strings corresponding to names of genes marked as having an ATAC-Seq peak

read_priors <- function(prior_dat, TF){
  fullPrior <- read.delim(prior_dat, header = TRUE)
  TF_Prior <- fullPrior[, c('X', TF)]
  pos_genes <- TF_Prior[which(TF_Prior[, TF] == 1) , 'X']
  return(pos_genes)
}
