#' Estimate transcription factor activity
#' 
#' Not all target genes necessarily have entries in the prior, so the intersection of these two is taken.
#' 
#' @param Data     List containing TPM and target genes; see process-for-lasso.R.
#' @param Prior    Data table prior, as returned by `Read()`. Each row is a target gene, each column is a TF.
#' @param weight   Weight given to prior
#' 
#' @return A matrix of estimated transcription factor activities. Each row is a sample, each column is a TF.

estimate_activity <- function(Data, Prior, weight=0) {
  if (missing(Prior)) {
    id <- sort(intersect(Data$Exp$Gene, Data$tf))
  } else {
    id <- sort(intersect(Data$Exp$Gene, names(Prior)[-1]))
  }
  B <- dtx(Data$Exp, id)
  if (weight==0) {
    return(t(B))
  } else {
    X <- dtx(Data$TPM, id)
    A <- MASS::ginv(P) %*% X
    rownames(A) <- colnames(P)
    return(weight*t(A) + (1-weight)*t(B))
  }
}
