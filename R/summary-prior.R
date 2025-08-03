#' Summarize a prior
#' 
#' @param object   An object of class 'prior', as created by `import_prior()`
#' @param thresh   A connection has support from the prior if its entry is `>= thresh` (default: 1). For a 0-1 prior, there is no point in changing this parameter, but it is useful for fractional/continuous priors (see example).
#' @param ...      For S3 compatibility; not used
#' 
#' @examples
#' Data <- import_expr()
#' Prior <- import_prior('cmb', Data)
#' summary(Prior)       # Support from both H3K and ATAC
#' summary(Prior, 0.5)  # Support from either H3K or ATAC
#' @method summary prior
#' @export

summary.prior <- function(object, thresh=1, ...) {
  P <- object[,-1]
  Edges <- apply(P >= thresh, 2, sum)
  tf <- names(object)
  TF <- apply(P[object$Gene %in% tf,] >= thresh, 2, sum)
  Target <- apply(P[!(object$Gene %in% tf),] >= thresh, 2, sum)
  out <- data.table(
    Gene=names(Edges),
    Edges, TF, Target,
    edgRank = prank(Edges), trfRank = prank(TF), trgRank = prank(Target))
  setkey(out, 'Gene')
  out
}
