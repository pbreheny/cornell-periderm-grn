#' Extract feature matrix from central data.table
#' 
#' @param DT          Expression values (data.table with genes as key)
#' @param rowid       Rows to extract, returned object will be ordered accordingly. Default: all. (character vector)
#' @param colid       Columns to extract, returned object will be ordered accordingly. Default: all. (character vector)
#' @param transpose   Transpose matrix before returning?
#' 
#' @return A matrix
#' 
#' @export

dtx <- function(DT, rowid, colid, transpose=FALSE) {
  if (missing(rowid)) rowid <- DT[[1]]
  if (missing(colid)) colid <- names(DT)[-1]
  X <- DT[rowid][, ..colid] |> as.matrix
  rownames(X) <- rowid
  colnames(X) <- colid
  if (transpose) X <- t(X)
  X
}
