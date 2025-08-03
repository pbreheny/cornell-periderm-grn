#' Combine gold standards from multiple TFs
#' 
#' @param gs_list   A list of gold standards, as produced by `make_gs()`
#' 
#' 
#' @returns A data table of class 'goldStandard'
#' 
#' @examples
#' irf6_chip100kb_gs <- make_gs('data/chip/irf6_chip_100kb.rds')
#' tfap2a_chip100kb_gs <- make_gs('data/chip/tfap2a_chip_100kb.rds')
#' grhl3_chip100kb_gs <- make_gs('data/chip/grhl3_chip_100kb.rds')
#' 
#' many_tf_gs(list(irf6_chip100kb_gs, grhl3_chip100kb_gs, tfap2a_chip100kb_gs))
#' @export

many_tf_gs <- function(gs_list) {
  n <- length(gs_list)
  x <- vector('list', n)
  for (i in 1:n) {
    x[[i]] <- copy(gs_list[[i]])
    x[[i]][, id := paste(SYMBOL, unique(x[[i]]$TF))]
  }
  gs <- rbindlist(x)
  class(gs) <- c('gold standard', class(gs))
  return(gs)
}
