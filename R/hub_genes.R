#' Extract list of hub genes from a gene regulatory network
#'
#' @param g           GRN graph, as produced by `make_graph()`
#' @param n           Number of genes to include; if omitted, all genes are included
#' @param focus       TFs to focus on; these are colored blue while all other TFs are colored red
#'
#' @examples
#' make_graph('network/genie/fns-genie.rds') |> hub_genes
#' make_graph('network/gl/kna-h3k-s4.rds') |> hub_genes(15)
#' make_graph('network/gl/kna-h3k-s4.rds') |> hub_genes(15, focus=c('irf6', 'grhl3'))
#' @export

hub_genes <- function(g, n, focus) {
  x <- eigen_centrality(g)$vector |> sort(decreasing=TRUE)
  top <- data.table(Gene=names(x), X_=x)
  setnames(top, 'X_', 'Centrality')
  top[, Prank := prank(Centrality)]
  if (!missing(n)) top <- head(top, n)

  if (!missing(focus)) {
    focus <- data.table(Gene=focus, X_=x[focus], Rank=rank(x)[focus], pRank=prank(x)[focus])
    setnames(focus, 'X_', 'Centrality')
    list(top=top, focus=focus)
  } else {
    top
  }
}
