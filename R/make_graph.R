#' Construct igraph object from a gene regulatory network
#'
#' Takes a GRN (gene regulatory network), as produced by `fit_stability_network()`,
#' and produces an `igraph` object. Requires the `igraph` package to be installed.
#' Only plots transcription factors. If a list of TFs is not supplied, attempts to
#' read from rds/data-raw.rds.
#'
#' @param grn             A gene regulatory network, as produced by `fit_stability_network()` (either filename or object itself)
#' @param genes           A list of genes to include (default: all TFs)
#' @param remove_loners   By default, TFs with no edges are omitted from the graph; to see them, set this to `TRUE`
#' @param weighted        By default, a weighted graph is constructed; to create an unweighted graph set this to `FALSE`
#' @param ...             Additional arguments to `igraph::graph_from_adjacency_matrix()`
#'
#' @examples
#' g <- make_graph('network/genie/kna10-mc2-s025.rds')
#' g <- make_graph('network/genie/fns-genie.rds', genes=c('grhl3', 'klf17', 'aagab', 'aar2'))
#' @export

make_graph <- function(grn, genes, remove_loners=TRUE) {

  # Setup
  if (is.character(grn)) {
    stages <- (strsplit(grn, '/')[[1]][2] %>% strsplit('-'))[[1]][1]
    net <- (strsplit(grn, '/')[[1]][2] %>% strsplit('-'))[[1]][2]
    strength <- (strsplit(grn, '/')[[1]][2] %>% strsplit('-'))[[1]][3]%>%substr(start = 2, stop = 2)
    grn <- readRDS(grn)
  } else {
    stages <- net <- strength <- NULL
  }
  if (missing(genes)) genes <- colnames(grn)

  # create adjacency matrix
  S <- matrix(0, length(genes), length(genes))  # Signed version
  gtmap <- match(genes, rownames(grn))
  if (any(is.na(gtmap))) stop('Some genes do not appear in network', call.=FALSE)
  gfmap <- match(genes, colnames(grn))
  S[,!is.na(gfmap)] <- grn[gtmap,gfmap[!is.na(gfmap)]]
  U <- abs(S)  # Unsigned version
  colnames(S) <- rownames(S) <- colnames(U) <- rownames(U) <- genes

  # Construct graph
  if (remove_loners) {
    n_edge <- apply(U, 1, sum)
    interest <- n_edge[n_edge > 0] %>% names
    S <- S[interest, interest]
    U <- U[interest, interest]
  }
  g <- graph_from_adjacency_matrix(U, mode = 'undirected', weighted = TRUE)
  s1 <- t(S)[rev(upper.tri(S))]
  s2 <- S[lower.tri(S)]
  s <- vapply(seq_along(s1), function(i) {
    if (abs(s1[i]) > abs(s2[i])) {
      return(s1[i])
    } else {
      return(s2[i])
    }
  }, numeric(1))

  sgn <- sign(s)
  E(g)$sign <- factor(sgn[s!=0])
  attr(g, 'stages') <- stages
  attr(g, 'network') <- net
  attr(g, 'strength') <- strength
  g
}
