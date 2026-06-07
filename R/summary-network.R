#' Summarize a gene regulatory network
#'
#' @param object      An object of class 'network', as created by `fit_stability_network()` or `fit_cor_network()`
#' @param resolution  For clustering; `igraph::cluster_walktrap()` by default, but `igraph::cluster_louvain()` if `resolution` is supplied.
#' @param tf_only     Return information on TFs only? (default: true)
#'
#' @method summary network
#'
#' @examples
#' readRDS('network/wgcna/kna10-mc2-cor-s100.rds') |> summary()
#' @export

summary.network <- function(object, resolution=NULL, tf_only=TRUE, ...) {
  # Set a seed
  original_seed <- .GlobalEnv$.Random.seed
  on.exit(.GlobalEnv$.Random.seed <- original_seed)
  set.seed(1)

  # Edge sum
  e <- abs(object) |> rowSums()
  edge <- data.table(Gene=names(e), Edgesum=e, key='Gene')

  # Cluster
  if (tf_only) {
    g <- make_graph(object)
  } else {
    g <- make_graph(object, rownames(object))
  }
  if (is.null(resolution)) {
    m <- cluster_walktrap(g) |> membership
  } else {
    m <- cluster_louvain(g, resolution=resolution) |> membership
  }
  cluster <- data.table(
    Gene=names(m),
    Cluster=factor(m),
    TF=names(m) %in% colnames(object),
    key='Gene')

  # Raw expression
  meta <- attr(object, 'meta')
  Data <- import_expr(meta$source, meta$preprocess, meta$subset)
  pct <- data.table(Gene=names(Data$pct), Pct=Data$pct, key='Gene')

  # Enrichment
  enr <- import_enr()[, .(Gene, Enr=avg_log2FC)]

  # Gene panel
  panel <- import_gp()
  gp <- data.table(Gene = rownames(object), key='Gene')
  gp[, Panel := factor(Gene %in% panel$Gene)]

  # Return
  hub_genes(g) |>
    .[, .(Gene, Centrality)] |>
    merge(edge) |>
    merge(cluster) |>
    merge(pct) |>
    merge(enr) |>
    merge(gp)
}
