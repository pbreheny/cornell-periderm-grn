#' Plot a gene regulatory network
#'
#' @param g             GRN graph, as produced by `make_graph()`
#' @param n             Number of genes to label (default: 10)
#' @param percentile    a number in (0, 1) indicating the percentile of edge weights to keep
#' @param max_width     maximum edge width
#' @param max_node      maximum node size
#'
#' @examples
#' make_graph('data/network/genie/kna10-gne-s100.rds') %>% plot_graph()
#' @export

plot_graph <- function(g, n=10, percentile = 0.7, max_width = 0.1, max_node = 6) {
  weights <- E(g)$weight
  w <- quantile(weights, percentile)
  importance <- eigen_centrality(g)$vector
  V(g)$size <- importance
  ggraph(g, layout = "stress") +
    geom_edge_link0(aes(alpha=weights, width=weights, filter = weight >= w)) +
    scale_edge_width(range = c(0.01, max_width)) +
    theme_graph() +
    geom_node_point(aes(size = size), col='#3FA9FF', shape=19) +
    scale_size(range = c(1, max_node)) +
    theme(legend.position = "none") +
    geom_node_shadowtext(aes(filter = size >= sort(size, decreasing=TRUE)[n], size=size/1.5, label = name), family = "sans", fontface='bold', color='black', bg.color='white')
}

geom_node_shadowtext <- function (mapping = NULL, data = NULL, position = "identity",
                                  parse = FALSE, nudge_x = 0, nudge_y = 0, check_overlap = FALSE,
                                  show.legend = NA, repel = FALSE, ...) {
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      stop("Specify either `position` or `nudge_x`/`nudge_y`",
           call. = FALSE)
    }
    position <- position_nudge(nudge_x, nudge_y)
  }
  geom <- shadowtext::GeomShadowText
  params <- list(parse = parse, na.rm = FALSE, ...)
  mapping <- ggraph:::aes_intersect(mapping, aes(x = x, y = y))
  layer(data = data, mapping = mapping, stat = StatFilter,
        geom = geom, position = position, show.legend = show.legend,
        inherit.aes = FALSE, params = params)
}

