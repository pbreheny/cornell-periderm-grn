#' Plot a gene regulatory network with modules
#'
#' @param g             GRN graph, as produced by `make_graph()`
#' @param genes         optional gene list for genes to include in the network plot
#' @param n_lab         Number of hub genes to label per module (default: 5)
#' @param percentile    a number in (0, 1) indicating the percentile of edge weights to keep (default: 0.8)
#' @param resolution    A number in [0, Inf). Controls how many clusters will be returned. Lower values yield fewer, larger clusters.
#' @param nodeType      Possible values "color" (default) or "shape". Controls graph aesthetics.
#' @param label_size    Size of labels (default: 5)
#' @param label_border  Size of border around label (default: 0.1)
#' @param label_repel   Force of repulsion to keep labels from overlapping (default: 0.25)
#' @param primary       Use primary colors? (default: `FALSE`)
#' @param layout        a ggraph layout
#'
#' @examples
#' make_graph("network/wgcna/kna10-cor-s100.rds") |> plot_modules(nodeType = "color")
#' make_graph("network/wgcna/kna10-cor-s100.rds") |> plot_modules(nodeType = "shape")
#' @export

plot_modules <- function(
  g,
  genes,
  resolution = NULL,
  percentile = 0.8,
  n_lab = 5,
  nodeType = "color",
  label_size = 5,
  label_border = 0.1,
  label_repel = 0.25,
  primary = FALSE,
  layout = "fr",
  background = "dark"
) {
  # Set a seed
  original_seed <- .GlobalEnv$.Random.seed
  on.exit(.GlobalEnv$.Random.seed <- original_seed)
  set.seed(1)

  if (is.null(resolution)) {
    clust <- cluster_walktrap(g)
  } else {
    clust <- cluster_louvain(g, resolution = resolution)
  }

  res <- lapply(1:length(clust), function(x) {
    data.table(gname = clust[[x]], cluster = x)
  }) |>
    rbindlist()
  if (!missing(genes)) {
    ind <- which(names(V(g)) %in% genes)
    res <- res[gname %in% genes]
    g <- induced_subgraph(g, ind, "copy_and_delete")
    V(g)$module <- res$cluster
    V(g)$size <- eigen_centrality(g)$vector * 10
  } else {
    V(g)$module <- clust$membership
    V(g)$size <- eigen_centrality(g)$vector * 10
  }

  weights <- E(g)$weight
  w <- quantile(weights, percentile)

  clrs <- MetBrewer::met.brewer(
    "Tsimshian",
    length(unique(clust$membership)),
    "discrete"
  )
  shapes <- c(15, 17, 19, 18, 7, 8, 9, 10, 11, 12, 13, 14)

  label <- get.vertex.attribute(g) |>
    as.data.table() |>
    _[order(module, -size)] |>
    _[, head(.SD, n_lab), by = module]

  v_lab <- ifelse(V(g)$name %in% label$name, V(g)$name, "")

  if (nodeType == "shape") {
    p <- ggraph(g, layout = layout) +
      geom_edge_link0(aes(
        alpha = weights * 0.8,
        width = E(g)$width,
        filter = weights >= w,
        color = sign
      )) +
      scale_edge_color_discrete(h = c(128, 330) + 12, c = 100, l = 50) +
      theme_graph() +
      geom_node_point(
        aes(size = V(g)$size, alpha = 0.6),
        col = "black",
        shape = shapes[V(g)$module]
      ) +
      scale_size(range = c(1, 7))
  } else if (nodeType == "color") {
    p <- ggraph(g, layout = layout) +
      geom_edge_link0(aes(
        alpha = weights * 0.8,
        width = E(g)$width,
        filter = weights >= w,
        color = sign
      ))
    if (background == "dark") {
      p <- p +
        scale_edge_color_grey(start = 1, end = 0) +
        theme_graph(background = "#5F5F5F")
    } else if (background == "white") {
      p <- p +
        scale_edge_color_discrete(h = c(128, 330) + 12, c = 100, l = 50) +
        theme_graph(background = "white")
    }
    p <- p + scale_size(range = c(1, 7))
    if (primary) {
      p <- p +
        geom_node_point(
          aes(size = V(g)$size, alpha = 0.6, col = factor(V(g)$module)),
          shape = 19
        )
    } else {
      p <- p +
        geom_node_point(
          aes(size = V(g)$size, alpha = 0.6),
          col = clrs[V(g)$module],
          shape = 19
        )
    }
    p <- p +
      geom_node_text(
        aes(label = v_lab),
        repel = TRUE,
        size = label_size,
        color = "black",
        bg.color = "#FFFFFF80",
        bg.r = label_border,
        force = label_repel
      ) +
      guides(
        size = "none",
        alpha = "none",
        width = "none",
        edge_color = guide_legend("Direction of Association"),
        edge_alpha = "none"
      )
  }
  p
}
