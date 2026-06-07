#' Make a (potentially interactive) plot of TF cluster assignments
#'
#' @param net         A `network` object
#' @param panel       Show
#' @param resolution  As in `summary.network()`
#' @param plotly      Convert to plotly (enables hover tooltips; default=`FALSE`)
#' @param boring      Use boring colors? (default: `FALSE`)
#'
#' @examples
#' readRDS('network/wgcna/kna10-cor-s100.rds') |> cluster_plot()
#' \dontrun{
#' readRDS('network/wgcna/kna10-cor-s100.rds') |> cluster_plot(boring=TRUE) + hrbrthemes::theme_ft_rc()
#' }
#' @export

cluster_plot <- function(net, panel=FALSE, resolution=NULL, plotly=FALSE, boring=FALSE) {
  DT <- summary(net, resolution=resolution)
  if (panel) {
    p <- ggplot(DT, aes(Pct, Enr, color=Cluster, size=Panel, shape=Panel, text=Gene)) +
      scale_size_manual(values=c(2.5,6)) +
      theme(legend.position="none")
  } else {
    p <- ggplot(DT, aes(Pct, Enr, color=Cluster, text=Gene))
  }
  p <- p + geom_point() +
    theme_minimal() +
    xlab("Expression (% of EVL cells)") +
    ylab("Enrichment (log2 fold change, EVL vs other tissues)")
  if (boring) p <- p + scale_colour_manual(values = MetBrewer::met.brewer("Thomas", uniqueN(DT$Cluster), 'discrete'))
  if (plotly) {
    plotly::ggplotly(p, tooltip="text")
  } else {
    p
  }
}
