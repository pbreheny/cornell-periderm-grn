#' Plot MCC vs confidence for a network
#' 
#' @param rds         File name of saved network, as fit by `fit_stability_network()`, for example
#' @param gs          Gold standard, as returned by `make_gs()`
#' @param tf          Name of the TF that is being considered (default: obtained from gs object)
#'
#' @examples
#' irf6_gs <- make_gs('data/chip/irf6_chip_10kb.rds')
#' plot_mcc('network/gl/kna-atac-s0.rds', irf6_gs)
#' @export

plot_mcc <- function(rds, gs, tf) {
  if (missing(tf)) tf <- gs$TF[1]

  # Create data.table
  network <- readRDS(rds)
  DT <- data.table(
    SYMBOL = rownames(network),
    Edge = network[, tf]) %>% 
    merge(gs, by='SYMBOL')

  # Plot
  x <- seq(min(DT$Edge), max(DT$Edge), len=101)[c(-1, -101)]
  y <- vapply(x, function(a) cor(DT$Edge >= a, DT$truth), 1.0)
  qplot(x, y) + xlab('Confidence') + ylab('MCC') +
    ggtitle(paste0('Maximum MCC: ', sprintf('%.2f', max(y)), ', edges: ', sum(DT$Edge >= x[which.max(y)])))
}
