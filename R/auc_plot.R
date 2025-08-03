#' Plot AUC values
#'
#' @param DT          Output from `auc_table()` (data.table)
#' @param group       What to group by (string)
#' @param facet       What to facet by (string, optional)
#' @param auc_metric  What for of AUC to plot on the y-axis (string)
#' @param ncol        How many columns should be in the facet
#'
#' @returns   A `ggplot2` object
#'
#' @examples
#' gs <- list(make_gs('data/knockout/grhl3_KO.rds', 'lfsr < 0.01'),
#'            make_gs('data/knockout/irf6_KO.rds', 'lfsr < 0.01'))
#' f <- list.files('data/network/gl/', pattern='fns-atac*', full.names=TRUE)
#' DT <- combine_pr(f, gs)
#' auc_table(DT) %>% auc_plot('Strength', 'TF')
#' @export

auc_plot <- function(DT, group, facet, auc_metric, ncol = NULL) {
  DT <- copy(DT)

  p <- ggplot(DT, aes(x = strength, y = get(auc_metric), group = get(group), color = get(group))) +
    labs(color = group)

  # Facet if necessary
  if (!missing(facet)) {
    if (is.null(ncol)){
      ncol = length(unique(DT[[facet]]))
    }
    p <- p + geom_line() + facet_wrap(facet, scales = 'free_x', ncol = ncol)
  } else{
    p <- p + geom_line()
  }

  ylabel <- str_replace(auc_metric, "_", " ")
  p <- p + ylab(ylabel) + xlab("Prior Strength")
  p
}
