#' Plot (overlay) multiple precision-recall curves
#'
#' @param DT        Output from `combine_pr()` (data.table)
#' @param group     What to group by (string)
#' @param facet     What to facet by (string, optional)
#' @param col       If faceting, the number of columns (default: 2)
#'
#' @returns   A `ggplot2` object
#'
#' @examples
#' gs <- list(make_gs('data/knockout/grhl3_KO.rds', 'lfsr  < 0.01'),
#'            make_gs('data/knockout/irf6_KO.rds', 'lfsr  < 0.01'))
#' combine_pr(c('network/mir/fns-mall-s0.rds', 'network/mir/fns-mall-s1.rds',
#'              'network/mir/fns-mall-s2.rds'), gs[[1]]) |> pr_plot('Strength')
#' combine_pr(c('network/gl/fns-mall-s0.rds', 'network/gl/fns-mall-s1.rds',
#'              'network/gl/fns-mall-s2.rds'), gs[[1]]) |> pr_plot('Strength')
#' combine_pr(c('network/mir/fns-mall-s0.rds', 'network/mir/fns-mall-s1.rds'), gs) |>
#'   pr_plot('Strength', 'TF')
#' @export

pr_plot <- function(DT, group, facet, col = 2) {

  # attach a table of auc total, partial, and strength to this plot
  if (!is.factor(DT[[group]])) {
    DT <- copy(DT)
    set(DT, j=group, value=factor(DT[[group]]))
  }
  p <- ggplot(DT, aes_string('Recall', 'Precision', Group=group, color=group)) +
    geom_line()
  # + ggtitle(title)
  if (missing(facet)) {
    p <- p + geom_hline(aes(yintercept=(Random)), lty=2, col='darkgray')
  } else {
    facet <- unlist(facet)
    RND <- DT[, .(Random=(Random)), c(facet)]
    p <- p + geom_hline(aes(yintercept=Random), RND, lty=2, col='darkgray') +
      facet_wrap(facet, ncol = col)
  }
  p
}
