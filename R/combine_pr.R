#' Compare multiple precision-recall curves
#'
#' @param files     List of `.rds` network files to combine
#' @param gs        A gold standard, or list of gold standards, as produced by `make_gs()`
#' @param rlim      Maximum value of recall to consider for partial auc (default: 0.25)
#' @param collapse  Indicates if the PR table should collapse across TFs or keep them separate
#'
#' @seealso `pr_plot()`
#'
#' @returns A data.table with columns for Recall, Specificity, Precision, Subset, Prior, and Strength
#' Also contains the attribute `random`, which gives the precision of random guessing.
#'
#' @examples
#' grhl3_gs <- make_gs('data/knockout/grhl3_KO.rds', 'lfsr  < 0.01')
#' combine_pr('network/gl/fns-mall-s0.rds', grhl3_gs)
#' combine_pr(c('network/gl/fns-mall-s0.rds', 'network/gl/fns-mall-s2.rds'), grhl3_gs)
#'
#' gs <- list(grhl3_gs, make_gs('data/knockout/irf6_KO.rds', 'lfsr  < 0.01'))
#' combine_pr('network/gl/fns-mall-s0.rds', gs)
#' combine_pr(c('network/gl/fns-mall-s0.rds', 'network/gl/fns-mall-s2.rds'), gs)
#' @export

combine_pr <- function(files, gs, rlim=0.25, collapse=FALSE) {

  if (inherits(gs, 'goldStandard')) gs <- list(gs)

  if (collapse) {
    gs_many <- suppressWarnings(many_tf_gs(gs))
    out <- combine_pr_manyTF(files, gs_many, rlim)
  } else {
    pr_list <- NULL
    for (g in gs) {
      pr <- lapply(files, function(rds) pr_network(rds, g, rlim, plot=FALSE)[-1])
      for (i in 1:length(pr)) {
        pr[[i]] <- pr[[i]][is.finite(Precision) & Recall != 0]
      }
      pr_list <- append(pr_list, pr)
    }
    out <- rbindlist(pr_list, fill=TRUE)
  }
  out
}
