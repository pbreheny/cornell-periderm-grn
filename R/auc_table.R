#' Print a table of auc values for specified pr curves
#'
#' @param DT        A data.table as returned by `combine_pr()`
#'
#' @returns         A `data.table` object
#'
#' @examples
#' gs <- list(make_gs('data/knockout/grhl3_KO.rds', 'lfsr', 'maxval', 0.01),
#'            make_gs('data/knockout/irf6_KO.rds', 'lfsr','maxval', 0.01))
#' f <- list.files('data/network/wgcna/', pattern='kna10-mc2-cor*', full.names=TRUE)
#' DT <- combine_pr(f, gs)
#' auc_table(DT)
#' @export

auc_table <- function(DT) {
  if (inherits(DT$TF, 'list')) DT$TF <- vapply(DT$TF, paste, '', collapse='/')
  out <- DT[, first(.SD), .(Network, TF, method)]
  out[, Recall := NULL]
  out[, Specificity := NULL]
  out[, Precision := NULL]
  out
}
