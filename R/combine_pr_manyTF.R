#' Compare PR curves across networks for many TFs at once
#'
#' @param files       List of `.rds` network files to combine
#' @param gs          A list of gold standard objects
#' @param rlim        Maximum value of recall to consider for partial auc (default: 0.25)
#'
#' @return A data.table of Recall/Specificity/Precision.
#' Also contains the attribute `random`, which gives the precision of random guessing.
#'
#' @examples
#' irf6_chip10kb_gs <- make_gs('data/chip/irf6_chip_10kb.rds')
#' tfap2a_chip10kb_gs <- make_gs('data/chip/tfap2a_chip_10kb.rds')
#' grhl3_chip10kb_gs <- make_gs('data/chip/grhl3_chip_10kb.rds')
#'
#' allChip_100kb_gs <- suppressWarnings(many_tf_gs(list(irf6_chip10kb_gs, grhl3_chip10kb_gs, tfap2a_chip10kb_gs)))
#' combine_pr_manyTF(c('network/gl/fns-atac-s1.rds', 'network/gl/fns-atac-s4.rds'), allChip_100kb_gs)
#' @export

combine_pr_manyTF <- function(files, gs, rlim=0.25) {
  out <- lapply(seq_along(files), function(i) pr_network(files[i], gs, rlim, plot=FALSE)[-1]) %>%
    rbindlist(fill = TRUE)
  out <- out[is.finite(Precision) & Recall != 0]
  attr(out, 'random') <- mean(gs$truth)

  out
}

