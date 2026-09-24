#' Makes a gold standard object from processed KO or ChIP data
#'
#' @param data            RDS file of processed ChIP or KO data
#' @param metric          argument for KO data indicating the metric and threshold
#' @param selectby method of setting gold standard either by a metric maximum (maxval) or total number (N)
#' @param threshold threshold value for the 'metric'
#' @param targets optional vector of targets to include in the gold standard. If not supplied, all available genes will be included
#' @param tie     For `metric = 'lfsr'` and `selectby = 'N'`: lfsr values below `tie` are treated as tied, and
#'   those genes are ranked by absolute posterior mean instead (default: 1e-4)
#'
#' @returns   A data.table of class 'goldStandard'
#'
#' @examples
#' make_gs('data/knockout/grhl3_KO.rds', 'lfsr', 'maxval', 0.01)
#' make_gs('data/knockout/grhl3_KO.rds', 'lfsr', 'N', 250)
#' @export

make_gs <- function(data, metric = c('lfsr', 'svalue', 'lfdr', 'qvalue'), selectby = c('maxval', 'N'), threshold = NULL, targets = NULL, tie = 1e-4) {
  metric <- match.arg(metric)
  selectby <- match.arg(selectby)
  dat <- readRDS(data)

  if(!is.null(targets)){
    dat <- dat[SYMBOL %in% targets]
  }

  if ("KO" %in% class(dat)) {
    if (selectby == 'maxval') {
      stopifnot("threshold must be between 0 and 1" = between(threshold, 0, 1))
      dat[, truth := ifelse(dat[[metric]] < threshold, 1, 0)]
    } else if (selectby == 'N') {
      stopifnot("threshold must be a positive integer" = (threshold %% 1 == 0) && (threshold >= 1))
      if (metric == 'lfsr') {
        # lfsr is effectively tied near zero; rank those genes by |posterior mean|, the rest by lfsr
        tied <- dat$lfsr < tie
        dat <- dat[order(!tied, ifelse(tied, -abs(PosteriorMean), lfsr))]
      } else {
        dat <- dat[order(dat[[metric]])]
      }
      dat[, truth := c(rep(1, threshold), rep(0, nrow(dat) - threshold))]
    }
  } else if("chipSeq" %in% class(dat)) {
    setnames(dat, "True_Pos", "truth")
    metric <- NULL
  }

  flank <- attr(dat, 'Flank')
  gs <- dat[, c('SYMBOL', 'truth', ..metric)]

  if(attr(dat, 'TF') == 'tfap2a2c'){
    gs[, TF := 'tfap2a']
  }else{
    gs[, TF := attr(dat, 'TF')]
  }

  gs[, Flank := if (is.null(flank)) NA else flank]
  gs[, Threshold := if (is.null(threshold)) NA else threshold]
  class(gs) <- c("goldStandard", class(gs))
  return(gs)
}
