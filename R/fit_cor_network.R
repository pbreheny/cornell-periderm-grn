#' Fit correlation network
#'
#' @param Data       Object created by `import_expr()`
#' @param Prior      An object of class 'prior', as created by `import_prior()`
#' @param filename   If provided, will store the resulting `network` object as an `rds` file rather than return it
#' @param strength   Weight given to prior in the analysis; 0 means no weight (all features penalized equally)
#' @param tom        Use topological overlap to construct network (default=FALSE)
#'
#' @returns An object of class `network`
#'
#' @examples
#' Data <- import_expr(filter=0.25)
#' Prior <- import_prior('m100', Data)
#' Res <- fit_cor_network(Data, Prior, strength=1)
#' @export

fit_cor_network <- function(Data, Prior, filename, strength=1, type=c('cor', 'tom')) {
  type <- match.arg(type)

  # Fit
  cr <- dtx(Data$Exp, transpose = TRUE) %>% cor()
  diag(cr) <- 0

  # TOM
  if (type == 'tom') {
    ts <- WGCNA::TOMsimilarity(abs(cr), verbose=0)
    rownames(ts) <- colnames(ts) <- rownames(cr)
    cr <- ts
  }

  # Boost
  x <- cr[, Data$tf]
  if (strength > 0) {
    P <- dtx(Prior)
    x <- x + strength * P * x
  }

  # Save / return
  out <- network(x)
  attr(out, "meta")$strength <- strength
  if (missing(filename)) {
    return(invisible(out))
  } else {
    saveRDS(out, filename)
  }
}
