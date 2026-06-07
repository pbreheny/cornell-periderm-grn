#' Fit network with stability subsampling to quantify edge confidence
#'
#' @param Data       Object created by `import_expr()`
#' @param Prior      An object of class 'prior', as created by `import_prior()`
#' @param filename   If provided, will store the resulting `network` object as an `rds` file rather than return it
#' @param subset     Which samples to include? (character vector of row names) Default: all.
#' @param N          Number of subsamples
#' @param strength   Weight given to prior in the analysis; 0 means no weight (all features penalized equally)
#' @param weight     Passed to `estimate_activity()`
#' @param thin       If `thin=10`, only fit network for 10 random target genes; for debugging purposes.
#'
#' @returns An object of class `network`
#'
#' @examples
#' Data <- import_expr()
#' Prior <- import_prior('atac', Data)
#' Res <- fit_stability_network(Data, Prior, N=10, thin=10)
#' @export

fit_stability_network <- function(Data, Prior, filename, subset, N=100, strength=2, weight=0, thin) {
  A <- estimate_activity(Data, Prior, weight)
  Y <- dtx(Data$Exp, transpose=TRUE)
  if (!missing(subset)) {
    A <- A[subset, ]
    Y <- Y[subset, ]
  }
  pb <- txtProgressBar(0, ncol(Y), style=3)
  tt <- 1:ncol(Y)
  if (!missing(thin)) tt <- sample(tt, thin)

  # Fit network
  Stability <- matrix(NA, ncol(Y), ncol(A), dimnames=list(colnames(Y), colnames(A)))
  Beta <- matrix(NA, ncol(Y), ncol(A), dimnames=list(colnames(Y), colnames(A)))
  for (j in tt) {
    stopifnot(identical(names(Prior)[-1], colnames(A)))
    prior <- Prior[Gene==colnames(Y)[j], -1] |> as.numeric
    if (colnames(Y)[j] %in% colnames(A)) prior[which(colnames(A) == colnames(Y)[j])] <- -Inf
    if (sd(Y[,j]) < 1e-10) {
      Stability[j,] <- 0
      Beta[j,] <- 0
    } else {
      fit <- stability(A, Y[,j], prior, strength=strength, N=N, bar=FALSE)
      Stability[j,] <- fit$stb
      Beta[j,] <- fit$beta
    }
    setTxtProgressBar(pb, j)
  }
  close(pb)

  x <- 0.5*Stability + 0.5*abs(Beta)
  x[x < 0] <- 0

  # Return / save
  out <- network(x)
  if (missing(filename)) {
    return(invisible(out))
  } else {
    saveRDS(out, filename)
  }
}
