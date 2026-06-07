#' Fit basic lasso network (graphical lasso)
#'
#' @param Data       Object created by `import_expr()`
#' @param Prior      An object of class 'prior', as created by `import_prior()`
#' @param filename   If provided, will store the resulting `network` object as an `rds` file rather than return it
#' @param strength   Weight given to prior in the analysis; 0 means no weight (all features penalized equally)
#' @param alpha      Elastic net parameter (default=1, lasso)
#' @param thin       If `thin=10`, only fit network for 10 random target genes; for debugging purposes.
#'
#' @returns An object of class `network`
#'
#' @examples
#' Data <- import_expr()
#' Prior <- import_prior('mc2', Data)
#' Res <- fit_glasso_network(Data, Prior, thin=10)
#' @export

fit_glasso_network <- function(Data, Prior, filename, strength=2, alpha=1, thin) {
  A <- scale(estimate_activity(Data, Prior))
  Y <- dtx(Data$Exp, transpose=TRUE)
  pb <- txtProgressBar(0, ncol(Y), style=3)
  tt <- 1:ncol(Y)
  if (!missing(thin)) tt <- sample(tt, thin)

  # Fit network
  Res <- matrix(0, ncol(Y), ncol(A), dimnames=list(colnames(Y), colnames(A)))
  for (j in tt) {
    stopifnot(identical(names(Prior)[-1], colnames(A)))
    prior <- Prior[Gene==colnames(Y)[j], -1] |> as.numeric()
    if (sd(Y[,j]) < 1e-10) {
      Res[j,] <- 0
    } else {
      y <- (Y[,j] - mean(Y[,j]))/sd(Y[,j])
      if (colnames(Y)[j] %in% colnames(A)) prior[which(colnames(A) == colnames(Y)[j])] <- -Inf
      m <- strength_formula(prior, strength)
      fit <- glmnet(A, y, penalty.factor=m, alpha=alpha)
      l <- apply(coef(fit)[-1, ] != 0, 1, which.max); l[l==1] <- Inf
      p <- apply(coef(fit)[-1, ] > 0, 1, which.max); p[p==1] <- Inf
      n <- apply(coef(fit)[-1, ] < 0, 1, which.max); n[n==1] <- Inf
      lam <- fit$lambda[pmin(p,n)] * ifelse(p < n, 1, -1)
      lam[is.na(lam)] <- 0
      Res[j,] <- lam
    }
    setTxtProgressBar(pb, j)
  }
  close(pb)

  # Return / save
  out <- network(Res)
  if (missing(filename)) {
    return(invisible(out))
  } else {
    saveRDS(out, filename)
  }
}
