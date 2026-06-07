#' Fit stability
#'
#' This is typically not called directly, but invoked by `fit_stability_network()`.
#'
#' @param X,y        As usual in linear regression
#' @param prior      Numeric vector of prior support Features with prior support (character vector)
#' @param strength   Ratio of penalty weights (no peak / peak). Default: 2.
#' @param N          Number of subsamples. Default: 100.
#' @param seed       For reproducibility.
#' @param bar        Display a progress bar? Default: TRUE.
#'
#' @return A list:
#' * `stb`: Percent of time a feature was selected at optimal lambda
#' * `beta`: Estimated coefficient at optimal lambda
#' * `m`: The penalty multipliers (from prior)
#'
#' @examples
#' X <- matrix(rnorm(50*10), 50, 10, dimnames=list(1:50, LETTERS[1:10]))
#' y <- rnorm(50)
#' prior <- rep(1:0, times=c(3,7))
#' stability(X, y, prior)$Max
#' stability(X, y, prior, 50)$Max
#' @export

stability <- function(X, y, prior, strength=1, N=100, seed=1, bar=TRUE) {

  # Set up
  X <- scale(X)
  y <- (y - mean(y))/sd(y)
  m <- strength_formula(prior, strength)
  n <- nrow(X)
  p <- ncol(X)

  # Fit full data model
  fit <- glmnet(X, y, penalty.factor=m)
  s <- sapply(predict(fit, type="nonzero"), length)
  lam <- fit$lambda[s <= (p/2)]
  if (bar) pb <- txtProgressBar(0, N, style=3)
  SS <- array(NA, dim=c(N, ncol(X), length(lam)), dimnames=list(1:N, colnames(X), lam))
  Y1 <- array(NA, dim=c(N, n, length(lam)), dimnames=list(1:N, 1:n, lam))
  Y2 <- array(NA, dim=c(N, n, length(lam)), dimnames=list(1:N, 1:n, lam))
  set.seed(seed)
  for (i in 1:N) {
    ind <- sample(1:n, ceiling(n/2))
    fit.i <- glmnet(X[ind,], y[ind], lambda=lam, penalty.factor=m)
    Y1[i,,] <- y - predict(fit.i, X)
    Y2[i,-ind,] <- y[-ind] - predict(fit.i, X[-ind,])
    SS[i,,] <- as.matrix(coef(fit.i)[-1,])
    if (bar) setTxtProgressBar(pb, i)
  }
  if (bar) close(pb)
  err <- 0.5*apply(Y1, 3, mean) + 0.5*(apply(Y2, 2:3, mean, na.rm=TRUE) |> apply(2, mean, na.rm=TRUE))
  Pct <- apply(SS != 0, 2:3, mean)
  return(list(stb=Pct[,which.min(err)], beta=coef(fit, as.numeric(names(which.min(err))))[-1], m=m))
}

strength_formula <- function(x, s) {
  ifelse(x == -Inf, Inf, 1 / (1 + x*s))
}
