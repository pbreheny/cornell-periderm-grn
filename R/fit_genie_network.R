#' Fit correlation network
#'
#' @param Data       Object created by `import_expr()`
#' @param Prior      An object of class 'prior', as created by `import_prior()`
#' @param filename   If provided, will store the resulting `network` object as an `rds` file rather than return it
#' @param strength   Weight given to prior in the analysis; 0 means no weight (all features penalized equally)
#'
#' @returns An object of class `network`
#'
#' @examples
#' Data <- import_expr(filter=0.6)
#' Prior <- import_prior('m100', Data)
#' Res <- fit_genie_network(Data, Prior, strength=0.5)
#' @export

fit_genie_network <- function(Data, Prior, filename, strength=1) {

  # Fit
  X <- dtx(Data$Exp) |> GENIE3(regulators=Data$tf, nCores=4)

  # Boost
  x <- t(X)
  if (strength > 0) {
    P <- dtx(Prior)
    x <- x + strength * P * x
  }

  # Save / return
  out <- network(x)
  if (missing(filename)) {
    return(invisible(out))
  } else {
    saveRDS(out, filename)
  }
}
