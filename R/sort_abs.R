#' Sort by absolute value
#'
#' @param x   A numeric vector to be sorted
#'
#' @examples
#' sort_abs(c(3, 0, -2, 1))
#' @export

sort_abs <- function(x, decreasing=TRUE) {
  ord <- order(abs(x), decreasing=decreasing)
  x[ord]
}
