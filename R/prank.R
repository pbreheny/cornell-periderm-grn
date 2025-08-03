#' Compute percentile rank
#' 
#' @param x           Either a numeric vector or a data.table
#' @param column      If a data.table is supplied for `x`, name of the column to calculate ranks for. If that data.table has a key, the elements of the returned vector will be named according to the key.
#' @param return_dt   By default, prank() returns a vector. If `return_dt = TRUE`, returns a data.table with one column for the key and one column for the percentile rank.
#' 
#' @examples
#' prank(airquality$Ozone)
#' DT <- as.data.table(airquality)
#' DT[, Date := paste(Month, Day, sep='-')]
#' setkey(DT, Date)
#' prank(DT, 'Ozone')
#' prank(DT, 'Ozone', return_dt=TRUE)
#' @export

prank <- function(x, column, return_dt=FALSE) {
  if (is.data.table(x)) {
    xx <- x[[column]]
    if (return_dt & haskey(x)) {
      out <- data.table(Key = x[[key(x)]],
                        V1 = rank(xx)/length(xx))
      setnames(out, 'V1', column)
      setnames(out, 'Key', key(x))
      setkeyv(out, key(x))
    } else {
      out <- rank(xx)/length(xx)
      names(out) <- x[[key(x)]]
    }
  } else {
    out <- rank(x)/length(x)
  }
  out
}
