#' Box plot but whiskers show quantiles instead of outliers
#'
#' @param df
#' @param x
#' @param g
#' @param orientation
#'
#' @examples
#'
#' @export

qbox <- function(df, x, g, orientation=c('horizontal', 'vertical')) {
  orientation <- match.arg(orientation)
  dat <- eval(substitute(as.data.table(df)[, .(
    min = quantile(x, 0.1),
    low = quantile(x, 0.25),
    mid = quantile(x, 0.5),
    hi  = quantile(x, 0.75),
    max = quantile(x, 0.9)),
    g]))
  if (orientation == 'horizontal') {
      ggplot(dat, aes(y={{g}}, fill={{g}})) +
      geom_boxplot(aes(xmin=min, xlower=low, xmiddle=mid, xupper=hi, xmax=max), stat = "identity") +
      theme_minimal() +
      theme(legend.position = "none")
  } else {
    ggplot(dat, aes(x = {{g}}, fill = {{g}})) +
      geom_boxplot(aes(ymin=min, lower=low, middle=mid, upper=hi, ymax=max), stat = "identity") +
      theme_minimal() +
      theme(legend.position = "none")
  }
}
