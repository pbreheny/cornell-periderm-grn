#' Calculate area under PR curve for a network
#'
#' @param precision       Calculated and ordered precision values from a GS object
#' @param recall          Calculated and ordered recall values from a GS object
#' @param random          Height of the random guess line for a PR curve
#' @param upper           Upper limit of integration (default: maximum value of recall prior to 1)
#' @param rlim            Maximum value of recall to consider for partial auc (default: 0.25)
#'
#' @return a list with total AUC and partial AUC (AUC - area under random guess)
#'
#' @examples
#' p <- sort(runif(30), decreasing=TRUE)
#' r <- sort(runif(30))
#' plot(r, p, las=1); abline(h=0.1, lty=2)
#' aupr(p, r, 0.1)
#' @export

aupr <- function(precision, recall, random, upper=max(recall[recall != 1]), rlim=0.25) {

  partial <- which(recall <= rlim)
  pupper <- max(recall[partial])
  f1 <- approxfun(recall, precision, ties=mean)
  f2 <- approxfun(recall, precision - random, ties=mean)
  f3 <- approxfun(recall[partial], precision[partial], ties=mean)
  f4 <- approxfun(recall[partial], precision[partial] - random[partial], ties=mean)
  auc_tot <- integrate(f1, min(recall), upper, stop.on.error=FALSE)$value
  auc_gain <- integrate(f2, min(recall), upper, stop.on.error=FALSE)$value
  pauc_tot <- integrate(f3, min(recall), pupper, stop.on.error=FALSE)$value
  pauc_gain <- integrate(f4, min(recall), pupper, stop.on.error=FALSE)$value

  return(list('Total AUC' = auc_tot, 'AUC Gain' = auc_gain, 'pAUC' = pauc_tot, 'pAUC Gain' = pauc_gain))
}
