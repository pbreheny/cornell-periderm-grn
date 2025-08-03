#' Compare a prior to gold standard and print a summary
#' 
#' @param prior   The prior
#' @param gs      The gold standard
#' @param tf      The transcription factor; by default, the TF in the gs metadata
#' @param thresh  If prior is continuous (not 0/1)
#' @param quiet   Suppress printing
#' 
#' @return Invisibly returns a vector of percentages and ratios; needs to be labeled better
#' 
#' @examples
#' Data <- import_expr('kenny', 'none', 'all')
#' h3k <- import_prior('h3k', Data)
#' gs <- make_gs('data/chip/irf6_chip_10kb.rds')
#' compare_prior_gs(h3k, gs)
#' @export

compare_prior_gs <- function(prior, gs, tf, thresh=0.5, quiet=FALSE) {
  if (missing(tf)) {
    if (is.list(gs$TF[1])) {
      tf <- gs$TF[[1]][1]
    } else {
      tf <- gs$TF[1]
    }
  }
  p <- prior[, c('Gene', tf), with=FALSE]
  setnames(p, 'Gene', 'SYMBOL')
  merged <- merge(gs, p, by='SYMBOL')
  setnames(merged, 'SYMBOL', 'Gene')
  setnames(merged, tf, 'Prior')
  setnames(merged, 'truth', 'GS')
  merged[, Prior := 1 * (Prior > thresh)]
  Tab <- with(merged, table(Prior, GS))
  p0 <- (Tab[2,1] + Tab[2,2])/sum(Tab)
  p1 <- Tab[2,2]/(Tab[2,1] + Tab[2,2])
  p2 <- Tab[1,2]/(Tab[1,1] + Tab[1,2])
  p3 <- Tab[2,2]/(Tab[1,2] + Tab[2,2])
  p4 <- Tab[2,1]/(Tab[1,1] + Tab[2,1])
  if (!quiet) {
    print(Tab)
    cat('Fraction of genes w/ accessible motif: ', sprintf('%.1f', 100*p0), '%\n', sep='')
    cat('Fraction of genes w/ accessible motif that are bound:  ', sprintf('%.1f', 100*p1), '%\n', sep='')
    cat('Fraction of genes w/o accessible motif that are bound: ', sprintf('%.1f', 100*p2), '%\n', sep='')
    cat('Fraction of genes w/ binding sites that have accessible motif:  ', sprintf('%.1f', 100*p3), '%\n', sep='')
    cat('Fraction of genes w/o binding sites that have accessible motif: ', sprintf('%.1f', 100*p4), '%\n', sep='')
    cat('Prior enrichment: ', sprintf('%.2f', p1/p2), '\n', sep='')
  }
  return(invisible(c(p1/p2, p0, p1, p2, p3, p4)))
}
