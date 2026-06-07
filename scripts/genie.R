#' Fit genie networks
#'
#' Input: Raw data via import_expr()
#' Output: `network/genie/*`
source("scripts/setup.R")
library(GENIE3)

# Fit networks
# cr = correlation
# ts = topographical similarity
for (filter in c(0.10, 0.15, 0.20, 0.25)) {
  kna <- import_expr('kenny', 'none', 'all', filter=filter)
  id <- paste0('kna', round(filter*100))
  mc2 <- import_prior('mc2', kna)
  for (s in seq(0, 4, 0.25)) {
      file <- paste0('network/genie/', id, '-', 'gne', '-s', sprintf('%0.3d', as.integer(100*s)), '.rds')
      fit_genie_network(kna, mc2, file, strength=s)
  }
}
