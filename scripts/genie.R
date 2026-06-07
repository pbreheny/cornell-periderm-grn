suppressPackageStartupMessages(devtools::load_all())
library(GENIE3)

# Farrell (none)
Data <- import_expr('farrell', 'none', 'shield')
X <- dtx(Data$Exp)
Res <- GENIE3(X, regulators=Data$tf, nCores=4)
attr(Res, 'method') <- 'genie'
attr(Res, 'nTrees') <- 1000
attr(Res, 'samples') <- colnames(X)
class(Res) <- c('network', class(Res))
saveRDS(t(Res), 'network/fns-genie.rds')

# Farrell (ccr)
Data <- import_expr('farrell', 'ccr', 'shield')
X <- dtx(Data$Exp)
Res <- GENIE3(X, regulators=Data$tf, nCores=4)
attr(Res, 'method') <- 'genie'
attr(Res, 'nTrees') <- 1000
attr(Res, 'samples') <- colnames(X)
class(Res) <- c('network', class(Res))
saveRDS(t(Res), 'network/fcs-genie.rds')

# Kenny
Data <- import_expr('kenny')
X <- dtx(Data$Exp)
Res <- GENIE3(X, regulators=Data$tf, nCores=4)
attr(Res, 'method') <- 'genie'
attr(Res, 'nTrees') <- 1000
attr(Res, 'samples') <- colnames(X)
class(Res) <- c('network', class(Res))
saveRDS(t(Res), 'network/kna-genie.rds')


# new way
#' Fit genie networks
#'
#' Input: Raw data via import_expr()
#' Output: `network/genie/*`
devtools::load_all(quiet=TRUE)

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
