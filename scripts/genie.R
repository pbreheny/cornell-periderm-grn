#' Fit genie networks
#'
#' Input: Raw data via import_expr()
#' Output: `network/genie/*`
source("scripts/setup.R")
library(GENIE3)
filter <- 0.10
dir.create("network/genie", recursive = TRUE, showWarnings = FALSE)

# Fit networks
kna <- import_expr('kenny', 'none', 'all', filter = filter)
mc2 <- import_prior('mc2', kna)
for (s in seq(0, 4, 0.25)) {
  file <- glue(
    'network/genie/\\
        kna{sprintf("%0.2d", round(100*filter))}\\
        -mc2-gne-\\
        s{sprintf("%0.3d", as.integer(100*s))}.rds'
  )
  fit_genie_network(kna, mc2, file, strength = s)
}
