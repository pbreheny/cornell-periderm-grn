#' Fit wgcna networks
#'
#' Input:
#'   - Raw data via import_expr()
#'   - Prior via import_prior()
#' Output: `network/wgcna/*`
devtools::load_all(quiet = TRUE)

# Fit networks
# cor = correlation
# tom = topographical similarity
for (filter in c(0.10, 0.15, 0.20, 0.25)) {
  kna <- import_expr("kenny", "none", "all", filter = filter)
  mc2 <- import_prior("mc2", kna)
  for (s in seq(0, 4, 0.25)) {
    file <- glue(
      'network/wgcna/\\
        kna{sprintf("%0.2d", round(100*filter))}\\
        -mc2-cor-\\
        s{sprintf("%0.3d", as.integer(100*s))}.rds'
    )
    fit_cor_network(kna, mc2, file, s, type = "cor")
  }
}
