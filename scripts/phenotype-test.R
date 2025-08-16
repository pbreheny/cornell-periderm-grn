#' Produces phenotype
#'   grhl1 (1)
#'   grhl3 (2)
#'   klf6a (3)
#'   klf17 (4)
#'   tead3b (7)
#'   irf6 (20) - this one is debatable for several reasons
#' Does not produce phenotype
#'   cebpb (5)
#'   ikzf2 (6)
#'   tfap2a (10)
#'   gata3 (12)
#'   zbtb10 (38)
#'   znf750 (51)

# Excluding IRF6
wilcox.test(c(1, 2, 3, 4, 7), c(5, 6, 10, 12, 38, 51))

# Including IRF6
wilcox.test(c(1, 2, 3, 4, 7, 20), c(5, 6, 10, 12, 38, 51))

