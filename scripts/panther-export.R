#' Export data for panther testing
#' http://www.pantherdb.org/
#'   - Danio rerio
#'   - Statistical enrichment test
#'   - PANTHER GO-Slim Biological Process
#'
#' Overrepresentation: Binary, unsatisfying
#' Enrichment: Continuous, better
source("scripts/setup.R")

# Targets by cluster
net <- readRDS("network/wgcna/mc2/kna10-cor-s100.rds")
g <- transform(net) |> make_graph()
gc <- cluster_walktrap(g)
m <- membership(gc)
sizes(gc)

# Overrepresentation
if (interactive()) {
  a <- apply(net[, names(m[m == 1])] > 0.25, 1, any)
  b <- apply(net[, names(m[m == 2])] > 0.25, 1, any)
  c <- apply(net[, names(m[m == 3])] > 0.25, 1, any)
  cat(names(a), sep = "\n", file = "ref.txt")
  names(a[a]) |> clipr::write_clip()
  names(b[b]) |> clipr::write_clip()
  names(c[c]) |> clipr::write_clip()
}

# Enrichment
export <- function(j) {
  x <- apply(net[, names(m[m == j])], 1, max)
  data.table(Gene = names(x), Score = x) |>
    fwrite(glue::glue("data/go/cluster-{letters[j]}.txt"), sep = "\t")
}
lapply(1:3, export) |> invisible()
