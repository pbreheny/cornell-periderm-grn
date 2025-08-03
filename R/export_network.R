#' Export csv version of network
#'
#' Truncates digits to save space.
#'
#' @param obj   A `network` object
#'
#' @examples
#' \dontrun{
#' readRDS('data/network/wgcna/mc2/kna10-cor-s100.rds') %>%
#'   export('net.csv', FALSE)
#' }
#' @export

export.network <- function(obj, file, as_list=TRUE, ...) {
  if (as_list) {
    g <- make_graph(obj, rownames(obj))
    out <- data.table(as_edgelist(g), weight=sprintf('%0.5f', E(g)$weight))
    names(out)[1:2] <- c('gene1', 'gene2')
  } else {
    out <- as.data.table(obj)
    for (i in seq_along(out)) {
      if (i == 1) next
      out[[i]] <- sprintf('%0.5f', out[[i]])
    }
  }
  if (missing(file)) {
    out
  } else {
    fwrite(out, file=file)
  }
}

#' @export
export <- function(obj, ...) UseMethod('export')
