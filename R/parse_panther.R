#' Parse JSON output from panther
#'
#' @param json_file   JSON filename
#'
#' @seealso `scripts/panther-export.r`
#'
#' @examples
#' parse_panther('data/go/cluster-a.json')
#' @export

parse_panther <- function(json_file) {
  jsonlite::fromJSON(json_file)[[1]]$group$result |>
    lapply(fix_panther_single) |>
    rbindlist()
}

fix_panther_single <- function(x) {
  if (is.data.frame(x)) {
    x$mapped_id_list <- NULL
    as.data.table(x)
  } else {
    if (length(x$term) != 3) {
      x$term <- list(level='', id='', label='')
    }
    data.table(number_in_list = x$number_in_list,
               fdr            = x$fdr,
               pValue         = x$pValue,
               term.level     = x$term$level,
               term.id        = x$term$id,
               term.label     = x$term$label,
               plus_minus     = x$plus_minus)
  }
}
