#' Import gene panel information
#'
#' @examples
#' gp <- import_gp()
#' @export

import_gp <- function() {
  f <- "data/raw/CleftGenePanel_comprehensive_with _Danio_IDs.xlsx"
  out <- readxl::read_xlsx(f, progress = FALSE) |> as.data.table()
  setnames(out, "Danio_orthologousGene_name", "Gene")
  out[, Gene := tolower(Gene)]
  setkey(out, "Gene")
  out
}
