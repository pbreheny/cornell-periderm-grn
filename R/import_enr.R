#' Import enrichment information
#'
#' @examples
#' enr <- import_enr()
#' @export

import_enr <- function() {
  #f <- 'data/raw/ShieldStageMultiome_aggrCK5CK6_DE_EVLvsAll.xlsx'
  f <- 'data/raw/TargetGeneEnrichment_in_ShieldStageMultiome_aggrCK5CK6_1Sep23.xlsx'
  out <- readxl::read_xlsx(f, progress=FALSE, .name_repair='minimal') %>%
    as.data.table()
  out[, Gene := tolower(Gene)]
  setkey(out, 'Gene')
  out
}
