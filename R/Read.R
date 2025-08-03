#' Import text file to data table, where row names are genes
#' 
#' @param input   File to read (character string)
#' 
#' @return A data table, keyed to 'Gene'. There is a row for each target gene, and a column for each transcription factor (plus one column for the key / gene name).
#' @keywords internal

Read <- function(input) {
  out <- fread(input)[V1 != '']
  names(out) <- tolower(names(out))
  setnames(out, 'v1', 'Gene')
  out[, Gene := tolower(Gene)]
  setkey(out, 'Gene')
  out
}
