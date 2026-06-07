#' Makes a gold standard object by combining two gold standard objects
#'
#' @param gs1     An object of class 'goldStandard'
#' @param gs2     An object of class 'goldStandard'
#'
#' @returns   A data.table of class 'goldStandard
#'
#' @examples
#' combine_gs(
#'   make_gs('data/knockout/grhl3_KO.rds', 'lfsr', 'maxval', 0.01),
#'   make_gs('data/chip/grhl3_chip_100kb.rds'))
#' @export

combine_gs <- function(gs1, gs2){
  gs <- merge(gs1, gs2, by = c('SYMBOL'))
  gs[, truth := ifelse((gs$truth.x & gs$truth.y == 1), 1, 0)]

  Flank <- c(gs[['Flank.x']][1], gs[['Flank.y']][1]) |>
    subset(!is.na(.)) |>
    unique() |>
    paste(collapse = ',')

  Threshold <- c(gs[['Threshold.x']][1], gs[['Threshold.y']][1]) |>
    subset(!is.na(.)) |>
    unique() |>
    paste(collapse = ',')

  TF <- c(gs[['TF.x']][1], gs[['TF.y']][1]) |>
    subset(!is.na(.)) |>
    unique() |>
    paste(collapse = ',')

  gs[, Flank := Flank]
  gs[, Threshold := Threshold]
  gs[, TF := TF]

  gs <- gs[, c('SYMBOL', 'truth', 'TF', 'Flank', 'Threshold')]

  class(gs) <- class(gs)

  return(gs)
}
