#' Subset an expr object
#'
#' @param x        An object of class 'expr', as created by `import_expr()`
#' @param subset   One of the following:
#'   * A logical vector indicating whether to keep a sample or not
#'   * A vector of row indices
#'   * A vector of sample names to keep
#'   * An expression that evaluates to TRUE/FALSE when evaluated in x$Tissue
#' @param ...      For S3 compatibility; not used
#'
#' @examples
#' Data <- import_expr()
#' subset(Data, 1:100)
#' subset(Data, names(Data$Exp)[11:51])
#' subset(Data, rep(c(TRUE, FALSE), c(30, ncol(Data$Exp)-31)))
#' Data <- import_expr('farrell')
#' subset(Data, Tissue == 'EVL')  # Only works with farrell
#' @method subset expr
#' @export

subset.expr <- function(x, subset, ...) {
  subset <- try(eval(match.call()$subset, envir=parent.frame()), silent=TRUE)
  if (inherits(subset, 'try-error')) subset <- try(eval(match.call()$subset, envir=x$Tissue), silent=TRUE)
  if (inherits(subset, 'try-error')) stop("Can't evaluate subset expression", call. = FALSE)
  if (is.integer(subset)) subset <- (2:ncol(x$Exp) - 1) %in% subset
  if (is.character(subset)) subset <- names(x$Exp)[-1] %in% subset
  if (!is.logical(subset)) stop('Subset must be logical, integer, or character', call. = FALSE)
  if (!is.null(x$Tissue)) {
    copy(x$Tissue)
    x$Tissue <- subset(x$Tissue, subset)
  }
  copy(x$Exp)
  x$Exp <- subset(x$Exp, select=c(TRUE, subset))
  x
}
