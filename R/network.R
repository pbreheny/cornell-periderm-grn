#' Create `network` object
#'
#' Intended to be called at the end of `fit_cor_network()`, `fit_gl_network()`, etc.
#'
#' @param x       A matrix of edge confidence weights; can be positive or negative
#'
#' @returns A network object, which is the matrix of edge confidence weights but
#' with an attribute `meta` that contains metadata on the inputs and arguments
#' used to fit the network:
#'   * `method`: If called by `fit_cor_network()`, this equals `cor` and so on
#'   * `samples`: Sample ids used to fit network (character vector)
#'   * `source`: As in `import_expr()`
#'   * `preprocess`: As in `import_expr()`
#'   * `subset`: As in `import_expr()`
#'   * `filter`: As in `import_expr()`
#'   * `prior`: `id` used in `import_prior()`
#'   * Additional arguments that were passed to `fit_x_network()`
#'
#' @examples
#' Data <- import_expr(filter=0.25)
#' Prior <- import_prior('m100', Data)
#' s <- 1
#' Res <- fit_cor_network(Data, Prior, strength=s)
#' Res[1:10,]
#'
#' # Handles language objects as arguments
#' Res <- fit_cor_network(Data, Prior, strength=s, tom=('t' %in% 'tom'))
#' Res[1:10,]
#' @export

network <- function(x) {
  Data <- get('Data', parent.frame())
  Prior <- get('Prior', parent.frame())
  args <- list()
  if (sys.nframe() == 1) {
    warning('network() is intended to be called inside of fit_xxx_network()')
    method <- '???'
  } else {
    call <- sys.calls()[[sys.nframe()-1]] |> as.list()
    fmls <- call[[1]] |> deparse() |> formals()
    for (item in setdiff(names(fmls), c('Data', 'Prior', 'filename'))) {
      if (item %in% names(call)) {
        args[[item]] <- eval(call[[item]])
      } else if (is.symbol(fmls[[item]])) {
        args[[item]] <- ''
      } else {
        args[[item]] <- fmls[[item]]
      }
    }
    method <- deparse(call[[1]]) |>
      str_extract('_(.*)', 1) |>
      str_replace('_network', '')
  }
  meta <- append(list(
    method = method,
    samples = names(Data$Exp)[-1],
    source = Data$source,
    preprocess = Data$preprocess,
    subset = Data$subset,
    filter = Data$filter,
    prior = attr(Prior, 'id')
  ), args)
  structure(x, class = 'network', meta = meta)
}

#' @export
print.network <- function(x, ...) {
  glue('Matrix of dimension {nrow(x)} x {ncol(x)} (targets x transcription factors)\n\n') |> cat()
  meta <- attr(x, 'meta')
  for (item in setdiff(names(meta), c('samples'))) {
    if (is.symbol(meta[[item]])) next
    sprintf('%-11s: %s\n', item, meta[[item]]) |> cat()
  }
}

#' @export
`[.network` <- function(x, i, j) {
  sub <- unclass(x)[i,j]
  if (is.array(sub)) {
    structure(sub, class = 'network', meta = append(attr(x, 'meta'), list(subnet=TRUE)))
  } else {
    sub
  }
}

#' @export
`as.data.table.network` <- function(x, long=FALSE) {
  if (long) {
    as.data.table(unclass(x), keep.rownames='target') |>
      melt(id.vars='target',
           variable.name = "tf",
           value.name = "confidence")
  } else {
    as.data.table(unclass(x), keep.rownames='target')
  }
}

#' @export
`as.matrix.network` <- function(x) {
  matrix(as.numeric(x), nrow=nrow(x), ncol=ncol(x), dimnames=dimnames(x))
}

#' @export
transform.network <- function(x, thresh=0.1, power=2, ...) {
  structure(
    soft(unclass(x), thresh)^power,
    class = 'network',
    meta = attr(x, 'meta'))
}
