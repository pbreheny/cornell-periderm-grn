#' Fit inferelator network
#'
#' @param Data       Object created by `import_expr()`
#' @param Prior      An object of class 'prior', as created by `import_prior()`
#' @param filename   If provided, will store the resulting `network` object as an `rds` file rather than return it
#' @param type       'bbsr' or 'enet' [string]
#' @param strength   Weight given to prior in the analysis (bbsr only); 0 means no weight (all features penalized equally)
#' @param alpha      Elastic net mixing parameter (enet only)
#' @param reps       Number of replications (different cv seeds; default: 5)
#' @param clean      Remove raw files (from python) after fitting? (default: TRUE)
#'
#' @returns An object of class `network`
#'
#' @examples
#' Data <- import_expr(filter = 0.5)
#' Prior <- import_prior("mc2", Data)
#' Res <- fit_inferelator_network(Data, Prior, type = "bbsr", strength = 1, reps = 2)
#' @export

fit_inferelator_network <- function(
  Data,
  Prior,
  filename,
  type = c("bbsr", "enet"),
  strength = NULL,
  alpha = NULL,
  reps = 5,
  clean = TRUE
) {
  type <- match.arg(type)
  prefix <- stringi::stri_rand_strings(1, 12)

  # Write data for python
  dir.create(glue("network/inferelator/{prefix}"), recursive = TRUE, showWarnings = FALSE)
  fwrite(Data$Exp, glue("network/inferelator/{prefix}-exprs.tsv"), sep = "\t")
  writeLines(Data$tf, glue("network/inferelator/{prefix}-tf.txt"))
  fwrite(Prior, glue("network/inferelator/{prefix}-prior.tsv"), sep = "\t")

  # Run python script
  if (type == "bbsr") {
    cmd <- glue("scripts/fit-inferelator.py {prefix} bbsr {strength} {reps}")
  } else {
    cmd <- glue("scripts/fit-inferelator.py {prefix} enet {alpha} {reps}")
  }
  status <- system(cmd)
  if (status != 0) stop("Error running python code: ", status)

  # Process inferelator output
  x <- matrix(0, length(Data$Exp$Gene), length(Data$tf), dimnames = list(Data$Exp$Gene, Data$tf))
  s <- matrix(0, length(Data$Exp$Gene), length(Data$tf), dimnames = list(Data$Exp$Gene, Data$tf))
  for (rep in 1:reps) {
    DT <- fread(glue("network/inferelator/{prefix}/random_seed_{rep}/network.tsv.gz"))
    DT <- DT[target %in% rownames(x) & regulator %in% colnames(x) & !is.na(combined_confidences)]
    for (j in 1:nrow(DT)) {
      x[DT$target[j], DT$regulator[j]] <- x[DT$target[j], DT$regulator[j]] +
        DT[j]$combined_confidences
      s[DT$target[j], DT$regulator[j]] <- s[DT$target[j], DT$regulator[j]] + DT[j]$beta.sign.sum
    }
  }
  out <- network(sign(s) * x / reps)

  # Clean up
  if (clean) {
    unlink(glue("network/inferelator/{prefix}"), recursive = TRUE)
    unlink(glue("network/inferelator/{prefix}-exprs.tsv"))
    unlink(glue("network/inferelator/{prefix}-prior.tsv"))
    unlink(glue("network/inferelator/{prefix}-tf.txt"))
  }

  # Save / return
  if (missing(filename)) {
    invisible(out)
  } else {
    saveRDS(out, filename)
    cat(crayon::green$bold("Output saved:", filename, "\n"))
  }
}
