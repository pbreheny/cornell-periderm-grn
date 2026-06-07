#' Calculate PR curve for network
#'
#' Optionally, plot the PR curve as well, although see also `pr_plot()`.
#'
#' @param rds         File name of saved network object
#' @param gs          an object of class 'goldStandard' make from `many_tf_gs()`
#' @param rlim        Maximum value of recall to consider for partial auc (default: 0.25)
#' @param plot        Make plot? Default: FALSE
#'
#' @return A data.table of Recall/Specificity/Precision, along with metadata (Subset/Prior/Strength)
#' to facilitate merging of multiple `pr_network()` calls with `combine_pr()`. Also
#' contains the attribute `random`, which gives the precision of random guessing.
#'
#' @examples
#' irf6_chip10kb_gs <- make_gs("data/chip/irf6_chip_10kb.rds")
#' tfap2a_chip10kb_gs <- make_gs("data/chip/tfap2a_chip_10kb.rds")
#' grhl3_chip10kb_gs <- make_gs("data/chip/grhl3_chip_10kb.rds")
#'
#' pr_network("network/mir/fns-mall-s0.rds", grhl3_chip10kb_gs) |> head()
#'
#' allChip_10kb_gs <- many_tf_gs(list(irf6_chip10kb_gs, grhl3_chip10kb_gs, tfap2a_chip10kb_gs))
#' pr_network("network/mir/fns-atac-s1.rds", allChip_10kb_gs) |> head()
#' @export

pr_network <- function(rds, gs, rlim = 0.25, plot = FALSE) {
  network <- readRDS(rds)

  # Create data.table
  tf <- unique(unlist(gs[["TF"]]))
  if (any(grepl("tfap2", tf))) {
    tf <- gsub("tfap2aMO", "tfap2a", tf)
    tf <- gsub("tfap2c", "tfap2a", tf)
    tf <- gsub("tfap2a,tfap2a", "tfap2a", tf)
  }

  cols <- c("target", tf)
  DT <- as.data.table(network)[, ..cols] |>
    melt(id.vars = "target") |>
    _[, value := abs(value)]

  # Merge GS and Network
  if (length(tf) == 1) {
    # DT <- merge(DT, gs, by = 'target')
    DT <- merge(DT, gs, by.x = "target", by.y = "SYMBOL")
  } else {
    DT <- DT[, id := paste(target, variable)] |>
      merge(gs, by = "id")
    # DT[,target.y := NULL]
    # setnames(DT, 'target.x', 'target')
  }
  DT <- DT[target != TF]

  # Calculate sensitivity, specificity, precision
  threshold <- sort(unique(DT$value))
  TP <- sapply(threshold, function(x) {
    sum(DT$value >= x & DT$truth == 1)
  })
  TN <- sapply(threshold, function(x) {
    sum(DT$value < x & DT$truth == 0)
  })
  FP <- sapply(threshold, function(x) {
    sum(DT$value >= x & DT$truth == 0)
  })
  FN <- sapply(threshold, function(x) {
    sum(DT$value < x & DT$truth == 1)
  })
  random <- mean(DT$truth)

  # Assemble into DT and collapse to unique recall points
  dir <- fifelse(
    basename(dirname(rds)) == "network",
    "",
    basename(dirname(rds))
  )
  out <- data.table(
    Network = rds,
    Recall = TP / (FN + TP),
    Specificity = TN / (FP + TN),
    Precision = TP / (TP + FP),
    Random = mean(DT$truth)
  )[, last(.SD), Recall]
  auc <- aupr(
    precision = out$Precision,
    recall = out$Recall,
    random = out$Random
  )
  out[, AUC := auc[[1]]]
  # out[, AUC_Partial := auc[[2]]]
  out[, AUC_Gain := auc[[2]]]
  out[, pAUC := auc[[3]]]
  out[, pAUC_Gain := auc[[4]]]
  if (length(tf) == 1) {
    out[, TF := tf]
  } else {
    TF <- tf |>
      subset(!is.na(.)) |>
      unique() |>
      paste(collapse = ",")
    out[, TF := TF]
  }

  # Plot
  if (plot) {
    x <- out$Recall[-1]
    y <- out$Precision[-1]
    ylim <- range(c(y, out$Random))
    plot(
      x,
      y,
      type = "l",
      xlab = "Recall",
      ylab = "Precision",
      bty = "n",
      las = 1,
      lwd = 2,
      ylim = ylim
    )
    abline(h = out$Random, v = rlim, col = "gray", lty = 2)
  }

  # Format metadata
  ncol_init <- ncol(out)
  meta <- attr(network, "meta")
  for (j in 1:length(meta)) {
    if (length(meta[[j]]) != 1) {
      next
    }
    if (inherits(meta[[j]], "name")) {
      next
    }
    set(out, j = names(meta)[j], value = meta[[j]])
  }
  if (!("alpha" %in% names(out))) {
    out$alpha <- 1
  }
  if (!("type" %in% names(out))) {
    out$type <- ""
  }
  if (!("strength" %in% names(out))) {
    out$strength <- 0
  }
  new_order <- c(1:ncol_init, ncol_init + order(names(out)[-(1:ncol_init)]))
  setcolorder(out, new_order)
  return(invisible(out))
}
