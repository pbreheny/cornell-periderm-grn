#' Processes expression data for use in lasso/genie/wcgna analyses
#'
#' @param source       Source of expression data
#' @param preprocess   Use a preprocessed version of the data? One of:
#'   * 'none': No preprocessing
#'   * 'imp'; MAGIC used to impute 0s away
#'   * 'ccr': Cell cycle regression used to correct for mitosis effects
#' @param subset       For 'farrell', use shield stage or all samples?
#' @param filter       Expression filter (default: 0.1, meaning genes with 10% expression or less are removed)
#' @param qc           Print QC checks?
#'
#' @returns An object of class `expr`, with elements \describe{
#' \item{Exp}{A matrix of expression values}
#' \item{Tissue}{If `source="farrell"`, the tissue assigned to each sample}
#' \item{tf}{List of transcription factors (character vector)}
#' \item{source}{See above}
#' \item{subset}{See above}
#' \item{filter}{See above}
#' \item{preprocess}{See above}
#' \item{pct}{Named vector containing percent of cells expressing each gene}
#' }
#'
#' @examples
#' import_expr()
#' import_expr("farrell", subset = "shield")
#' @export

import_expr <- function(
  source = c("kenny", "farrell", "kennyOld"),
  preprocess = c("none", "imp", "ccr", "both"),
  subset = c("all", "shield"),
  filter = 0.1,
  qc = FALSE
) {
  source <- match.arg(source)
  preprocess <- match.arg(preprocess)
  subset <- match.arg(subset)

  tf <- fread("data/raw/zf_All_TF_names_Curated.txt")[[1]] |>
    str_to_lower() |>
    sort()
  if (source == "farrell") {
    # Tissue
    DT1 <- fread("data/raw/EarlyEVL.txt", header = FALSE)
    DT1[, Tissue := "EarlyEVL"]
    DT2 <- fread("data/raw/EVL.txt", header = FALSE)
    DT2[, Tissue := "EVL"]
    DT3 <- fread("data/raw/LateEVL.txt", header = FALSE)
    DT3[, Tissue := "LateEVL"]
    DT4 <- fread("data/raw/Epidermis.txt", header = FALSE)
    DT4[, Tissue := "Epidermis"]
    Tissue <- rbind(DT1, DT2, DT3, DT4)
    setnames(Tissue, "V1", "Sample")
    setkey(Tissue, "Sample")

    # Merge with time
    Time <- readxl::read_excel(
      "data/raw/Time_series_data4GeneExpression.xlsx"
    ) |>
      as.data.table()
    setnames(Time, "Sample Name", "Sample")
    setnames(Time, "TimePoint", "Time")
    Tissue <- merge(Tissue, Time, by = "Sample")
    Tissue[, Sample := tolower(Sample)]
    setkey(Tissue, "Sample")

    # Exp
    if (preprocess == "imp") {
      Exp <- Read("data/raw/zf-rna-tpm.txt")
    } else if (preprocess == "none") {
      Exp <- Read("data/raw/evlSC_rawExpressionMatrix.txt")
    } else if (preprocess == "ccr") {
      Exp <- Read("data/raw/evlSC_CellCycleRegressed_rawData.txt")
    } else if (preprocess == "both") {
      Exp <- Read("data/raw/MAGIC_processedkNN8_evlSC_CellCycleRegressed.txt")
    }
  } else if (source == "kenny") {
    Tissue <- NULL
    # Exp <- fread('data/raw/TPM_aggrCK1CK2_ShieldRNA.txt')
    Exp <- readxl::read_excel(
      "data/raw/InputDataAll_12July23_Sunil.xlsx",
      sheet = "TPM_GeneExpression",
      .name_repair = "minimal"
    ) |>
      as.data.table()
    samples <- names(Exp)[-1] |>
      str_replace("\\..*", "") |>
      str_to_lower()
    names(Exp) <- c("Gene", samples)
    Exp[, Gene := str_to_lower(Gene)]
    Exp[, (samples) := lapply(.SD, log1p2), .SDcols = samples]
    setkey(Exp, "Gene")
  } else if (source == "kennyOld") {
    Tissue <- NULL
    Exp <- fread("data/raw/TPM_aggrCK1CK2_ShieldRNA.txt")
    samples <- names(Exp)[-1] |>
      str_replace("\\..*", "") |>
      str_to_lower()
    names(Exp) <- c("Gene", samples)
    Exp[, Gene := str_to_lower(Gene)]
    Exp[, (samples) := lapply(.SD, log1p2), .SDcols = samples]
    setkey(Exp, "Gene")
  }
  # Fix 'sep' genes
  Exp[, Gene := str_replace(Gene, "^sep(\\d+)$", "sept\\1")]

  # Subset
  if (subset == "shield") {
    samples <- Tissue[Time > 4.5 & Time < 6.01]$Sample
    Exp <- Exp[, c("Gene", samples), with = FALSE]
    Tissue <- Tissue[samples]
    setkey(Tissue, "Sample")
  }

  # Apply expression filter
  DT <- data.table(
    gene = Exp$Gene,
    pct = apply(Exp[, -1] > 0, 1, mean),
    tf = Exp$Gene %in% tf
  )
  pct <- DT$pct
  names(pct) <- DT$gene
  # tf_include <- DT[tf == TRUE & avg > 0]$gene
  # target_include <- DT[tf == FALSE & avg > 0][order(avg, decreasing=TRUE)][1:n]$gene
  pct <- DT$pct
  names(pct) <- DT$gene
  tf_include <- DT[tf == TRUE & pct > filter]$gene
  target_include <- DT[tf == FALSE & pct > filter]$gene

  # QC (do any of this stuff?)
  if (qc) {
    print("qc not set up yet")
    # sum(!(tf %in% Exp$Gene))
    # setdiff(tf, Exp$Gene)
    # sum(!(target %in% Exp$Gene))
    # setdiff(target, Exp$Gene)
    # target <- intersect(union(target, tf), Exp$Gene)
    # tf <- intersect(tf, Exp$Gene)
    # dim(Tissue)
    # dim(Exp)
    # any(duplicated(Tissue$Gene))
  }

  # Return
  Exp <- Exp[Gene %in% c(tf_include, target_include)]
  setkey(Exp, "Gene")
  structure(
    list(
      Exp = Exp,
      Tissue = Tissue,
      tf = sort(tf_include),
      source = source,
      preprocess = preprocess,
      subset = subset,
      filter = filter,
      pct = pct
    ),
    class = "expr"
  )
}

#' @export
print.expr <- function(x, ...) {
  cat("Expression matrix contains:\n")
  cat(
    "  ",
    nrow(x$Exp),
    " targets: ",
    paste(x$Exp$Gene[1:3], collapse = ", "),
    ", ...\n",
    sep = ""
  )
  cat(
    "  ",
    ncol(x$Exp) - 1,
    " samples: ",
    paste(names(x$Exp)[2:3], collapse = ", "),
    ", ...\n",
    sep = ""
  )
  cat(
    "  ",
    length(x$tf),
    " transcription factors: ",
    paste(x$tf[1:3], collapse = ", "),
    ", ...\n",
    sep = ""
  )
  if (!is.null(x$Tissue)) {
    cat(
      "Metadata for: ",
      paste(names(x$Tissue)[-1], collapse = ", "),
      "\n",
      sep = ""
    )
  }
  cat("source      : ", x$source, "\n", sep = "")
  cat("preprocess  : ", x$preprocess, "\n", sep = "")
  cat("subset      : ", x$subset, "\n", sep = "")
}

tpm <- function(x, size, log2 = FALSE) {
  y <- 1e6 * (x / size) / sum(x / size)
  if (log2) log1p(y) / log2 else y
}

log1p2 <- function(x) {
  log1p(x) / log(2)
}
