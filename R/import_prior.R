#' Reads in a prior and zero-pads if any target genes are not included
#'
#' Prior should be scaled so that maximum value (maximum prior support) is 1.
#'
#' @param id      Data that the prior should be based on:
#'   * `atac`: Target has an accessible motif binding site
#'   * `h3k`: Target has an accessible motif binding site AND an active enhancer mark (intersects an H3K27ac peak)
#'   * `cmb`: Combines `atac` and `h3k` into a 0/0.5/1 prior
#'   * `exp`: Prior is based on expression, not motif (and thus not target-specific)
#'   * `ovr`: Prior based on overexpression in periderm
#'   * `scat`: Single-cell ATAC-Seq based on Kenny multiome data
#'   * `mall`: Multiome shield stage (all EVL)
#'   * `mspc`: Multiome shield stage (specific to EVL)
#'   * `m100`: Multiome shield stage (all EVL, any peak within 100kb, not just closest)
#'   * `mc2`: Multiome shield stage from JASPAR version 9 (all EVL, any peak within 100kb, not just closest)
#' @param Data    From `import_expr()`; prior will expand/subset to match target genes
#' @param imputeTF If TRUE, mean imputation is used for missing TFs. Otherwise, insert 0.
#'
#' @examples
#' Data <- import_expr()
#' m100 <- import_prior('m100', Data)
#' dim(m100); dim(Data$Exp)  # Rows match
#' m100
#' H3K <- import_prior('h3k', Data)
#' summary(m100)['irf6',]
#' summary(H3K)['irf6',]
#' EXP <- import_prior('exp', Data)
#' summary(EXP, 0.1)['irf6']
#' OVR <- import_prior('ovr', Data)
#' OVR[, c('Gene', 'arid3a', 'atf4a', 'znfl2a', 'zzz3')]
#' @export

import_prior <- function(id = c('mc2', 'atac', 'h3k', 'cmb', 'exp', 'ovr', 'scat', 'mall', 'mspc', 'm100'), Data, imputeTF = TRUE) {
  match.arg(id)
  target <- Data$Exp$Gene

  # Read in data
  atac_file <- 'data/raw/100kbBinaryShieldStageATAC_InferredAll_danRer10.tsv'
  h3k_file <- 'data/raw/corrected_BinaryShieldStageATACnShieldStageH3K27ac_FIMO_InferredAll.tsv'
  if (id == 'atac') {
    Prior <- Read(atac_file)
  } else if (id == 'h3k') {
    Prior <- Read(h3k_file)
  } else if (id == 'cmb') {
    P1 <- Read(atac_file)
    P2 <- Read(h3k_file)
    cmb_target <- intersect(P1$Gene, P2$Gene)
    Prior <- data.table(Gene=cmb_target, (P1[cmb_target][, -1] + P2[cmb_target][, -1])/2)
    setkey(Prior, 'Gene')
  } else if (id == 'exp') {
    f <- function(x, offset=0.01) (log(x+offset) - log(offset))/(log(1+offset) - log(offset))
    p <- apply(Data$Exp[Data$tf][, -1] != 0, 1, mean)
    names(p) <- Data$tf
    Prior <- data.table(Gene=target, as.data.table(as.list(f(p))))
    setkey(Prior, 'Gene')
  } else if (id == 'ovr') {
    cluster_data <- readxl::read_excel('data/raw/Cluster_Home_Genes_multiome.xls', 'Cluster31') |> as.data.table()
    ovr <- rep(0, length(Data$tf))
    names(ovr) <- Data$tf
    ovr[Data$tf %in% cluster_data$gene] <- 1
    Prior <- data.table(Gene=target, as.data.table(as.list(ovr)))
    setkey(Prior, 'Gene')
  } else if (id == 'scat') {
    Prior <- readxl::read_excel(
      'data/raw/CK1CK2_2SS_MultiomeData4TRN_PeridermCells.xlsx',
      sheet='Prior_Allperiderm', .name_repair='minimal') |> as.data.table()
    names(Prior)[1] <- 'Gene'
    names(Prior)[-1] <- tolower(names(Prior)[-1])
    setkey(Prior, 'Gene')
  } else if (id == 'mall') {
    Prior <- fread('data/raw/100kbBinaryShieldStageMultiome_allEVL_danRer11.txt')
    names(Prior) <- c('Gene', tolower(names(Prior)[-1]))
    setkey(Prior, 'Gene')
  } else if (id == 'mspc') {
    Prior <- fread('data/raw/100kbBinaryShieldStageMultiome_onlyEVL_danRer11.txt')
    names(Prior) <- c('Gene', tolower(names(Prior)[-1]))
    setkey(Prior, 'Gene')
  } else if (id == 'm100') {
  Prior <- fread('data/raw/Binary_ShieldStageMultiome_EVLallProximal100Kb_TF2GeneNR.csv')
  names(Prior) <- c('Gene', tolower(names(Prior)[-1]))
  setkey(Prior, 'Gene')
  } else if (id == 'mc2') {
    Prior <- readxl::read_excel(
      'data/raw/InputDataAll_12July23_Sunil.xlsx',
      sheet='PriorMatrix_binary', .name_repair = 'minimal') |> as.data.table()
    names(Prior)[1] <- 'Gene'
    names(Prior)[-1] <- tolower(names(Prior)[-1])
    setkey(Prior, 'Gene')
  }

  # Match to data (pad missing with zeros)
  Prior <- pad(Prior[intersect(target, Gene)], setdiff(target, Prior$Gene))
  missing_tf <- setdiff(Data$tf, names(Prior))
  means <- Prior[, rowMeans(.SD), .SDcols = names(Prior)[-1]]
  if (length(missing_tf) > 0) {
    if (imputeTF == TRUE) {
      Prior[, (missing_tf) := means]
      } else{
        Prior[, (missing_tf) := 0]
      }
  }
  tf <- intersect(names(Prior)[-1], Data$tf) |> sort()
  Prior <- cbind(Prior[, 1], Prior[, ..tf])

  # Return
  setkey(Prior, 'Gene')
  attr(Prior, 'id') <- id
  class(Prior) <- c('prior', class(Prior))
  Prior
}

#' Helper for read_prior; pads empty rows in a prior with zeroes
#' @keywords internal
pad <- function(DT, new_genes) {
  if (length(new_genes) == 0) return(DT)
  New <- matrix(0, length(new_genes), ncol(DT)-1)
  colnames(New) <- names(DT)[-1]
  New <- as.data.table(New)
  New$Gene <- new_genes
  rbind(DT, New)[order(Gene)]
}

#' @export
print.prior <- function(x, ...) {
  cols <- 1:min(ncol(x), 11)
  data.table:::print.data.table(x[, ..cols])
  if (ncol(x) > 11) cat('(additional columns suppressed)\n')
  cat('\n')
  cat('id: ', attr(x, 'id'), '\n', sep='')
  cat('targets: ', nrow(x), '\n', sep='')
  cat('TFs: ', ncol(x)-1, '\n', sep='')
}
