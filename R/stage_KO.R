#' Match output from 'DESeq2' Data (KO experiment) with target genes
#'
#' @param DESeq2_result       Result from DESeq2 after calling 'DESeq2::results'
#' @param TF                  Transcription Factor you are building a gold standard for
#' @param target_genes        List of the target genes of interest
#' @param restrict_targets    Restrict the analysis to target genes? Default: TRUE
#'
#' @return A data.table of class 'gold standard'
#'
#' @examples
#' colData <- data.frame("condition" = factor(c(rep("WT", 4), rep("Mutant", 4))),
#' "batch" = factor(rep(c("A", "B", "C", "D"), 2)),
#' row.names = c("WT1", "WT2", "WT3", "WT4",
#'            "IRF6_Mut1", "IRF6_Mut2", "IRF6_Mut3",
#'            "IRF6_Mut4"))
#' countData <- fread("data/raw/MatrixSTAR_counts.txt", skip=2)
#' M <- as.matrix(countData[,-1])
#' colnames(M) <- rownames(colData)
#' rownames(M) <- countData$V1
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData = M,
#'                            colData = colData,
#'                            design = ~ condition + batch)
#' dds <- DESeq2::DESeq(dds, quiet=TRUE)
#' irf6_res <- DESeq2::results(dds, name = "condition_WT_vs_Mutant")
#' irf6_res$ENSEMBL <- rownames(irf6_res)
#' Data <- import_expr()
#' targets <- Data$Exp$Gene
#' irf6_gs <- stage_KO(DESeq2_result = irf6_res, TF = 'irf6', target_genes = targets, restrict = TRUE)
#' @export

stage_KO <- function(DESeq2_result, TF, target_genes, restrict_targets = TRUE){
  ash <- ashr::ash(betahat = DESeq2_result$log2FoldChange,
                   sebetahat = DESeq2_result$lfcSE)
  if(rownames(DESeq2_result)[1] != "1"){
    rownames(ash$result) <- rownames(DESeq2_result)
    KO_ash <- ash$result

    if(!is.null(DESeq2_result$ENSEMBL)){
      # Convert ensemble IDs to gene symbols
      df <-  suppressMessages(AnnotationDbi::select(
        org.Dr.eg.db::org.Dr.eg.db,
        keys = as.character(rownames(KO_ash)),
        columns = c('SYMBOL', 'ENSEMBL'),
        keytype = 'ENSEMBL')) |>
        data.table()

      # take only the first instance of an ENSEMBL ID
      unique_names <- df[!(duplicated(ENSEMBL))]

      # add gene symbol info into the GS_as df
      KO_ash <- data.table(cbind(unique_names, KO_ash))
    }else{
      KO_ash$SYMBOL <- rownames(KO_ash)
      KO_ash <- as.data.table(KO_ash)
    }
  }
  else{
    KO_ash <- ash$result
    if(!is.null(DESeq2_result$ENSEMBL)){
      KO_ash$ENSEMBL <- DESeq2_result$ENSEMBL
      KO_ash$SYMBOL <- DESeq2_result$SYMBOL
    }
    KO_ash <- as.data.table(KO_ash)
  }

  if(restrict_targets){
    KO_ash <- KO_ash[SYMBOL %in% as.character(target_genes)]
  }

  # make sure not to include duplicated gene entries
  KO_ash <- KO_ash[!duplicated(SYMBOL)]

  # change lfsr to 1 if betahat is missing
  KO_ash[is.na(betahat), lfsr := 1]

  ### add attributes and class to gold standard object ###
  attr(KO_ash, 'TF') <- TF
  class(KO_ash) <- c('KO', class(KO_ash))

  return(KO_ash)
}
