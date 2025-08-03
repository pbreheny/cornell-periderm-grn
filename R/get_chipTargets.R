#' Add flanking regions to a list of target genes
#'
#' @param target_genes       bedfile of target genes
#' @param flank_size         number of base pairs that should flank either side of the gene start site
#'
#' @examples
#' Data <- import_expr()
#' get_chipTargets(target_genes = Data$Exp$Gene, flank_size = 100)
#' @export

get_chipTargets <- function(target_genes, flank_size = 0){
  options(scipen = 9999)
  # convert gene symbols to ensebl IDs
  suppressMessages(ensembl_targets <- AnnotationDbi::mapIds(
    org.Dr.eg.db::org.Dr.eg.db,
    keys = target_genes,
    column = 'ENSEMBL',
    keytype = 'SYMBOL'))

  ensembl_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")
  total_bed <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'chromosome_name',
                                             'start_position', 'end_position'), mart = ensembl_mart)
  total_bed <- total_bed[total_bed$chromosome_name != "MT", ]

  target_bed <- total_bed[which(total_bed$ensembl_gene_id %in% ensembl_targets), ]
  setcolorder(target_bed, c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"))
  target_bed$chromosome_name <- paste0("chr", target_bed$chromosome_name)
  names(target_bed) <- c("chr", "start", "end", "ensembl_gene_id")

  N <- length(unique(target_bed$chr))
  n <- N - 25
  for (i in 1:n){
    target_bed[which(target_bed$chr %in% unique(target_bed$chr)[26:N]), ]$chr[i] <-
      paste0("chrUn_", unlist(strsplit(target_bed[which(target_bed$chr %in% unique(target_bed$chr)[26:N]), ]$chr[i], split = "[.]|chr"))[2], "v1")
  }

  # export bedfiles
  target_bed[order(target_bed$chr, target_bed$start), ] %>%
    write.table(file = "data/raw/target_genes.bed",
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  if (flank_size > 0) {
    system(paste0("bedtools flank -i data/raw/target_genes", ".bed -g data/raw/danRer11.chrom.sizes -b ",
                  flank_size, " > data/raw/targets_flank", flank_size, ".bed"))
  }
}
