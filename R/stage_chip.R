#' Match peaks from a ChIP-Seq experiment with target genes
#' 
#' Requires bedtools to be installed and capable of being called via `system()`
#' 
#' @param target_bed      File path to a bedfile of target genes
#' @param chip_bed        File path to a bedfile of ChIP-Seq peaks 
#' @param basePairs       The number of overlapping base pairs required to associate a ChIP-Seq peak with a target gene
#' 
#' @returns   A data.table of class 'chipSeq'
#' 
#' @examples
#' \dontrun{
#' grhl3_chip <- stage_chip(target_bed = "data/raw/target_genes.bed", chip_bed = "data/raw/grhl3_chip.bed")
#' tfap2a_chip <- stage_chip(target_bed = "data/raw/target_genes.bed", chip_bed = "data/raw/tfap2a_chip.bed")
#' irf6_chip <- stage_chip(target_bed = "data/raw/target_genes.bed", chip_bed = "data/raw/irf6_chip.bed")
#' }
#' @export

stage_chip <- function(target_bed, chip_bed, basePairs = 1){
  # sort
  system(paste0("bedtools sort -i ", target_bed, " > ", 
                unlist(strsplit(target_bed, ".bed")), "_sort.bed"))
  system(paste0("bedtools sort -i ", chip_bed, " > ", 
                unlist(strsplit(chip_bed, ".bed")), "_sort.bed"))
  
  # intersect
  system(paste0("bedtools intersect -a ", unlist(strsplit(target_bed, ".bed")), 
                "_sort.bed -b ", unlist(strsplit(chip_bed, ".bed")), "_sort.bed -wao > ",
                unlist(strsplit(chip_bed, "_chip.bed")), "_int_targets.out"))
  
  int <- fread(paste0(unlist(strsplit(chip_bed, "_chip.bed")), "_int_targets.out"))
  names(int) <- c("chr_targ", "start_targ", "end_targ", "ENSEMBL", "chr_chip", 
                  "start_chip", "end_chip", "seq", "intersection")
  int[, SYMBOL := AnnotationDbi::mapIds(org.Dr.eg.db, keys = ENSEMBL,
                                        column = c('SYMBOL'), keytype = 'ENSEMBL')]
  int[, True_Pos := ifelse(int$intersection > (basePairs - 1), 1, 0), ]
  dat <- int[match(unique(int$ENSEMBL), int$ENSEMBL), c(1, 2, 3, 4, 10, 11)]
  
  ### add attributes and class to chip-seq object ###
  attr(dat, 'TF') <- strsplit(str_sub(chip_bed, 10), '_')[[1]][1]
  class(dat) <- c('chipSeq', class(dat))
  
  return(dat)
}
