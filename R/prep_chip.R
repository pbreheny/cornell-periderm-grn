#' Prepare ChIP-Seq data to make gold standard
#'
#' @param chip_data       path to ChIP-Seq data (or consensus peaks) of interest
#' @param TF              the transcription factor the ChIP-Seq experiment is for
#'
#' @return A bed file with properly formatted ChIP-Seq peaks
#'
#' @export
#'
#' @examples
#' prep_chip("data/raw/grhl3_consensus.bed", "grhl3")
#' prep_chip("data/raw/tfap2a_consensus.bed", "tfap2a")
#' prep_chip("data/raw/IRF6ChIP_danRer11_peaks.narrowPeak", "irf6")
#' @export

prep_chip <- function(chip_data, TF){
  chr_alias <- readRDS('data/raw/chr_alias.rds')
  chip <- fread(chip_data, select = c(1, 2, 3))
  names(chip) <- c('seqnames', 'start', 'end')

  for(i in 1:nrow(chip)){
    chip[i, chr := chr_alias[which(chr_alias$`UCSC database: danRer11` == chip$seqnames[i]), ]$sequenceName]
  }

  setcolorder(chip, c("chr", "start", "end", "seqnames"))
  write.table(chip, file=paste0("data/raw/", tolower(TF), "_chip.bed"), quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = FALSE)
}
