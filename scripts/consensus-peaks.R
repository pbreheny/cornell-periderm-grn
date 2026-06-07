#' Annika: What does this script do?

source("scripts/setup.R")
library(rtracklayer)
library(GenomicRanges)
library(org.Dr.eg.db)

# Find Consensus Peaks from ChIP-Seq and CUT&RUN data --------------------------

grhl3_files <- list.files('data', pattern = '^MACS2_callpeak_Grhl3*', full.names = TRUE)
tfap2a_files <- list.files('data', pattern = 'MACS2_callpeak_Tfap2a*', full.names = TRUE)

# grhl3 -----------------------------------------------------------------------
grhl3_granges <- lapply(grhl3_files, import, format = "narrowPeak")
grhl3_grangeslist <- GRangesList(grhl3_granges)
grhl3_coverage <- coverage(grhl3_grangeslist)
covered_ranges_grhl3 <- slice(grhl3_coverage, lower = 2, rangesOnly=T)
covered_ranges_grhl3 <- GRanges(covered_ranges_grhl3)
export(covered_ranges_grhl3, "data/grhl3_consensus.bed")

# tfap2a -----------------------------------------------------------------------
tfap2a_granges <- lapply(tfap2a_files, import, format = "narrowPeak")
tfap2a_grangeslist <- GRangesList(tfap2a_granges)
tfap2a_coverage <- coverage(tfap2a_grangeslist)
covered_ranges_tfap2a <- slice(tfap2a_coverage, lower = 2, rangesOnly=T)
covered_ranges_tfap2a <- GRanges(covered_ranges_tfap2a)
export(covered_ranges_tfap2a, "data/tfap2a_consensus.bed")


# Convert provided chrom names to standard chrom names
# chromosome alias
chr_alias <- fread("data/danRer11.chromAlias.txt", fill = TRUE, header = TRUE)
for (i in 1:nrow(chr_alias)) {
  chr_alias$sequenceName[i] <- strsplit(chr_alias$sequenceName[[i]], '_')[[1]]
  if(chr_alias$V4[i] != ""){
    chr_alias$`alias names`[i] <- chr_alias$`UCSC database: danRer11`[i]
    chr_alias$`UCSC database: danRer11`[i] <- chr_alias$V4[i]
  }
}

# grhl3 -----------------------------------------------------------------------
prep_chip("data/grhl3_consensus.bed", "grhl3")

# tfap2a -----------------------------------------------------------------------
prep_chip("data/tfap2a_consensus.bed", "tfap2a")

# irf6 -----------------------------------------------------------------------
prep_chip("data/IRF6ChIP_danRer11_peaks.narrowPeak", "irf6")


library(biomaRt)
# Target Genes
Data <- readRDS('rds/data.rds')
targets <- Data$target
# convert gene symbols to ensebl IDs
ensembl_targets <- AnnotationDbi::mapIds(org.Dr.eg.db, keys = targets,
                      column = c('ENSEMBL'), keytype = 'SYMBOL')

## note that there are 4363 unique target gene symbols but only 3933 unique ensembl gene IDs ##

ensembl_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")

total_bed <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name',
                                   'start_position', 'end_position'), mart = ensembl_mart)

target_bed <- total_bed[which(total_bed$ensembl_gene_id %in% ensembl_targets), ]
setcolorder(target_bed, c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"))

target_bed$chromosome_name <- paste0("chr", target_bed$chromosome_name)
names(target_bed) <- c("chr", "start", "end", "ensembl_gene_id")

# export bedfile
write.table(target_bed, file = "data/target_genes.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



