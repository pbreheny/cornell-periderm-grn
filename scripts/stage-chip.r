devtools::load_all()
library(rtracklayer)
library(GenomicRanges)
library(org.Dr.eg.db)
library(biomaRt)

# Find Consensus Peaks from ChIP-Seq and CUT&RUN data --------------------------
grhl3_files <- list.files('data/raw', pattern = '^MACS2_callpeak_Grhl3*', full.names = TRUE)
tfap2a_files <- list.files('data/raw', pattern = 'MACS2_callpeak_Tfap2a*', full.names = TRUE)

# grhl3 -----------------------------------------------------------------------
grhl3_granges <- lapply(grhl3_files, import, format = "narrowPeak")
grhl3_grangeslist <- GRangesList(grhl3_granges)
grhl3_coverage <- coverage(grhl3_grangeslist)
covered_ranges_grhl3 <- IRanges::slice(grhl3_coverage, lower = 2, rangesOnly=T)
covered_ranges_grhl3 <- GRanges(covered_ranges_grhl3)
export(covered_ranges_grhl3, "data/raw/grhl3_consensus.bed")

# tfap2a -----------------------------------------------------------------------
tfap2a_granges <- lapply(tfap2a_files, import, format = "narrowPeak")
tfap2a_grangeslist <- GRangesList(tfap2a_granges)
tfap2a_coverage <- coverage(tfap2a_grangeslist)
covered_ranges_tfap2a <- IRanges::slice(tfap2a_coverage, lower = 2, rangesOnly=T)
covered_ranges_tfap2a <- GRanges(covered_ranges_tfap2a)
export(covered_ranges_tfap2a, "data/raw/tfap2a_consensus.bed")


# Convert provided chrom names to standard chrom names -------------------------
# chromosome alias
chr_alias <- fread("data/raw/danRer11.chromAlias_ucsc.txt", fill = TRUE, header = TRUE)

for (i in 1:nrow(chr_alias)) {
  chr_alias$sequenceName[i] <- strsplit(chr_alias$sequenceName[[i]], '_')[[1]]
  if(chr_alias$V4[i] != ""){
    chr_alias$`alias names`[i] <- chr_alias$`UCSC database: danRer11`[i]
    chr_alias$`UCSC database: danRer11`[i] <- chr_alias$V4[i]
  }
}
saveRDS(chr_alias, 'data/raw/chr_alias.rds')

# Target Genes (combine all data sources)
fna <- import_expr('farrell', 'none', 'all')
fns <- import_expr('farrell', 'none', 'shield')
kna <- import_expr('kenny', 'none', 'all')

targets <- unique(c(fna$Exp$Gene, fns$Exp$Gene, kna$Exp$Gene))
# convert gene symbols to ensebl IDs
ensembl_targets <- AnnotationDbi::mapIds(org.Dr.eg.db, keys = targets,
                                         column = c('ENSEMBL'), keytype = 'SYMBOL')

## note that there are 4343 unique target gene symbols but only 3901 unique ensembl gene IDs ##

ensembl_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")
total_bed <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name',
                                  'transcription_start_site', 'transcript_end'), mart = ensembl_mart)
target_bed <- total_bed[which(total_bed$ensembl_gene_id %in% ensembl_targets), ]
setcolorder(target_bed, c("chromosome_name", "transcription_start_site", "transcript_end", "ensembl_gene_id"))
target_bed$chromosome_name <- paste0("chr", target_bed$chromosome_name)
names(target_bed) <- c("chr", "start", "end", "ensembl_gene_id")

# export bedfile
write.table(target_bed, file = "data/raw/target_genes.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Make target gene lists with flanking regions of 5kb, 10kb, and 100kb
## 1kb
get_chipTargets(target_genes = targets, flank_size = 1000)
## 5kb
get_chipTargets(target_genes = targets, flank_size = 5000)
## 10kb
get_chipTargets(target_genes = targets, flank_size = 10000)
## 100kb
get_chipTargets(target_genes = targets, flank_size = 100000)

# Make Gold Standards ----------------------------------------------------- ###

# grhl3 -----------------------------------------------------------------------
prep_chip("data/raw/grhl3_consensus.bed", "grhl3")
# no flank
grhl3_chip <- stage_chip(target_bed = "data/raw/target_genes.bed", chip_bed = "data/raw/grhl3_chip.bed")
attr(grhl3_chip, "Flank") <- '0kb'
saveRDS(grhl3_chip, 'data/chip/grhl3_chip.rds')
# 5kb flank
grhl3_chip_5kb <- stage_chip(target_bed = "data/raw/targets_flank5000.bed", chip_bed = "data/raw/grhl3_chip.bed")
attr(grhl3_chip_5kb, "Flank") <- '5kb'
saveRDS(grhl3_chip_5kb, 'data/chip/grhl3_chip_5kb.rds')
# 10kb flank
grhl3_chip_10kb <- stage_chip(target_bed = "data/raw/targets_flank10000.bed", chip_bed = "data/raw/grhl3_chip.bed")
attr(grhl3_chip_10kb, "Flank") <- '10kb'
saveRDS(grhl3_chip_10kb, 'data/chip/grhl3_chip_10kb.rds')
# 100kb flank
grhl3_chip_100kb <- stage_chip(target_bed = "data/raw/targets_flank100000.bed", chip_bed = "data/raw/grhl3_chip.bed")
attr(grhl3_chip_100kb, "Flank") <- '100kb'
saveRDS(grhl3_chip_100kb, 'data/chip/grhl3_chip_100kb.rds')


# tfap2a -----------------------------------------------------------------------
prep_chip("data/raw/tfap2a_consensus.bed", "tfap2a")
# no flank
tfap2a_chip <- stage_chip(target_bed = "data/raw/target_genes.bed", chip_bed = "data/raw/tfap2a_chip.bed")
attr(tfap2a_chip, "Flank") <- '0kb'
saveRDS(tfap2a_chip, 'data/chip/tfap2a_chip.rds')
# 5kb flank
tfap2a_chip_5kb <- stage_chip(target_bed = "data/raw/targets_flank5000.bed", chip_bed = "data/raw/tfap2a_chip.bed")
attr(tfap2a_chip_5kb, "Flank") <- '5kb'
saveRDS(tfap2a_chip_5kb, 'data/chip/tfap2a_chip_5kb.rds')
# 10kb flank
tfap2a_chip_10kb <- stage_chip(target_bed = "data/raw/targets_flank10000.bed", chip_bed = "data/raw/tfap2a_chip.bed")
attr(tfap2a_chip_10kb, "Flank") <- '10kb'
saveRDS(tfap2a_chip_10kb, 'data/chip/tfap2a_chip_10kb.rds')
# 100kb flank
tfap2a_chip_100kb <- stage_chip(target_bed = "data/raw/targets_flank100000.bed", chip_bed = "data/raw/tfap2a_chip.bed")
attr(tfap2a_chip_100kb, "Flank") <- '100kb'
saveRDS(tfap2a_chip_100kb, 'data/chip/tfap2a_chip_100kb.rds')


# irf6 -----------------------------------------------------------------------
prep_chip("data/raw/IRF6ChIP_danRer11_peaks.narrowPeak", "irf6")
# no flank
irf6_chip <- stage_chip(target_bed = "data/raw/target_genes.bed", chip_bed = "data/raw/irf6_chip.bed")
attr(irf6_chip, "Flank") <- '0kb'
saveRDS(irf6_chip, 'data/chip/irf6_chip.rds')
# 1kb flank
irf6_chip_1kb <- stage_chip(target_bed = "data/raw/targets_flank1000.bed", chip_bed = "data/raw/irf6_chip.bed")
attr(irf6_chip_1kb, "Flank") <- '1kb'
saveRDS(irf6_chip_1kb, 'data/chip/irf6_chip_1kb.rds')
# 5kb flank
irf6_chip_5kb <- stage_chip(target_bed = "data/raw/targets_flank5000.bed", chip_bed = "data/raw/irf6_chip.bed")
attr(irf6_chip_5kb, "Flank") <- '5kb'
saveRDS(irf6_chip_5kb, 'data/chip/irf6_chip_5kb.rds')
# 10kb flank
irf6_chip_10kb <- stage_chip(target_bed = "data/raw/targets_flank10000.bed", chip_bed = "data/raw/irf6_chip.bed")
attr(irf6_chip_10kb, "Flank") <- '10kb'
saveRDS(irf6_chip_10kb, 'data/chip/irf6_chip_10kb.rds')
# 100kb flank
irf6_chip_100kb <- stage_chip(target_bed = "data/raw/targets_flank100000.bed", chip_bed = "data/raw/irf6_chip.bed")
attr(irf6_chip_100kb, "Flank") <- '100kb'
saveRDS(irf6_chip_100kb, 'data/chip/irf6_chip_100kb.rds')


