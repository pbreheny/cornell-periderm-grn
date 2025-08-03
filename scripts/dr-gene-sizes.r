#' Calculate zebrafish gene sizes
#' Output: data/raw/dr-gene-sizes.rds (a data.table)

suppressPackageStartupMessages(library(org.Dr.eg.db))
suppressPackageStartupMessages(library(TxDb.Drerio.UCSC.danRer11.refGene))

# Get exon sizes (as a proxy for transcript size)
txdb <- TxDb.Drerio.UCSC.danRer11.refGene
GRList <- exonsBy(txdb, by = "gene")
w <- vapply(width(GRList), sum, 1)

# Manually add grhl3, since it is only a "predicted" gene
DT1 <- rbind(
  data.table(ENTREZID=names(w), size=w),
  data.table(ENTREZID='794613', size=2659)) # From https://zfin.org/ZDB-GENE-030131-9854

# Match to symbol
DT2 <- select(org.Dr.eg.db, keys=keys(org.Dr.eg.db, keytype='ENTREZID'), columns=c('SYMBOL', 'GENENAME'), keytype='ENTREZID') %>% as.data.table
DT2 <- DT2[, first(.SD), ENTREZID]
DT <- merge(DT1, DT2, by='ENTREZID')

# Save
setnames(DT, 'SYMBOL', 'Gene')
setcolorder(DT, 'Gene')
setkey(DT, 'Gene')
saveRDS(DT, 'data/raw/dr-gene-sizes.rds')
