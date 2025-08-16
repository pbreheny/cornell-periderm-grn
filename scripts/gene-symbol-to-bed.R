library(rtracklayer)
library(GenomicRanges)
library(org.Dr.eg.db)
library(biomaRt)

# Target Genes
gene_symbols <- readxl::read_xlsx('data/raw/InputDataAll_12July23_Sunil.xlsx', sheet = 'TargetGenes+TFs')
targets <- gene_symbols$TargetGenes

# convert gene symbols to ensebl IDs
ensembl_targets <- AnnotationDbi::mapIds(org.Dr.eg.db, keys = targets,
                                         column = c('ENSEMBL'), keytype = 'SYMBOL')

## note that there are 3787 unique target gene symbols but only 3520 unique ensembl gene IDs ##

ensembl_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")

total_bed <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name',
                                  'start_position', 'end_position'), mart = ensembl_mart)

target_bed <- total_bed[which(total_bed$ensembl_gene_id %in% ensembl_targets), ]
setcolorder(target_bed, c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"))

target_bed$chromosome_name <- paste0("chr", target_bed$chromosome_name)
names(target_bed) <- c("chr", "start", "end", "ensembl_gene_id")

# export bedfile
write.table(target_bed, file = "data/target_genes_jul.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
