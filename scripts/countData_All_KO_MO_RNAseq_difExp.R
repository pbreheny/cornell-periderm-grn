#' Annika: What does this script do?
#' And can it have a better filename??

devtools::load_all()
library(DESeq2)
library(org.Dr.eg.db)
library(readxl)

# Setup -------------------------------------------------------------------

# Analyze Irf6 knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("WT", 4), rep("Mutant", 4))),
                      "batch" = factor(rep(c("A", "B", "C", "D"), 2)),
                      row.names = c("WT1", "WT2", "WT3", "WT4",
                                    "irf6_matMutant1", "irf6_matMutant2", "irf6_matMutant",
                                    "irf6_matMutant4"))
countData <- readxl::read_excel("data/raw/countData_All_KO_MO_RNAseq_corrected_7AUG23.xlsx", sheet = "irf6MatMutant_vs_WT", skip = 1)
M <- as.matrix(countData[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- countData$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
irf6_res <- results(dds, name = "condition_WT_vs_Mutant")

# Analyze tfap2a knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("ControlMO1", "ControlMO2", "ControlMO3",
                                    "tfap2aMO1", "tfap2aMO2", "tfap2aMO3"))
tfap2_dat <- readxl::read_excel("data/raw/countData_All_KO_MO_RNAseq_corrected_7AUG23.xlsx", skip = 1, sheet = "tfap2aMO_vs_controlMO")
M <- as.matrix(tfap2_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- tfap2_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
tfap2a_res <- results(dds, name = "condition_MO_vs_CT")

# Analyze grhl3 knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("stdControlMO4", "stdControlMO5", "stdControlMO6",
                                    "grhl3MO1", "grhl3MO2", "grhl3MO3"))
grhl3_dat <- readxl::read_excel("data/raw/countData_All_KO_MO_RNAseq_corrected_7AUG23.xlsx", skip = 1, sheet = "grhl3MO_vs_controlMO")
M <- as.matrix(grhl3_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- grhl3_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
grhl3_res <- results(dds, name = "condition_MO_vs_CT")

###########
gkktz_dat <- readxl::read_excel("data/raw/countData_All_KO_MO_RNAseq_corrected_7AUG23.xlsx", skip = 1, sheet = "grhl1-klf17-klf6a-tead3b-zbtb10")
###########
# Analyze grhl1 knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("ControlMO1", "ControlMO2", "ControlMO3",
                                    "grhl1MO1", "grhl1MO2", "grhl1MO3"))
grhl1_dat <- gkktz_dat[, 1:7]
M <- as.matrix(grhl1_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- grhl1_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
grhl1_res <- results(dds, name = "condition_MO_vs_CT")

# Analyze klf6a knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("ControlMO1", "ControlMO2", "ControlMO3",
                                    "klf6aMO1", "klf6aMO2", "klf6aMO3"))
klf6a_dat <- gkktz_dat[, c(1, 2, 3, 4, 8, 9, 10)]
M <- as.matrix(klf6a_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- klf6a_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
klf6a_res <- results(dds, name = "condition_MO_vs_CT")

# Analyze klf17 knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("ControlMO1", "ControlMO2", "ControlMO3",
                                    "klf17MO1", "klf17MO2", "klf17MO3"))
klf17_dat <- gkktz_dat[, c(1, 2, 3, 4, 11, 12, 13)]
M <- as.matrix(klf17_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- klf17_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
klf17_res <- results(dds, name = "condition_MO_vs_CT")

# Analyze tead3b knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("ControlMO1", "ControlMO2", "ControlMO3",
                                    "tead3bMO1", "tead3bMO2", "tead3bMO3"))
tead3b_dat <- gkktz_dat[, c(1, 2, 3, 4, 14, 15, 16)]
M <- as.matrix(tead3b_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- tead3b_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
tead3b_res <- results(dds, name = "condition_MO_vs_CT")

# Analyze zbtb10 knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("ControlMO1", "ControlMO2", "ControlMO3",
                                    "zbtb10MO1", "zbtb10MO2", "zbtb10MO3"))
zbtb10_dat <- gkktz_dat[, c(1, 2, 3, 4, 17, 18, 19)]
M <- as.matrix(zbtb10_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- zbtb10_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
zbtb10_res <- results(dds, name = "condition_MO_vs_CT")

###########
# cebpbMO-gata3MO-znf750_vs_CtlMO
cgz_dat <- readxl::read_excel("data/raw/countData_All_KO_MO_RNAseq_corrected_7AUG23.xlsx", skip = 1, sheet = "cebpbMO-gata3MO-znf750_vs_CtlMO")
###########

# Analyze cebpb knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("stdControlMO1", "stdControlMO2", "stdControlMO3",
                                    "cebpbMO1", "cebpbMO2", "cebpbMO3"))
cebpb_dat <- cgz_dat[, c(1, 2, 3, 4, 5, 6, 7)]
M <- as.matrix(cebpb_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- cebpb_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
cebpb_res <- results(dds, name = "condition_MO_vs_CT")

# Analyze gata3 knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("stdControlMO1", "stdControlMO2", "stdControlMO3",
                                    "gata3MO1", "gata3MO2", "gata3MO3"))
gata3_dat <- cgz_dat[, c(1, 2, 3, 4, 8, 9, 10)]
M <- as.matrix(gata3_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- gata3_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
gata3_res <- results(dds, name = "condition_MO_vs_CT")

# Analyze znf750 knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("CT", 3), rep("MO", 3))),
                      "batch" = factor(rep(c("A", "B", "C"), 2)),
                      row.names = c("stdControlMO1", "stdControlMO2", "stdControlMO3",
                                    "znf750MO1", "znf750MO2", "znf750MO3"))
znf750_dat <- cgz_dat[, c(1, 2, 3, 4, 11, 12, 13)]
M <- as.matrix(znf750_dat[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- znf750_dat$Gene
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
znf750_res <- results(dds, name = "condition_MO_vs_CT")

# Make Gold Standards -------------------------------------------------------------------
# Target Genes
# targets <- fread("data/zfTargetGenes_S2982.txt", header = FALSE)
#targets <- fread("data/targetGenes4283.txt", header = FALSE)

#targets <- readxl::read_excel("data/raw/InputDataAll_12July23_Sunil.xlsx", sheet = "TargetGenes+TFs")
# Target Genes (combine all data sources)
#fna <- import_expr('farrell', 'none', 'all')
#fns <- import_expr('farrell', 'none', 'shield')
kna <- import_expr('kenny', 'none', 'all')

targets <- kna$Exp$Gene

### IRF6 ###
irf6 <- stage_KO(DESeq2_result = irf6_res, TF = 'irf6', target_genes = targets,
                 restrict = TRUE)

saveRDS(irf6, 'data/knockout/irf6_KO.rds')

### TFAP2a ###
tfap2a <- stage_KO(DESeq2_result = tfap2a_res, TF = 'tfap2a', target_genes = targets,
                   restrict = TRUE)
saveRDS(tfap2a, 'data/knockout/tfap2a_KO.rds')

### GRHL3 ###
grhl3_gs <- stage_KO(DESeq2_result = grhl3_res, TF = 'grhl3', target_genes = targets,
                     restrict = TRUE)
saveRDS(grhl3_gs, 'data/knockout/grhl3_KO.rds')

### GRHL1 ###
grhl1_gs <- stage_KO(DESeq2_result = grhl1_res, TF = 'grhl1', target_genes = targets,
                    restrict = TRUE)
saveRDS(grhl1_gs, 'data/knockout/grhl1_KO.rds')


### KLF17 ###
klf17_gs <- stage_KO(DESeq2_result = irf6_res, TF = 'klf17', target_genes = targets,
                    restrict = TRUE)
saveRDS(klf17_gs, 'data/knockout/klf17_KO.rds')

### KLF6a ###
klf6a_gs <- stage_KO(DESeq2_result = klf6a_res, TF = 'klf6a', target_genes = targets,
                    restrict = TRUE)
saveRDS(klf6a_gs, 'data/knockout/klf6a_KO.rds')

### TEAD3b ###
tead3b_gs <- stage_KO(DESeq2_result = tead3b_res, TF = 'tead3b', target_genes = targets,
                    restrict = TRUE)
saveRDS(tead3b_gs, 'data/knockout/tead3b_KO.rds')

### IRF6 ###
zbtb10_gs <- stage_KO(DESeq2_result = zbtb10_res, TF = 'zbtb10', target_genes = targets,
                    restrict = TRUE)
saveRDS(zbtb10_gs, 'data/knockout/zbtb10_KO.rds')

### IRF6 ###
cebpb_gs <- stage_KO(DESeq2_result = cebpb_res, TF = 'cebpb', target_genes = targets,
                    restrict = TRUE)
saveRDS(cebpb_gs, 'data/knockout/cebpb_KO.rds')

### GATA3 ###
gata3_gs <- stage_KO(DESeq2_result = gata3_res, TF = 'gata3', target_genes = targets,
                    restrict = TRUE)
saveRDS(gata3_gs, 'data/knockout/gata3_KO.rds')

### ZNF750 ###
znf750_gs <- stage_KO(DESeq2_result = znf750_res, TF = 'znf750', target_genes = targets,
                    restrict = TRUE)
saveRDS(znf750_gs, 'data/knockout/znf750_KO.rds')

















