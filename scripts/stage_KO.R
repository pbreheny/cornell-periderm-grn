source("scripts/setup.R")
library(DESeq2)
library(org.Dr.eg.db)
library(readxl)

# Setup -------------------------------------------------------------------

# Analyze Irf6 knockout experiment with DESeq2
colData <- data.frame("condition" = factor(c(rep("WT", 4), rep("Mutant", 4))),
                      "batch" = factor(rep(c("A", "B", "C", "D"), 2)),
                      row.names = c("WT1", "WT2", "WT3", "WT4",
                                    "IRF6_Mut1", "IRF6_Mut2", "IRF6_Mut3",
                                    "IRF6_Mut4"))
countData <- fread("./data/raw/MatrixSTAR_counts.txt", skip=2)
M <- as.matrix(countData[,-1])
colnames(M) <- rownames(colData)
rownames(M) <- countData$V1
dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = colData,
                              design = ~ condition + batch)
dds <- DESeq(dds)
irf6_res <- results(dds, name = "condition_WT_vs_Mutant")
irf6_res$ENSEMBL <- rownames(irf6_res)

# irf6_norm <- lfcShrink(dds, coef="condition_WT_vs_Mutant", type="normal")
# irf6_ashr <- lfcShrink(dds, coef = "condition_WT_vs_Mutant", type="ashr", returnList = TRUE)
# irf6_ashr2 <- ashr::ash(betahat = irf6_res$log2FoldChange,
#                        sebetahat = irf6_res$lfcSE)
# irf6_apeglm <- lfcShrink(dds, coef = "condition_WT_vs_Mutant", type="apeglm", returnList = TRUE)

# par(mfcol = c(1, 2))
# irf6_apeglm$res$log2FoldChange |> hist()
# irf6_ashr$res$log2FoldChange |> hist()

# irf6_apeglm$res$log2FoldChange |> boxplot()
# irf6_ashr$res$log2FoldChange |> boxplot()

# Prepare the TFAP2a experiment

tfap2_dat <- readxl::read_excel("data/raw/Dooley_tfap2a_vs_WT_deTCT.xlsx")
names(tfap2_dat) <- c("ENSEMBL", "SYMBOL", "log2FoldChange", "pvalue", "padj")
tfap2_dat$lfcSE <- -(abs(tfap2_dat$log2FoldChange) / qnorm(0.5 * tfap2_dat$pvalue))
tfap2a_res <- as.data.frame(tfap2_dat)

# Prepare the TFAP2a2c experiment

tfap2ac_dat <- readxl::read_excel("data/raw/Dooley_tfap2a2c_deTCT.xlsx",
                                  sheet = "4som_tfap2a2c_vs_WT",
                                  col_types = c("text", "text", "numeric", "numeric", "numeric"))
names(tfap2ac_dat) <- c("ENSEMBL", "SYMBOL", "log2FoldChange", "pvalue", "padj")
tfap2ac_dat$lfcSE <- -(abs(tfap2ac_dat$log2FoldChange) / qnorm(0.5 * tfap2ac_dat$pvalue))
tfap2a2c_res <- as.data.frame(tfap2ac_dat)

# prepare the TFAP2a MO experiment

colData_tfap2aMO <- data.frame("condition" = factor(c(rep("Control", 3), rep("Mutant", 3))),
                            "batch" = factor(rep(c("A", "B", "C"), 2)),
                            row.names = c("Control1", "Control2", "Control3",
                                          "tfap2aMO_Mut1", "tfap2aMO_Mut2", "tfap2aMO_Mut3"))
countData_tfap2aMO <- readxl::read_excel("data/raw/tfap2aMOvsControlMO_RNAseq.xlsx")
M_tfap2aMO <- as.matrix(countData_tfap2aMO[, -1])
colnames(M_tfap2aMO) <- rownames(colData_tfap2aMO)
rownames(M_tfap2aMO) <- countData_tfap2aMO$Geneid
dds_tfap2aMO <- DESeqDataSetFromMatrix(countData = M_tfap2aMO,
                                    colData = colData_tfap2aMO,
                                    design = ~ condition + batch)
dds_tfap2aMO <- DESeq(dds_tfap2aMO)
resultsNames(dds_tfap2aMO)
tfap2aMO_res <- results(dds_tfap2aMO, name = "condition_Mutant_vs_Control")

# Prepare the KLF17 experiment

klf17_dat <- read_excel("data/raw/KLF17MO_vs_ControlMO_MicroarrayResults.xlsx")
names(klf17_dat) <- c("SYMBOL", "log2FoldChange", "pvalue", "padj")
klf17_dat$lfcSE <- -(abs(klf17_dat$log2FoldChange) /
                       qnorm(0.5 * ifelse(klf17_dat$pvalue != 1,
                                          klf17_dat$pvalue,
                                          0.99)))
klf17_res <- as.data.frame(klf17_dat)
klf17_res$ENSEMBL <- mapIds(org.Dr.eg.db, keys = klf17_res$SYMBOL, keytype = "SYMBOL",
                            column = "ENSEMBL")

# no raw data for TFAP2a or TFAP2a2c?

# Prepare the KLF17 experiment
# microarraay

# Prepare the cebpb experiment
colData_cebpb <- data.frame("condition" = factor(c(rep("Control", 3), rep("Mutant", 3))),
                            "batch" = factor(rep(c("A", "B", "C"), 2)),
                            row.names = c("Control1", "Control2", "Control3",
                                          "cebpb_Mut1", "cebpb_Mut2", "cebpb_Mut3"))
countData_cebpb <- fread("data/raw/cebpbMOvsControlMO.csv", skip = 1)
M_cebpb <- as.matrix(countData_cebpb[, -1])
colnames(M_cebpb) <- rownames(colData_cebpb)
rownames(M_cebpb) <- countData_cebpb$V1
dds_cebpb <- DESeqDataSetFromMatrix(countData = M_cebpb,
                                    colData = colData_cebpb,
                                    design = ~ condition + batch)
dds_cebpb <- DESeq(dds_cebpb)
resultsNames(dds_cebpb)
cebpb_res <- results(dds_cebpb, name = "condition_Mutant_vs_Control")

# prepare the gata3 experiment
colData_gata3 <- data.frame("condition" = factor(c(rep("Control", 3), rep("Mutant", 3))),
                            "batch" = factor(rep(c("A", "B", "C"), 2)),
                            row.names = c("ControlMO1", "ControlMO2", "ControlMO3",
                                          "gata3M01", "gata3M02", "gata3M03"))
countData_gata3 <- fread("./data/raw/gata3MOvsControlMO.csv", skip = 1)
M_gata3 <- as.matrix(countData_gata3[, -1])
colnames(M_gata3) <- rownames(colData_gata3)
rownames(M_gata3) <- countData_gata3$V1
dds_gata3 <- DESeqDataSetFromMatrix(countData = M_gata3,
                                    colData = colData_gata3,
                                    design = ~ condition + batch)
dds_gata3 <- DESeq(dds_gata3)
resultsNames(dds_gata3)
gata3_res <- results(dds_gata3, name = "condition_Mutant_vs_Control")

# prepare the grhl3 experiment
colData_grhl3 <- data.frame("condition" = factor(c(rep("Control", 3), rep("Mutant", 3))),
                            "batch" = factor(rep(c("A", "B", "C"), 2)),
                            row.names = c("ControlMO1", "ControlMO2", "ControlMO3",
                                          "grhl3M01", "grhl3M02", "grhl3M03"))
countData_grhl3 <- fread("./data/raw/grhl3MOvsControlMO.csv", skip = 1)
M_grhl3 <- as.matrix(countData_grhl3[, -1])
colnames(M_grhl3) <- rownames(colData_grhl3)
rownames(M_grhl3) <- countData_grhl3$V1
dds_grhl3 <- DESeqDataSetFromMatrix(countData = M_grhl3,
                                    colData = colData_grhl3,
                                    design = ~ condition + batch)
dds_grhl3 <- DESeq(dds_grhl3)
resultsNames(dds_grhl3)
grhl3_res <- results(dds_grhl3, name = "condition_Mutant_vs_Control")

# prepare the znf750 experiment
colData_znf750 <- data.frame("condition" = factor(c(rep("Control", 3), rep("Mutant", 3))),
                             "batch" = factor(rep(c("A", "B", "C"), 2)),
                             row.names = c("ControlMO1", "ControlMO2", "ControlMO3",
                                           "Znf750M01", "Znf750M02", "Znf750M03"))
countData_znf750 <- fread("./data/raw/znf750MOvsControlMO.csv", skip = 1)
M_znf750 <- as.matrix(countData_znf750[, -1])
colnames(M_znf750) <- rownames(colData_znf750)
rownames(M_znf750) <- countData_znf750$V1
dds_znf750 <- DESeqDataSetFromMatrix(countData = M_znf750,
                                     colData = colData_znf750,
                                     design = ~ condition + batch)
dds_znf750 <- DESeq(dds_znf750)
resultsNames(dds_znf750)
znf750_res <- results(dds_znf750, name = "condition_Mutant_vs_Control")

# znf750_norm <- lfcShrink(dds, coef="condition_WT_vs_Mutant", type="normal")
# znf750_ashr <- lfcShrink(dds, coef = "condition_WT_vs_Mutant", type="ashr", returnList = TRUE)
# znf750_ashr2 <- ashr::ash(betahat = znf750_res$log2FoldChange,
#                        sebetahat = znf750_res$lfcSE)
# znf750_apeglm <- lfcShrink(dds, coef = "condition_WT_vs_Mutant", type="apeglm", returnList = TRUE)

# znf750_apeglm$res$log2FoldChange |> hist()
# znf750_ashr$res$log2FoldChange |> hist()

# znf750_apeglm$res$log2FoldChange |> boxplot()
# znf750_ashr$res$log2FoldChange |> boxplot()

# Make Gold Standards -------------------------------------------------------------------
# Target Genes
# targets <- fread("data/zfTargetGenes_S2982.txt", header = FALSE)
#targets <- fread("data/targetGenes4283.txt", header = FALSE)
#-----------------------------------------------------------------------------

# Target Genes (combine all data sources)
fna_w <- import_expr('farrell', 'none', 'all')
fns_w <- import_expr('farrell', 'none', 'shield')
kna_w <- import_expr('kenny', 'none', 'all')
fna_m <- import_expr('farrell', 'none', 'all', filter = 0.3)
fns_m <- import_expr('farrell', 'none', 'shield', filter = 0.3)
kna_m <- import_expr('kenny', 'none', 'all', filter = 0.3)
fna_s <- import_expr('farrell', 'none', 'all', filter = 0.5)
fns_s <- import_expr('farrell', 'none', 'shield', filter = 0.5)
kna_s <- import_expr('kenny', 'none', 'all', filter = 0.5)

# targets <- unique(c(fna$Exp$Gene, fns$Exp$Gene, kna$Exp$Gene))
weak_targets <- unique(c(fna_w$Exp$Gene, fns_w$Exp$Gene, kna_w$Exp$Gene))
med_targets <- unique(c(fna_m$Exp$Gene, fns_m$Exp$Gene, kna_m$Exp$Gene))
strong_targets <- unique(c(fna_s$Exp$Gene, fns_s$Exp$Gene, kna_s$Exp$Gene))

### IRF6 ###
irf6_w <- stage_KO(DESeq2_result = irf6_res, TF = 'irf6', target_genes = weak_targets,
                   restrict = TRUE)
irf6_m <- stage_KO(DESeq2_result = irf6_res, TF = 'irf6', target_genes = med_targets,
                   restrict = TRUE)
irf6_s <- stage_KO(DESeq2_result = irf6_res, TF = 'irf6', target_genes = strong_targets,
                   restrict = TRUE)
saveRDS(irf6_w, 'data/knockout/weak/irf6_KO.rds')
saveRDS(irf6_m, 'data/knockout/med/irf6_KO.rds')
saveRDS(irf6_s, 'data/knockout/strong/irf6_KO.rds')

irf6_gs <- stage_KO(DESeq2_result = irf6_res, TF = 'irf6', target_genes = targets,
                    restrict = TRUE)
saveRDS(irf6_gs, 'data/knockout/irf6_KO.rds')

### TFAP2a ###
tfap2a_gs <- stage_KO(DESeq2_result = tfap2a_res, TF = 'tfap2a', target_genes = targets,
                      restrict = TRUE)
saveRDS(tfap2a_gs, 'data/knockout/tfap2a_KO.rds')

### TFAP2a2c ###
tfap2a2c_gs <- stage_KO(DESeq2_result = tfap2a2c_res, TF = 'tfap2a2c', target_genes = targets,
                      restrict = TRUE)
saveRDS(tfap2a2c_gs, 'data/knockout/tfap2a2c_KO.rds')

### TFAP2aMO ###
tfap2aMO_gs <- stage_KO(DESeq2_result = tfap2aMO_res, TF = 'tfap2aMO', target_genes = targets,
                        restrict = TRUE)
saveRDS(tfap2aMO_gs, 'data/knockout/tfap2aMO_KO.rds')

### KLF17 ###
klf17_gs <- stage_KO(DESeq2_result = klf17_res, TF = 'klf17', target_genes = targets,
                     restrict = TRUE)
saveRDS(klf17_gs, 'data/knockout/klf17_KO.rds')

### CEBPB ###
cebpb_gs <- stage_KO(DESeq2_result = cebpb_res, TF = 'cebpb', target_genes = targets,
                     restrict = TRUE)
saveRDS(cebpb_gs, 'data/knockout/cebpb_KO.rds')

### GATA3 ###
gata3_gs <- stage_KO(DESeq2_result = gata3_res, TF = 'gata3', target_genes = targets,
                     restrict = TRUE)
saveRDS(gata3_gs, 'data/knockout/gata3_KO.rds')

### GRHL3 ###
grhl3_gs <- stage_KO(DESeq2_result = grhl3_res, TF = 'grhl3', target_genes = targets,
                     restrict = TRUE)
saveRDS(grhl3_gs, 'data/knockout/grhl3_KO.rds')

### ZNF750 ###
znf750_gs <- stage_KO(DESeq2_result = znf750_res, TF = 'znf750', target_genes = targets,
                     restrict = TRUE)
saveRDS(znf750_gs, 'data/knockout/znf750_KO.rds')

