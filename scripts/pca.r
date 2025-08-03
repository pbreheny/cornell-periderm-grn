Data <- readRDS('rds/data.rds')
X <- as.matrix(Data$TPM[,-1])
rownames(X) <- Data$TPM$Gene
breheny::PCAplot(t(X), grp=Data$Tissue$Tissue)
