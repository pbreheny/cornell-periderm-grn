P1 <- import_prior('all')
dim(P1)
sum(!(P1$Gene %in% TPM$Gene))
sum(!(TPM$Gene %in% P1$Gene))
length(intersect(TPM$Gene, P1$Gene))
# (nrow(P1) + nrow(TPM) - sum(!(P1$Gene %in% TPM$Gene)) - sum(!(TPM$Gene %in% P1$Gene)))/2
length(intersect(TPM$Gene, names(P1)[-1]))

P2 <- import_prior('h3k')
dim(P2)
sum(!(P2$Gene %in% TPM$Gene))
sum(!(TPM$Gene %in% P2$Gene))
length(intersect(TPM$Gene, P2$Gene))
# (nrow(P2) + nrow(TPM) - sum(!(P2$Gene %in% TPM$Gene)) - sum(!(TPM$Gene %in% P2$Gene)))/2
length(intersect(TPM$Gene, names(P2)[-1]))
