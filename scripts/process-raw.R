# Tissue
DT1 <- fread('data/EarlyEVL.txt', header=FALSE)
DT1[, Tissue := 'EarlyEVL']
DT2 <- fread('data/EVL.txt', header=FALSE)
DT2[, Tissue := 'EVL']
DT3 <- fread('data/LateEVL.txt', header=FALSE)
DT3[, Tissue := 'LateEVL']
DT4 <- fread('data/Epidermis.txt', header=FALSE)
DT4[, Tissue := 'Epidermis']
Tissue <- rbind(DT1, DT2, DT3, DT4)
setnames(Tissue, 'V1', 'Sample')
setkey(Tissue, 'Sample')

# Merge with time
Time <- readxl::read_excel('data/Time_series_data4GeneExpression.xlsx') |> as.data.table
setnames(Time, 'Sample Name', 'Sample')
setnames(Time, 'TimePoint', 'Time')
Tissue <- merge(Tissue, Time, by='Sample')
Tissue[, Sample := tolower(Sample)]
setkey(Tissue, 'Sample')

# TPM
TPM <- Read('data/evlSC_rawExpressionMatrix.txt')
target <- fread('data/targetGenes4283.txt', header=FALSE)[[1]]
tf <- names(import_prior('atac'))[-1]

# Fix 'sep' genes
TPM[, Gene := str_replace(Gene, '^sep(\\d+)$', 'sept\\1')]

# Integrity check
sum(!(target %in% TPM$Gene))
setdiff(target, TPM$Gene)
target <- intersect(union(target, tf), TPM$Gene)
tf <- intersect(tf, TPM$Gene)
dim(Tissue)
dim(TPM)
any(duplicated(Tissue$Gene))

# Save
target <- sort(target)
setkey(TPM, 'Gene')
saveRDS(list(Tissue=Tissue, TPM=TPM[target], target=target, tf=sort(tf)), 'rds/data-raw.rds')
