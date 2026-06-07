# Prepare raw data for inferelator
Data <- import_expr('kenny', 'none', 'all', 0.1)
Prior <- import_prior('m100', Data)

fwrite(Data$Exp, 'network/inferelator/kna10.tsv', sep='\t')
writeLines(Data$tf, 'network/inferelator/tf.txt')
fwrite(Prior, 'network/inferelator/m100.tsv', sep='\t')

# Run python script
system('python network/inferelator/run.py')

# Process inferelator output
Out <- fread('network/inferelator/out/random_seed_1/network.tsv.gz')
X <- matrix(0, length(Data$Exp$Gene), length(Data$tf), dimnames=list(Data$Exp$Gene, Data$tf))
X[DT$target, DT$regulator] <- DT$combined_confidences
attr(X, 'method') <- 'inferelator'
attr(X, 'samples') <- colnames(Data$Exp)
attr(X, 'prior') <- attr(Prior, 'id')
attr(X, 'strength') <- 0
attr(X, 'alpha') <- 0
class(X) <- c('network', class(X))
saveRDS(X, 'network/inferelator/kna10-bbsr-s100.rds')
