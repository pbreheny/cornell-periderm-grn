net <- readRDS('network/wgcna/kna10-cor-s100.rds')
s1 <- summary(net)
s2 <- summary(net, tf_only=FALSE)
top_25 <- s1[, .SD[order(-Centrality)][1:25], Cluster]
cluster_score <- function(i) {
  g <- top_25[Cluster == i]$Gene
  apply(abs(soft(net[,colnames(net) %in% g], 0.1)), 1, sum)
}
target_scores <- vapply(1:3, cluster_score, numeric(nrow(net)))
colnames(target_scores) <- paste0('m', 1:3)
dat <- merge(s2, as.data.table(target_scores, keep.rownames='Gene'), by='Gene')
wilcox.test(m1 ~ Panel, dat)
wilcox.test(m2 ~ Panel, dat)
wilcox.test(m3 ~ Panel, dat)

melt(dat[, .(Gene, Panel, m1, m2, m3)], measure.vars=c('m1', 'm2', 'm3')) |>
  ggplot(aes(Panel, value)) +
  geom_boxplot() +
  facet_grid(~variable) +
  ylab('Total connections to hub genes')
