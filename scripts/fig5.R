#' Create Figure 5 for paper: Illustrate modules / clustering for top network
net <- readRDS("data/network/wgcna/kna10-mc2-cor-s100.rds") |> transform()
make_graph(net) |>
  cluster_walktrap() |>
  sizes()

# a
pdf("scripts/out/fig5a.pdf", 5, 5)
cluster_plot(net)
dev.off()

# b (run scripts/panther-export.r first)
da <- parse_panther("data/go/cluster-a.json")
db <- parse_panther("data/go/cluster-b.json")
# dc <- parse_panther('data/go/cluster-c.json')
da$Cluster <- "1"
db$Cluster <- "2"
# dc$Cluster <- '3'
# go <- rbind(da, db, dc)
go <- rbind(da, db)
pdf("scripts/out/fig5b.pdf", 7, 4)
go[term.level == 0] |>
  _[plus_minus == "+"] |>
  _[str_detect(term.label, "RNA", negate = TRUE)] |>
  _[order(fdr), head(.SD, 5), Cluster] |>
  _[order(-Cluster, -fdr)] |>
  ggplot(aes(-log10(fdr), forcats::fct_inorder(term.label), fill = Cluster)) +
  geom_col() +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("")
dev.off()

# c
op <- options(digits = 2)
summary(net) |>
  _[order(-Centrality), head(.SD, 5), Cluster] |>
  _[order(Cluster, -Edgesum), .(Cluster, Gene, Centrality, Edgesum)] |>
  fwrite("scripts/out/fig5c.txt")
options(op)
