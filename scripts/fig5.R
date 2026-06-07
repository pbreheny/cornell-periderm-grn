#' Create Figure 5 for paper: Illustrate modules / clustering for top network
devtools::load_all(quiet = TRUE)
dir.create("figs", showWarnings = FALSE)

# Load network
net <- readRDS("network/wgcna/kna10-mc2-cor-s100.rds") |> transform()

# Sizes
make_graph(net) |>
  cluster_walktrap() |>
  sizes()

# a
pdf("figs/fig5a.pdf", 5, 5)
cluster_plot(net)
dev.off()

# b
# See scripts/panther-export.r to reproduce json files
da <- parse_panther("data/go/cluster-a.json")
db <- parse_panther("data/go/cluster-b.json")
da$Cluster <- "1"
db$Cluster <- "2"
go <- rbind(da, db)
pdf("figs/fig5b.pdf", 7, 4)
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
  fwrite("figs/fig5c.csv")
options(op)
