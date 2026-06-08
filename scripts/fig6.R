#' Create Figure 7: Is network structure associated with human cleft genes?
#'
#' This is similar in content to scripts/gene-panel.r, but creates figures out
#' of the most important findings
source("scripts/setup.R")

# Read in network summary
net <- readRDS("network/wgcna/kna10-mc2-cor-s100.rds") |> transform()
tf <- summary(net, tf_only = TRUE) # Don't want clustering skewed by targets
non_tf <- summary(net, tf_only = FALSE) |>
  _[!(Gene %in% colnames(net))]

# 7a: All TFs
p <- wilcox.test(Centrality ~ Panel, tf)$p.value
f6a <- qbox(tf, Centrality, Panel) +
  xlab("Centrality") +
  ggtitle(glue("All TFs: {formatP(p, label=TRUE)}")) +
  ylab("") +
  scale_y_discrete(
    labels = c(
      "TRUE" = glue("TFs in OFC\ngene panel ({nrow(tf[Panel == TRUE])})"),
      "FALSE" = glue("TFs not in OFC\ngene panel ({nrow(tf[Panel == FALSE])})")
    )
  )

# 7b: Non-TFs (effectors)
p <- wilcox.test(Edgesum ~ Panel, non_tf)$p.value
f6b <- qbox(non_tf, Edgesum, Panel) +
  xlab("Edge sum") +
  ggtitle(glue("Effector genes: {formatP(p, label=TRUE)}")) +
  ylab("") +
  scale_y_discrete(
    labels = c(
      "TRUE" = glue(
        "Effectors in OFC\ngene panel ({nrow(non_tf[Panel == TRUE])})"
      ),
      "FALSE" = glue(
        "Effectors not in OFC\ngene panel ({nrow(non_tf[Panel == FALSE])})"
      )
    )
  )

# 7c: Cell adhesion TFs only (module 1)
p <- wilcox.test(Centrality ~ Panel, tf[Cluster == 1])$p.value
f6c <- qbox(tf[Cluster == 1], Centrality, Panel) +
  xlab("Centrality") +
  ggtitle(str_glue("Cell adhesion TFs: {formatP(p, label=TRUE)}")) +
  ylab("") +
  scale_y_discrete(
    labels = c(
      "TRUE" = glue(
        "Module 1 TFs\nin OFC\ngene panel ({nrow(tf[Panel == TRUE & Cluster == 1])})"
      ),
      "FALSE" = glue(
        "Module 1 TFs\nnot in OFC\ngene panel ({nrow(tf[Panel == FALSE & Cluster == 1])})"
      )
    )
  )

# Assemble into a figure
gridExtra::grid.arrange(f6a, f6b, f6c) |>
  ggsave(filename = "figs/fig6.pdf", height = 6, width = 4)
