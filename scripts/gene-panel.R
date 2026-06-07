#' ---
#' title: Is network structure associated with human cleft genes?
#' author: ''
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#' params:
#'   filter: all
#' ---
#'
#+ include=FALSE, purl=FALSE

# Display options
source("scripts/setup.R")
library(knitr)
library(gtsummary)
knitr::opts_knit$set(aliases = c(h = "fig.height", w = "fig.width"))
knitr::opts_chunk$set(comment = "#", echo = FALSE, message = FALSE,
                      collapse = TRUE, cache = FALSE, tidy = FALSE, fig.align = "center")

# Main network
net <- readRDS('network/wgcna/kna10-cor-s100.rds')
g_full <- make_graph(net)
m <- cluster_walktrap(g_full) |> membership

# Subset network
evl_expr <- import_enr()
keep <- switch (params$filter,
                all = rownames(net),
                weak = evl_expr[!(p_val < 0.05 & avg_log2FC < 0)]$Gene,
                strong = evl_expr[p_val < 0.05 & avg_log2FC > 0]$Gene
)
net <- net[intersect(rownames(net), keep), intersect(colnames(net), keep)]

# Gene panel
gp <- import_gp()

#' # Summary statistics on gene panel
#'
#' Here is an overall breakdown of the genes in the gene panel by level of evidence and whether
#' they are expressed in EVL or not (and if so, whether we consider them a TF or not):

gp[, .(Evidence=factor(evidence, levels=c('Green', 'Amber', 'Red', 'NA')),
       evl=fcase(gene %in% colnames(net), 'TF',
                 gene %in% rownames(net), 'Non-TF',
                 default = 'Not expressed') |>
         factor(levels=c('Not expressed', 'Non-TF', 'TF')))] |>
  tbl_cross(percent='row') |>
  modify_spanning_header(all_stat_cols() ~ "**EVL expression**")

#' # Analysis of transcription factors

# Setup
g <- make_graph(net)
x <- apply(soft(net, 0.1), 2, sum)
tf <- hub_genes(g) |>
  merge(data.table(Gene=names(m), Cluster=m)) |>
  merge(data.table(Gene=names(x), EdgeSum=x))

#' ## Test cluster

# g   List of genes (character vector)
# cl  Cluster assignment
test_cluster <- function(g, dt, label) {
  dt[, .(Panel = factor(Gene %in% g, labels=c('No', 'Yes')), Cluster=factor(Cluster))] |>
    tbl_summary(by=Panel, percent='row') |>
    add_p() |>
    modify_header(label = '') |>
    modify_spanning_header(all_stat_cols() ~ paste0('**', label, '**'))
}

test_cluster(gp$gene, tf, 'In gene panel')
test_cluster(gp[evidence=='Green']$gene, tf, 'Evidence = Green')
test_cluster(gp[evidence=='Green' | evidence == 'Amber']$gene, tf, 'Evidence = Green or Amber')
test_cluster(gp[genepanel==1]$gene, tf, '"genepanel"')
test_cluster(gp[nhs==1]$gene, tf, '"nhs"')
test_cluster(gp[omim==1]$gene, tf, '"omim"')
test_cluster(gp[prevgen==1]$gene, tf, '"prevgen"')

#' ## Test centrality

# g   List of genes (character vector)
# cl  Cluster assignment
test_central <- function(g, dt, label) {
  dt[, .(Panel = factor(Gene %in% g, labels=c('No', 'Yes')), Centrality, EdgeSum)] |>
    tbl_summary(by=Panel) |>
    add_p() |>
    modify_header(label = '') |>
    modify_spanning_header(all_stat_cols() ~ paste0('**', label, '**'))
}

test_central(gp$gene, tf, 'In gene panel')
test_central(gp[evidence=='Green']$gene, tf, 'Evidence = Green')
test_central(gp[evidence %in% c('Green', 'Amber')]$gene, tf, 'Evidence = Green or Amber')

#' Above, we tested for an association with centrality across all genes; below, we repeat the same tests,
#' but only among the genes in cluster 1 (the cell adhesion cluster)

test_central(gp$gene, tf[Cluster==1], 'In gene panel')
test_central(gp[evidence=='Green']$gene, tf[Cluster==1], 'Evidence = Green')
test_central(gp[evidence %in% c('Green', 'Amber')]$gene, tf[Cluster==1], 'Evidence = Green or Amber')

#' # All targets
#'
#' These results are a bit exploratory, as they involve determining clusters for non-tf genes. How best to
#' do this is still a bit of an unanswered question.
#'
#' ## Test cluster

# Setup
target <- data.table(
  Gene    = rownames(net),
  Signal1 = apply(soft(net[,colnames(net) %in% names(m[m==1])], 0.1), 1, mean),
  Signal2 = apply(soft(net[,colnames(net) %in% names(m[m==2])], 0.1), 1, mean),
  Signal3 = apply(soft(net[,colnames(net) %in% names(m[m==3])], 0.1), 1, mean),
  EdgeSum = apply(soft(net, 0.1), 1, sum)
)
target[, Cluster := apply(.SD[, .(Signal1, Signal2, Signal3)], 1, which.max), Gene]

test_cluster(gp$gene, target, 'In gene panel')
test_cluster(gp[evidence=='Green']$gene, target, 'Evidence = Green')
test_cluster(gp[evidence=='Green' | evidence == 'Amber']$gene, target, 'Evidence = Green or Amber')

#' ## Test centrality

test_edgesum <- function(g, dt, label) {
  dt[, .(Panel = factor(Gene %in% g, labels=c('No', 'Yes')), EdgeSum)] |>
    tbl_summary(by=Panel) |>
    add_p() |>
    modify_header(label = '') |>
    modify_spanning_header(all_stat_cols() ~ paste0('**', label, '**')) |>
    modify_footnote(update = everything() ~ NA)
}

test_edgesum(gp$gene, target, 'In gene panel')
test_edgesum(gp[evidence=='Green']$gene, target, 'Evidence = Green')
test_edgesum(gp[evidence=='Green' | evidence == 'Amber']$gene, target, 'Evidence = Green or Amber')

#' Above, we tested for an association with centrality across all target genes;
#' below, we repeat the same tests, but only among the genes in cluster 1
#' (the cell adhesion cluster)

test_edgesum(gp$gene, target[Cluster==1], 'In gene panel')
test_edgesum(gp[evidence=='Green']$gene, target[Cluster==1], 'Evidence = Green')
test_edgesum(gp[evidence %in% c('Green', 'Amber')]$gene, target[Cluster==1], 'Evidence = Green or Amber')
