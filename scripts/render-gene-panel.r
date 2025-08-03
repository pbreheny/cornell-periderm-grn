render_gene_panel = function(filter, ...) {
  rmarkdown::render(
    "scripts/gene-panel.r",
    params = list(filter = filter),
    output_file = paste0("gene-panel","-", filter, ".html"),
    knit_root_dir=getwd(), ...)
}

render_gene_panel('all')
render_gene_panel('weak')
render_gene_panel('strong')
