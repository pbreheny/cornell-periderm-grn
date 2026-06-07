#' Plot heatmap of auc by prior strength and gene
#'
#' @param DT        Output from `auc_table()` (data.table)
#' @param metric    What type of AUC should be plotted? ('AUC', 'AUC_Gain', 'pAUC', or 'pAUC_Gain')
#' @param grid      Should grid lines be plotted on the heatmap? (TRUE or FALSE)
#'
#' @returns   A `ggplot2` object
#'
#' @examples
#' gs <- list(make_gs('data/knockout/grhl3_KO.rds', 'lfsr', selectby = "maxval", threshold = 0.01),
#'            make_gs('data/knockout/irf6_KO.rds', 'lfsr', selectby = "maxval", threshold = 0.01))
#' f <- c(list.files('network/wgcna/', pattern = "kna10-mc2-cor*", full.names = TRUE),
#'        list.files('network/genie/', pattern = "kna10-gne-*", full.names = TRUE))
#' DT <- combine_pr(f, gs, collapse = FALSE)
#' auc_table(DT) |> master_hm('pAUC_Gain')
#' @export

master_hm <- function(DT, metric, grid = FALSE){

  lab <- str_replace(metric, "_", " ")
  heat.pal <- c("#440154FF", "#3B528BFF",
                "#21908CFF", "#5DC863FF", "#FDE725FF")
  probs <- seq(0, 1, length = length(heat.pal))
  quantiles <- quantile(DT[, ..metric], probs, na.rm = T)
  heat.pal.values <- (quantiles - min(quantiles)) / (max(quantiles) - min(quantiles))

  hm <- ggplot(DT) +
    geom_raster(aes(x = as.factor(TF), y = as.factor(strength), fill = get(metric))) +
    scale_fill_gradientn(values = heat.pal.values,
                         colours = heat.pal,
                         name = lab) +
    facet_grid(method ~ ., switch = "y") +
    scale_x_discrete(expand = c(0, 0), name = "Genes") +
    scale_y_discrete(expand = c(0, 0), name = "Prior Strength")

  if(grid){
    hm <- hm +
      geom_vline(xintercept = 1:length(unique(DT[, alpha])) + 0.5) +
      geom_hline(yintercept = 1:length(unique(DT[, strength])) + 0.5)
  }

  hm <- hm +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(face = "bold"),
          panel.border = element_rect(fill=NA, linewidth = 1.5),
          panel.spacing = unit(0, "lines"))

  hm
}
