#' Plot heatmap of AUC by Strength and Alpha
#'
#' @param DT   Output from `auc_table()` (data.table)
#'
#' @returns   A `ggplot2` object
#'
#' @examples
#' irf6_chip100kb_gs <- make_gs('data/chip/irf6_EVL_chip_100kb.rds')
#' grhl3_chip100kb_gs <- make_gs('data/chip/grhl3_EVL_chip_100kb.rds')
#' tfap2a_chip100kb_gs <- make_gs('data/chip/tfap2a_EVL_chip_100kb.rds')
#' grhl3_KO_gs <- make_gs('data/knockout/grhl3_KO.rds', threshold = 0.01)
#' irf6_KO_gs <- make_gs('data/knockout/irf6_KO.rds', threshold = 0.01)
#' tfap2a_KO_gs <- make_gs('data/knockout/tfap2a2c_KO.rds', threshold = 0.01)
#'
#' irf6_combo_100kb_gs <- combine_gs(irf6_chip100kb_gs, irf6_KO_gs)
#' grhl3_combo_100kb_gs <- combine_gs(grhl3_chip100kb_gs, grhl3_KO_gs)
#' tfap2a_combo_100kb_gs <- combine_gs(tfap2a_chip100kb_gs, tfap2a_KO_gs)
#'
#' DT <- combine_pr(gs = list(irf6_combo_100kb_gs, grhl3_combo_100kb_gs, tfap2a_combo_100kb_gs),
#'   files = c('data/network/gl/kna10-mc2-s000-a001.rds',
#'             'data/network/gl/kna10-mc2-s100-a001.rds',
#'             'data/network/gl/kna10-mc2-s200-a001.rds',
#'             'data/network/gl/kna10-mc2-s300-a001.rds',
#'             'data/network/gl/kna10-mc2-s400-a001.rds',
#'             'data/network/gl/kna10-mc2-s000-a050.rds',
#'             'data/network/gl/kna10-mc2-s100-a050.rds',
#'             'data/network/gl/kna10-mc2-s200-a050.rds',
#'             'data/network/gl/kna10-mc2-s300-a050.rds',
#'             'data/network/gl/kna10-mc2-s400-a050.rds',
#'             'data/network/gl/kna10-mc2-s000-a100.rds',
#'             'data/network/gl/kna10-mc2-s100-a100.rds',
#'             'data/network/gl/kna10-mc2-s200-a100.rds',
#'             'data/network/gl/kna10-mc2-s300-a100.rds',
#'             'data/network/gl/kna10-mc2-s400-a100.rds'),
#'   collapse = FALSE) %>%
#'   auc_table()
#' auc_heat(DT[TF == 'grhl3'])
#' @export

auc_heat <- function(DT) {
  heat.pal <- c("#440154FF", "#3B528BFF",
                "#21908CFF", "#5DC863FF", "#FDE725FF")
  probs <- seq(0, 1, length = length(heat.pal))
  quantiles <- quantile(DT$AUC_Gain, probs, na.rm = T)
  heat.pal.values <- (quantiles - min(quantiles)) /
    (max(quantiles) - min(quantiles))

  if(max(DT[['strength']]) > 100){
    DT[, strength := strength / 100]
  }

  ggplot(DT, aes(x = alpha, y = strength, fill = AUC_Gain)) +
    geom_raster() +
    scale_fill_viridis() +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 12, face = 'plain', hjust = 0.5),
          axis.title = element_text(face = 'plain'),
          legend.position = "right") +
    ylab("Strength") +
    xlab("Alpha")
}
