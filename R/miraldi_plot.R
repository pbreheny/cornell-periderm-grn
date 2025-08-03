#' Create ChIP/KO/ChIP+KO plot in style of Miraldi2019, Fig2A
#'
#' @param chip_dat    either a ggplot object of chip-seq data as produced by `pr_plot()` or a list of chip-seq gold standard objects as produced by `make_gs()`
#' @param KO_dat      either a ggplot object of KO data as produced by `pr_plot()` or a list of KO gold standard objects as produced by `make_gs()`
#' @param combo_dat   either a ggplot object of combind KO + ChIP-Seq data as produced by `pr_plot()` or a list of ko + ChIP-Seq gold standards as produced by `combine_gs()`
#' @param networks    A vector of `.rds` file names pointing to fitted networks
#' @param title       Title for entire figure (outermost title)
#'
#' @seealso `pr_plot()`
#'
#' @returns a `ggplot2` object
#'
#' @examples
#' irf6_chip100kb_gs <- make_gs('data/chip/irf6_chip_100kb.rds')
#' grhl3_chip100kb_gs <- make_gs('data/chip/grhl3_chip_100kb.rds')
#' tfap2a_chip100kb_gs <- make_gs('data/chip/tfap2a_chip_100kb.rds')
#' cebpb_KO_gs <- make_gs('data/knockout/cebpb_KO.rds', metric = 'lfsr < 0.01')
#' gata3_KO_gs <- make_gs('data/knockout/gata3_KO.rds', metric = 'lfsr < 0.01')
#' grhl3_KO_gs <- make_gs('data/knockout/grhl3_KO.rds', metric = 'lfsr < 0.01')
#' irf6_KO_gs <- make_gs('data/knockout/irf6_KO.rds', metric = 'lfsr < 0.01')
#' klf17_KO_gs <- make_gs('data/knockout/klf17_KO.rds', metric = 'lfsr < 0.01')
#' tfap2a_KO_gs <- make_gs('data/knockout/tfap2a_KO.rds', metric = 'lfsr < 0.01')
#' znf750_KO_gs <- make_gs('data/knockout/znf750_KO.rds', metric = 'lfsr < 0.01')
#' irf6_combo_100kb_gs <- combine_gs(irf6_chip100kb_gs, irf6_KO_gs)
#' grhl3_combo_100kb_gs <- combine_gs(grhl3_chip100kb_gs, grhl3_KO_gs)
#' tfap2a_combo_100kb_gs <- combine_gs(tfap2a_chip100kb_gs, tfap2a_KO_gs)
#' chip_100_p_fns <- pr_plot(combine_pr(gs = list(irf6_chip100kb_gs,grhl3_chip100kb_gs,tfap2a_chip100kb_gs),
#'                           files = c('data/network/genie/fns-genie.rds', 'data/network/wgcna/kna10-cor-s100.rds', 'data/network/mir/fns-mall-s0.rds'),
#'                           collapse = TRUE), group = 'Network')
#' KO_p_fns <- pr_plot(combine_pr(gs = list(cebpb_KO_gs, gata3_KO_gs, grhl3_KO_gs, irf6_KO_gs, klf17_KO_gs, tfap2a_KO_gs, znf750_KO_gs),
#'                     files = c('data/network/genie/fns-genie.rds', 'data/network/wgcna/kna10-cor-s100.rds', 'data/network/mir/fns-mall-s0.rds'),
#'                     collapse = TRUE), group = 'Network')
#' combo_100_p_fns <- pr_plot(combine_pr(gs = list(irf6_combo_100kb_gs, grhl3_combo_100kb_gs, tfap2a_combo_100kb_gs),
#'                            files = c('data/network/genie/fns-genie.rds', 'data/network/wgcna/kna10-cor-s100.rds', 'data/network/mir/fns-mall-s0.rds'),
#'                            collapse = TRUE), group = 'Network')
#'
#' miraldi_plot(chip_dat = list(irf6_chip100kb_gs, grhl3_chip100kb_gs, tfap2a_chip100kb_gs),
#'              KO_dat = list(cebpb_KO_gs, gata3_KO_gs, grhl3_KO_gs, irf6_KO_gs, klf17_KO_gs, tfap2a_KO_gs, znf750_KO_gs),
#'              combo_dat = list(irf6_combo_100kb_gs, grhl3_combo_100kb_gs, tfap2a_combo_100kb_gs),
#'              networks = c('data/network/genie/fns-genie.rds', 'data/network/wgcna/kna10-cor-s100.rds', 'data/network/mir/fns-mall-s0.rds'),
#'              title = 'Farrell Shield Stage 100 kb Flanking Region')
#' miraldi_plot(chip_dat = chip_100_p_fns, KO_dat = KO_p_fns, combo_dat = combo_100_p_fns, title = 'Farrell Shield Stage 100kb Flanking Region')
#' @export
#'

miraldi_plot <- function(chip_dat, KO_dat, combo_dat, networks, title = "Comparison of GS types (Miraldi Panel A)"){
  if (inherits(chip_dat, 'gg')) {
    chip <- chip_dat
  } else {
    chip <- pr_plot(combine_pr(gs = chip_dat, files = networks, collapse = TRUE), group = 'Network')
  }
  if (inherits(KO_dat, 'gg')) {
    KO <- KO_dat
  } else {
    KO <- pr_plot(combine_pr(gs = KO_dat, files = networks, collapse = TRUE), group = 'Network')
  }
  if (inherits(combo_dat, 'gg')) {
    combo <- combo_dat
  } else {
    combo <- pr_plot(combine_pr(gs = combo_dat, files = networks, collapse = TRUE), group = 'Network')
  }
  p <- ggpubr::ggarrange(chip + ggtitle("ChIP-Seq"), KO + ggtitle("Knockout"), combo + ggtitle("ChIP-Seq + KO"),
            nrow = 1, common.legend = TRUE, legend = 'top')
  ggpubr::annotate_figure(p, top = ggpubr::text_grob(title, color = 'black', face = 'bold', size = 14))
}
