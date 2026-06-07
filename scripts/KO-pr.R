#' Exploration of the PR curves for KO and chip-seq data
#' looks first at gold standards generated based on lfsr < 0.01 then gold standards with the 1000 smallest lfsr
#'

irf6_chip100kb_gs <- make_gs('data/chip/irf6_chip_100kb.rds')
grhl3_chip100kb_gs <- make_gs('data/chip/grhl3_chip_100kb.rds')
tfap2a_chip100kb_gs <- make_gs('data/chip/tfap2a_chip_100kb.rds')

# lfsr < 0.01 based KO gold standards
cebpb_KO_gs <- make_gs('data/knockout/cebpb_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
gata3_KO_gs <- make_gs('data/knockout/gata3_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
grhl3_KO_gs <- make_gs('data/knockout/grhl3_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
irf6_KO_gs <- make_gs('data/knockout/irf6_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
klf17_KO_gs <- make_gs('data/knockout/klf17_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
tfap2a_KO_gs <- make_gs('data/knockout/tfap2a_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
znf750_KO_gs <- make_gs('data/knockout/znf750_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
zbtb10_KO_gs <- make_gs('data/knockout/zbtb10_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
tead3b_KO_gs <- make_gs('data/knockout/tead3b_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
klf6a_KO_gs <- make_gs('data/knockout/klf6a_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)
grhl1_KO_gs <- make_gs('data/knockout/grhl1_KO.rds', metric = 'lfsr', selectby = 'maxval', threshold = 0.01)

irf6_combo_100kb_gs <- combine_gs(irf6_chip100kb_gs, irf6_KO_gs)
grhl3_combo_100kb_gs <- combine_gs(grhl3_chip100kb_gs, grhl3_KO_gs)
tfap2a_combo_100kb_gs <- combine_gs(tfap2a_chip100kb_gs, tfap2a_KO_gs)

chip_100_p_fns <- pr_plot(combine_pr(gs = list(irf6_chip100kb_gs,grhl3_chip100kb_gs,tfap2a_chip100kb_gs),
                                     files = c('network/wgcna/mc2/kna10-cor-s100.rds'),
                                     collapse = TRUE), group = 'Network')
KO_p_kna10 <- pr_plot(combine_pr(gs = list(cebpb_KO_gs, gata3_KO_gs, grhl3_KO_gs,
                                           irf6_KO_gs, klf17_KO_gs, tfap2a_KO_gs,
                                           znf750_KO_gs, zbtb10_KO_gs, tead3b_KO_gs,
                                           klf6a_KO_gs, grhl1_KO_gs),
                                 files = c('network/wgcna/mc2/kna10-cor-s100.rds'),
                                 collapse = FALSE), group = 'TF', facet = 'TF')

# smallest 1000 lfsr based KO gold standards
cebpb_KO_gs <- make_gs('data/knockout/cebpb_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
gata3_KO_gs <- make_gs('data/knockout/gata3_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
grhl3_KO_gs <- make_gs('data/knockout/grhl3_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
irf6_KO_gs <- make_gs('data/knockout/irf6_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
klf17_KO_gs <- make_gs('data/knockout/klf17_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
tfap2a_KO_gs <- make_gs('data/knockout/tfap2a_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
znf750_KO_gs <- make_gs('data/knockout/znf750_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
zbtb10_KO_gs <- make_gs('data/knockout/zbtb10_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
tead3b_KO_gs <- make_gs('data/knockout/tead3b_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
klf6a_KO_gs <- make_gs('data/knockout/klf6a_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)
grhl1_KO_gs <- make_gs('data/knockout/grhl1_KO.rds', metric = 'lfsr', selectby = 'N', threshold = 1000)

irf6_combo_100kb_gs <- combine_gs(irf6_chip100kb_gs, irf6_KO_gs)
grhl3_combo_100kb_gs <- combine_gs(grhl3_chip100kb_gs, grhl3_KO_gs)
tfap2a_combo_100kb_gs <- combine_gs(tfap2a_chip100kb_gs, tfap2a_KO_gs)

chip_100_p_fns <- pr_plot(combine_pr(gs = list(irf6_chip100kb_gs,grhl3_chip100kb_gs,tfap2a_chip100kb_gs),
                                     files = c('network/wgcna/mc2/kna10-cor-s100.rds'),
                                     collapse = TRUE), group = 'Network')
KO_p_kna10 <- pr_plot(combine_pr(gs = list(cebpb_KO_gs, gata3_KO_gs, grhl3_KO_gs,
                                           irf6_KO_gs, klf17_KO_gs, tfap2a_KO_gs,
                                           znf750_KO_gs, zbtb10_KO_gs, tead3b_KO_gs,
                                           klf6a_KO_gs, grhl1_KO_gs),
                                 files = c('network/wgcna/mc2/kna10-cor-s100.rds'),
                                 collapse = FALSE), group = 'TF', facet = 'TF')

