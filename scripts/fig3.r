# Figure 3
devtools::load_all(quiet = TRUE)
library(superheat)
library(ggplotify)
library(ggpubr)
library(MetBrewer)

enr <- import_enr()
targets <- readxl::read_excel('./data/raw/TargetGenes_EVL_TRN_Full-Weak-StrongFiltered.xlsx', range = "A1:A3788") %>%
  as.vector() %>%
  unlist() %>%
  unique()
N <- round(0.25 * length(targets))


# KO gold standards
cebpb_KO_gs <- make_gs('data/knockout/cebpb_KO.rds', metric = 'lfsr',
                       selectby = 'N', threshold = N, targets = targets)
gata3_KO_gs <- make_gs('data/knockout/gata3_KO.rds', metric = 'lfsr',
                       selectby = 'N', threshold = N, targets = targets)
grhl3_KO_gs <- make_gs('data/knockout/grhl3_KO.rds', metric = 'lfsr',
                       selectby = 'N', threshold = N, targets = targets)
irf6_KO_gs <- make_gs('data/knockout/irf6_KO.rds', metric = 'lfsr',
                      selectby = 'N', threshold = N, targets = targets)
klf17_KO_gs <- make_gs('data/knockout/klf17_KO.rds', metric = 'lfsr',
                       selectby = 'N', threshold = N, targets = targets)
tfap2a_KO_gs <- make_gs('data/knockout/tfap2a_KO.rds', metric = 'lfsr',
                        selectby = 'N', threshold = N, targets = targets)
znf750_KO_gs <- make_gs('data/knockout/znf750_KO.rds', metric = 'lfsr',
                        selectby = 'N', threshold = N, targets = targets)
zbtb10_KO_gs <- make_gs('data/knockout/zbtb10_KO.rds', metric = 'lfsr',
                        selectby = 'N', threshold = N, targets = targets)
tead3b_KO_gs <- make_gs('data/knockout/tead3b_KO.rds', metric = 'lfsr',
                        selectby = 'N', threshold = N, targets = targets)
klf6a_KO_gs <- make_gs('data/knockout/klf6a_KO.rds', metric = 'lfsr',
                       selectby = 'N', threshold = N, targets = targets)
grhl1_KO_gs <- make_gs('data/knockout/grhl1_KO.rds', metric = 'lfsr',
                       selectby = 'N', threshold = N, targets = targets)
# chip gold standards
irf6_chip100kb_gs <- make_gs('data/chip/irf6_chip_100kb.rds', targets = targets)
grhl3_chip100kb_gs <- make_gs('data/chip/grhl3_chip_100kb.rds', targets = targets)
tfap2a_chip100kb_gs <- make_gs('data/chip/tfap2a_chip_100kb.rds', targets = targets)
# chip + ko combined gold standards
irf6_combo_100kb_gs <- combine_gs(irf6_chip100kb_gs, irf6_KO_gs)
grhl3_combo_100kb_gs <- combine_gs(grhl3_chip100kb_gs, grhl3_KO_gs)
tfap2a_combo_100kb_gs <- combine_gs(tfap2a_chip100kb_gs, tfap2a_KO_gs)

## PR data tables ##
### combined chip+ko data ###
# entry for every network in gl-a001, wgcna, inferelator, and genie, collapsed over all TFs
pr_combo <- combine_pr(
  gs = list(
    irf6_combo_100kb_gs,
    grhl3_combo_100kb_gs,
    tfap2a_combo_100kb_gs
  ),
  files = c(
    list.files(
      path = 'data/network/gl',
      pattern = '^kna10-mc2.*a001\\.rds$',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/wgcna/',
      pattern = 'kna10-mc2-cor-*',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/inferelator',
      pattern = 'kna10-mc2-*',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/genie/',
      pattern = 'kna10-mc2-*',
      full.names = TRUE)
    ),
  collapse = TRUE
)
# entry for every network in gl-a001, wgcna, inferelator, and genie, separate for each TF
pr_combo_sep <- combine_pr(
  gs = list(
    irf6_combo_100kb_gs,
    grhl3_combo_100kb_gs,
    tfap2a_combo_100kb_gs
  ),
  files = c(
    list.files(
      path = 'data/network/gl',
      pattern = '^kna10-mc2.*a001\\.rds$',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/wgcna/',
      pattern = 'kna10-mc2-cor-*',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/inferelator',
      pattern = 'kna10-mc2-*',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/genie/',
      pattern = 'kna10-mc2-*',
      full.names = TRUE)
    ),
  collapse = FALSE
)
pr_master <- rbind(pr_combo, pr_combo_sep)
pr_master[, TF := factor(TF, levels = c('irf6,grhl3,tfap2a', 'grhl3', 'irf6', 'tfap2a'))]

### ko only data ###
pr_combo_ko <- combine_pr(
  gs = list(
    cebpb_KO_gs, gata3_KO_gs, grhl3_KO_gs, irf6_KO_gs, klf17_KO_gs, tfap2a_KO_gs,
    znf750_KO_gs, zbtb10_KO_gs, tead3b_KO_gs, klf6a_KO_gs, grhl1_KO_gs
  ),
  files = c(
    list.files(
      path = 'data/network/gl',
      pattern = '^kna10-mc2.*a001\\.rds$',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/wgcna/',
      pattern = 'kna10-mc2-cor-*',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/inferelator',
      pattern = 'kna10-mc2-*',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/genie/',
      pattern = 'kna10-mc2-*',
      full.names = TRUE
    )
  ),
  collapse = TRUE
)
pr_combo_sep_ko <- combine_pr(
  gs = list(
    cebpb_KO_gs, gata3_KO_gs, grhl3_KO_gs, irf6_KO_gs, klf17_KO_gs, tfap2a_KO_gs,
    znf750_KO_gs, zbtb10_KO_gs, tead3b_KO_gs, klf6a_KO_gs, grhl1_KO_gs
  ),
  files = c(
    list.files(
      path = 'data/network/gl',
      pattern = '^kna10-mc2.*a001\\.rds$',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/wgcna/',
      pattern = 'kna10-mc2-cor-*',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/inferelator',
      pattern = 'kna10-mc2-*',
      full.names = TRUE
    ),
    list.files(
      path = 'data/network/genie/',
      pattern = 'kna10-mc2-*',
      full.names = TRUE
    )
  ),
  collapse = FALSE
)
pr_master_ko <- rbind(pr_combo_ko, pr_combo_sep_ko)
tf_levels <- c(
  'cebpb,gata3,grhl3,irf6,klf17,tfap2a,znf750,zbtb10,tead3b,klf6a,grhl1',
  'cebpb', 'gata3', 'grhl3', 'irf6', 'klf17', 'tfap2a', 'znf750', 'zbtb10', 'tead3b', 'klf6a', 'grhl1'
)
tf_labels <- c(
  'combined', 'cebpb', 'gata3', 'grhl3', 'irf6', 'klf17', 'tfap2a', 'znf750', 'zbtb10', 'tead3b', 'klf6a', 'grhl1'
)
pr_master_ko[, TF := factor(TF, levels = tf_levels, labels = tf_labels)]

# --- plotting --------------------------------------------

pr_plot(pr_master[strength == 1], group = 'Network', facet = 'TF', col = 4)

# AUC line plots
tb_auc <- auc_table(pr_master)
tb_auc_ko <- auc_table(pr_master_ko)

auc_plot(tb_auc, 'method', 'TF', 'pAUC_Gain', ncol = 4)
auc_plot(tb_auc_ko, 'method', 'TF', 'pAUC_Gain', ncol = 4)


# --- heat maps -------------------------------------------

# combined chip+ko data

# all prior strengths
master_hm(pr_master, metric = "pAUC_Gain", grid = TRUE)

# only prior strengths 1, 2, 3, 4
master_hm(pr_master[strength == 1], metric = "pAUC_Gain", grid = TRUE)

# KO only
master_hm(pr_master_ko[strength == 0.5], metric = "pAUC_Gain", grid = TRUE)

# --- mini networks/subgraphs -----------------------------

genes <- c(
  'cebpb', 'gata3', 'grhl1', 'grhl3', 'irf6', 'klf6a', 'klf17',
  'tead3b', 'tfap2a', 'zbtb10', 'znf750'
)

# network plot
readRDS('data/network/wgcna/kna10-mc2-cor-s100.rds') |>
  transform() |>
  make_graph() |>
  plot_modules(
    genes = genes,
    percentile = 0.1,
    layout = "circle",
    n_lab = 11,
    background = "white"
  )

readRDS('data/network/wgcna/kna10-mc2-cor-s100.rds') |>
  transform() |>
  make_graph() |>
  plot_modules(percentile = 0.1, n_lab = 11, background = "white")
