#' Create figure 3
source("scripts/setup.R")

# --- Import enrichment and gold standards -----
enr <- import_enr()
targets <- "./data/raw/TargetGenes_EVL_TRN_Full-Weak-StrongFiltered.xlsx" |>
  readxl::read_excel(range = "A1:A3788") |>
  as.vector() |>
  unlist() |>
  unique()
N <- round(0.25 * length(targets))
source("scripts/load-gs.R")

# --- Precision-recall calculations: Gold standard -----

# Combine (collapsing over all TFs)
pr_combo <- combine_pr(
  gs = list(irf6_combo_100kb_gs, grhl3_combo_100kb_gs, tfap2a_combo_100kb_gs),
  files = c(
    list.files(path = "network/gl", pattern = "^kna10-mc2.*a001\\.rds$", full.names = TRUE),
    list.files(path = "network/wgcna", pattern = "kna10-mc2-cor-*", full.names = TRUE),
    list.files(path = "network/inferelator", pattern = "kna10-mc2-*", full.names = TRUE),
    list.files(path = "network/genie", pattern = "kna10-mc2-*", full.names = TRUE)
  ),
  collapse = TRUE
)

# Separate for each TF
pr_combo_sep <- combine_pr(
  gs = list(irf6_combo_100kb_gs, grhl3_combo_100kb_gs, tfap2a_combo_100kb_gs),
  files = c(
    list.files(path = "network/gl", pattern = "^kna10-mc2.*a001\\.rds$", full.names = TRUE),
    list.files(path = "network/wgcna", pattern = "kna10-mc2-cor-*", full.names = TRUE),
    list.files(path = "network/inferelator", pattern = "kna10-mc2-*", full.names = TRUE),
    list.files(path = "network/genie", pattern = "kna10-mc2-*", full.names = TRUE)
  ),
  collapse = FALSE
)
pr_master <- rbind(pr_combo, pr_combo_sep)
pr_master[,
  TF := factor(
    TF,
    levels = c("irf6,grhl3,tfap2a", "grhl3", "irf6", "tfap2a"),
    labels = c("combined", "grhl3", "irf6", "tfap2a")
  )
]

# --- Precision-recall calculations: KO -----
pr_ko <- combine_pr(
  gs = list(
    cebpb_KO_gs,
    gata3_KO_gs,
    grhl3_KO_gs,
    irf6_KO_gs,
    klf17_KO_gs,
    tfap2a_KO_gs,
    znf750_KO_gs,
    zbtb10_KO_gs,
    tead3b_KO_gs,
    klf6a_KO_gs,
    grhl1_KO_gs
  ),
  files = c(
    list.files(path = "network/gl", pattern = "^kna10-mc2.*a001\\.rds$", full.names = TRUE),
    list.files(path = "network/wgcna", pattern = "kna10-mc2-cor-*", full.names = TRUE),
    list.files(path = "network/inferelator", pattern = "kna10-mc2-*", full.names = TRUE),
    list.files(path = "network/genie", pattern = "kna10-mc2-*", full.names = TRUE)
  ),
  collapse = TRUE
)
pr_ko_sep <- combine_pr(
  gs = list(
    cebpb_KO_gs,
    gata3_KO_gs,
    grhl3_KO_gs,
    irf6_KO_gs,
    klf17_KO_gs,
    tfap2a_KO_gs,
    znf750_KO_gs,
    zbtb10_KO_gs,
    tead3b_KO_gs,
    klf6a_KO_gs,
    grhl1_KO_gs
  ),
  files = c(
    list.files(path = "network/gl", pattern = "^kna10-mc2.*a001\\.rds$", full.names = TRUE),
    list.files(path = "network/wgcna", pattern = "kna10-mc2-cor-*", full.names = TRUE),
    list.files(path = "network/inferelator", pattern = "kna10-mc2-*", full.names = TRUE),
    list.files(path = "network/genie", pattern = "kna10-mc2-*", full.names = TRUE)
  ),
  collapse = FALSE
)
pr_master_ko <- rbind(pr_ko, pr_ko_sep)
tf_levels <- c(
  "cebpb,gata3,grhl3,irf6,klf17,tfap2a,znf750,zbtb10,tead3b,klf6a,grhl1",
  "cebpb",
  "gata3",
  "grhl3",
  "irf6",
  "klf17",
  "tfap2a",
  "znf750",
  "zbtb10",
  "tead3b",
  "klf6a",
  "grhl1"
)
tf_labels <- c(
  "combined",
  "cebpb",
  "gata3",
  "grhl3",
  "irf6",
  "klf17",
  "tfap2a",
  "znf750",
  "zbtb10",
  "tead3b",
  "klf6a",
  "grhl1"
)
pr_master_ko[, TF := factor(TF, levels = tf_levels, labels = tf_labels)]

# --- Fig 3a -----

f3a <- pr_master[strength == 1] |>
  _[Network == "network/gl/kna10-mc2-s100-a001.rds", Network := "gLASSO"] |>
  _[Network == "network/genie/kna10-mc2-gne-s100.rds", Network := "GENIE3"] |>
  _[Network == "network/inferelator/kna10-mc2-bbsr-s100.rds", Network := "Inferelator 3.0"] |>
  _[Network == "network/wgcna/kna10-mc2-cor-s100.rds", Network := "WGCNA"] |>
  pr_plot(group = "Network", facet = "TF", col = 4) +
  guides(color = guide_legend(title = ""))

# --- Fig 3b -----

f3b <- pr_master[TF != "Combined"] |>
  _[strength %in% 0:4] |>
  _[method == "glasso", method := "gLASSO"] |>
  _[method == "genie", method := "GENIE3"] |>
  _[method == "inferelator", method := "Inferelator"] |>
  _[method == "cor", method := "WGCNA"] |>
  master_hm(metric = "pAUC_Gain", grid = FALSE)

# --- Fig 3c -----

f3c <- pr_master_ko[!(TF %in% c("combined", "irf6", "grhl3", "tfap2a"))] |>
  _[strength %in% 0:4] |>
  _[method == "glasso", method := "gLASSO"] |>
  _[method == "genie", method := "GENIE3"] |>
  _[method == "inferelator", method := "Inferelator"] |>
  _[method == "cor", method := "WGCNA"] |>
  master_hm(metric = "pAUC_Gain", grid = FALSE)

ggsave("figs/fig3.pdf", f3a / f3b / f3c, width = 8, height = 10)
