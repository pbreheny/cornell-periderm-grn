# --- KO gold standards -----
cebpb_KO_gs <- make_gs(
  "data/knockout/cebpb_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
gata3_KO_gs <- make_gs(
  "data/knockout/gata3_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
grhl3_KO_gs <- make_gs(
  "data/knockout/grhl3_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
irf6_KO_gs <- make_gs(
  "data/knockout/irf6_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
klf17_KO_gs <- make_gs(
  "data/knockout/klf17_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
tfap2a_KO_gs <- make_gs(
  "data/knockout/tfap2a_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
znf750_KO_gs <- make_gs(
  "data/knockout/znf750_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
zbtb10_KO_gs <- make_gs(
  "data/knockout/zbtb10_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
tead3b_KO_gs <- make_gs(
  "data/knockout/tead3b_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
klf6a_KO_gs <- make_gs(
  "data/knockout/klf6a_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)
grhl1_KO_gs <- make_gs(
  "data/knockout/grhl1_KO.rds",
  metric = "lfsr",
  selectby = "N",
  threshold = N,
  targets = targets
)

# --- ChIP gold standards -----
irf6_chip100kb_gs <- make_gs("data/chip/irf6_chip_100kb.rds", targets = targets)
grhl3_chip100kb_gs <- make_gs("data/chip/grhl3_chip_100kb.rds", targets = targets)
tfap2a_chip100kb_gs <- make_gs("data/chip/tfap2a_chip_100kb.rds", targets = targets)

# --- ChIP + KO combined gold standards -----
irf6_combo_100kb_gs <- combine_gs(irf6_chip100kb_gs, irf6_KO_gs)
grhl3_combo_100kb_gs <- combine_gs(grhl3_chip100kb_gs, grhl3_KO_gs)
tfap2a_combo_100kb_gs <- combine_gs(tfap2a_chip100kb_gs, tfap2a_KO_gs)
