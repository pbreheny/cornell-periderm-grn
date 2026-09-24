#' Export AUC/pAUC gains for creating Fig 3/4 heatmaps
v <- c("TF", "method", "AUC", "AUC_Gain", "pAUC", "pAUC_Gain", "strength")

all_chip <- pr_master |>
  _[strength %in% seq(0, 4, by = 0.5)] |>
  _[method == "glasso", method := "gLASSO"] |>
  _[method == "genie", method := "GENIE3"] |>
  _[method == "inferelator", method := "Inferelator"] |>
  _[method == "cor", method := "WGCNA"] |>
  _[, first(.SD), .(TF, method, strength)] |>
  _[, ..v]

all_ko <- pr_master_ko |>
  _[strength %in% seq(0, 4, by = 0.5)] |>
  _[method == "glasso", method := "gLASSO"] |>
  _[method == "genie", method := "GENIE3"] |>
  _[method == "inferelator", method := "Inferelator"] |>
  _[method == "cor", method := "WGCNA"] |>
  _[, first(.SD), .(TF, method, strength)] |>
  _[, ..v]

fwrite(all_chip, file = "~/chip.csv")
fwrite(all_ko, file = "~/ko.csv")
