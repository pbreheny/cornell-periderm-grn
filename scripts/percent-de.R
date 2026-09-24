# Explore percent differentially expressed
ff <- list.files("data/knockout", full.names = TRUE)
for (f in ff) {
  tf <- basename(f) |> str_split_i("_", 1)
  m <- mean(readRDS(f)$lfsr < 0.01)
  print(str_glue("{tf}: {m |> round(2)}"))
}
