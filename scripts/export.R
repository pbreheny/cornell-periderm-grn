#' Export networks, as requested by Sunil

readRDS('network/wgcna/kna10-mc2-cor-s100.rds') %>%
  export_network('~/tmp/sunil/wgcna-s100-list.csv', as_list=TRUE)
readRDS('network/wgcna/kna10-mc2-cor-s100.rds') %>%
  export_network('~/tmp/sunil/wgcna-s100-matrix.csv', as_list=FALSE)
readRDS('network/inferelator/kna10-mc2-bbsr-s100.rds') %>%
  export_network('~/tmp/sunil/bbsr-s100-list.csv', as_list=TRUE)
readRDS('network/inferelator/kna10-mc2-bbsr-s100.rds') %>%
  export_network('~/tmp/sunil/bbsr-s100-matrix.csv', as_list=FALSE)
readRDS('network/gl/kna10-mc2-s100-a050.rds') %>%
  export_network('~/tmp/sunil/gl-s100-a050-list.csv', as_list=TRUE)
readRDS('network/gl/kna10-mc2-s100-a050.rds') %>%
  export_network('~/tmp/sunil/gl-s100-a050-matrix.csv', as_list=FALSE)
readRDS('network/genie/kna10-mc2-s100.rds') %>%
  export_network('~/tmp/sunil/genie-s100-list.csv', as_list=TRUE)
readRDS('network/genie/kna10-mc2-s100.rds') %>%
  export_network('~/tmp/sunil/genie-s100-matrix.csv', as_list=FALSE)
