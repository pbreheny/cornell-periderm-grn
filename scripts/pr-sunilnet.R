devtools::load_all()

# Setup -------------------------------------------------------------------
## Retrieve Gold Standard Networks ##
# Gold Standard for irf6
irf6_gs <- readRDS('rds/irf6_gs.rds')

# Gold Standard for tfap2a
tfap2a_gs <- readRDS('rds/tfap2a_gs.rds')

# Gold Standard for tfap2a2c
tfap2a2c_gs <- readRDS('rds/tfap2a2c_gs.rds')

# Gold Standard for tfap2aMO
tfap2aMO_gs <- readRDS('rds/tfap2aMO_gs.rds')


# Gold Standard for klf17
klf17_gs <- readRDS('rds/klf17_gs.rds')

# Gold Standard for cebpb
cebpb_gs <- readRDS('rds/cebpb_gs.rds')

# Gold Standard for gata3
gata3_gs <- readRDS('rds/gata3_gs.rds')

# Gold Standard for grhl3
grhl3_gs <- readRDS('rds/grhl3_gs.rds')

# Gold Standard for znf750
znf750_gs <- readRDS('rds/znf750_gs.rds')

# Target Genes
# targets <- fread("data/zfTargetGenes_S2982.txt", header = FALSE)
targets <- fread("data/targetGenes4283.txt", header = FALSE)

## prepare Sunil's LASSO-StARS network for PR and ROC curve construction
# moderate prior
network_mod <- data.table::fread("./data/GPAE_FIMO_ModPrior_TFmRNA.tsv")
network_irf6_mod <- prep_network_StARS(GS_ash = irf6_gs, network_input = network_mod, 
                                       tf = "irf6", target_genes = targets$V1)
network_tfap2a_mod <- prep_network_StARS(GS_ash = tfap2a_gs, network_input = network_mod, 
                                         tf = "tfap2a", target_genes = targets$V1)
network_tfap2a2c_mod <- prep_network_StARS(GS_ash = tfap2a2c_gs, network_input = network_mod, 
                                           tf = c("tfap2a", "tfap2c"), target_genes = targets$V1)
network_tfap2aMO_mod <- prep_network_StARS(GS_ash = tfap2aMO_gs, network_input = network_mod,
                                           tf = "tfap2a", target_genes = targets$V1)
network_klf17_mod <- prep_network_StARS(GS_ash = klf17_gs, network_input = network_mod, 
                                           tf = "klf17", target_genes = targets$V1)
network_cebpb_mod <- prep_network_StARS(GS_ash = cebpb_gs, network_input = network_mod, 
                                        tf = "cebpb", target_genes = targets$V1)
network_gata3_mod <- prep_network_StARS(GS_ash = gata3_gs, network_input = network_mod, 
                                        tf = "gata3", target_genes = targets$V1)
network_grhl3_mod <- prep_network_StARS(GS_ash = grhl3_gs, network_input = network_mod, 
                                        tf = "grhl3", target_genes = targets$V1)
network_znf750_mod <- prep_network_StARS(GS_ash = znf750_gs, network_input = network_mod, 
                                        tf = "znf750", target_genes = targets$V1)

# high prior
network_high <- data.table::fread("./data/GPAE_FIMO_HighPrior_TFmRNA.tsv")
network_irf6_high <- prep_network_StARS(GS_ash = irf6_gs, network_input = network_high, tf = "irf6",
                                        target_genes = targets$V1)
network_tfap2a_high <- prep_network_StARS(GS_ash = tfap2a_gs, network_input = network_high, tf = "tfap2a",
                                          target_genes = targets$V1)
network_tfap2a2c_high <- prep_network_StARS(GS_ash = tfap2a2c_gs, network_input = network_high, 
                                           tf = c("tfap2a", "tfap2c"), target_genes = targets$V1)
network_tfap2aMO_high <- prep_network_StARS(GS_ash = tfap2aMO_gs, network_input = network_high,
                                           tf = "tfap2a", target_genes = targets$V1)
network_klf17_high <- prep_network_StARS(GS_ash = klf17_gs, network_input = network_high, 
                                        tf = "klf17", target_genes = targets$V1)
network_cebpb_high <- prep_network_StARS(GS_ash = cebpb_gs, network_input = network_high, 
                                        tf = "cebpb", target_genes = targets$V1)
network_gata3_high <- prep_network_StARS(GS_ash = gata3_gs, network_input = network_high, 
                                        tf = "gata3", target_genes = targets$V1)
network_grhl3_high <- prep_network_StARS(GS_ash = grhl3_gs, network_input = network_high, 
                                        tf = "grhl3", target_genes = targets$V1)
network_znf750_high <- prep_network_StARS(GS_ash = znf750_gs, network_input = network_high, 
                                         tf = "znf750", target_genes = targets$V1)

# no prior
network_no <- data.table::fread("./data/GPAE_FIMO_NoPrior_TFmRNA.tsv")
network_irf6_no <- prep_network_StARS(GS_ash = irf6_gs, network_input = network_no, tf = "irf6",
                                      target_genes = targets$V1)
network_tfap2a_no <- prep_network_StARS(GS_ash = tfap2a_gs, network_input = network_no, tf = "tfap2a",
                                        target_genes = targets$V1)
network_tfap2a2c_no <- prep_network_StARS(GS_ash = tfap2a2c_gs, network_input = network_no, 
                                           tf = c("tfap2a", "tfap2c"), target_genes = targets$V1)
network_tfap2aMO_no <- prep_network_StARS(GS_ash = tfap2aMO_gs, network_input = network_no,
                                           tf = "tfap2a", target_genes = targets$V1)
network_klf17_no <- prep_network_StARS(GS_ash = klf17_gs, network_input = network_no, 
                                        tf = "klf17", target_genes = targets$V1)
network_cebpb_no <- prep_network_StARS(GS_ash = cebpb_gs, network_input = network_no, 
                                        tf = "cebpb", target_genes = targets$V1)
network_gata3_no <- prep_network_StARS(GS_ash = gata3_gs, network_input = network_no, 
                                        tf = "gata3", target_genes = targets$V1)
network_grhl3_no <- prep_network_StARS(GS_ash = grhl3_gs, network_input = network_no, 
                                        tf = "grhl3", target_genes = targets$V1)
network_znf750_no <- prep_network_StARS(GS_ash = znf750_gs, network_input = network_no, 
                                         tf = "znf750", target_genes = targets$V1)

# export the processed Sunil network
write.table(network_irf6_mod, "./results/sunilNet.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)

# export genes that are predicted positive in Sunil's network
PredPos_SunilNet <- network_irf6[which(network_irf6$Pred_Pos == 1),
                                 'SYMBOL']
write.table(PredPos_SunilNet, "./results/PosGenes_SunilNet.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)

# moderate prior
ROC_mod_dat_irf6 <- build_ROC(network_irf6_mod)
ROC_mod_dat_tfap2a <- build_ROC(network_tfap2a_mod)
ROC_mod_dat_tfap2a2c <- build_ROC(network_tfap2a2c_mod)
ROC_mod_tafp2aMO <- build_ROC(network_tfap2aMO_mod)
ROC_mod_dat_klf17 <- build_ROC(network_klf17_mod)
ROC_mod_dat_cebpb <- build_ROC(network_cebpb_mod)
ROC_mod_dat_gata3 <- build_ROC(network_gata3_mod)
ROC_mod_dat_grhl3 <- build_ROC(network_grhl3_mod)
ROC_mod_dat_znf750 <- build_ROC(network_znf750_mod)

# high prior
ROC_high_dat_irf6 <- build_ROC(network_irf6_high)
ROC_high_dat_tfap2a <- build_ROC(network_tfap2a_high)
ROC_high_dat_tfap2a2c <- build_ROC(network_tfap2a2c_high)
ROC_high_tafp2aMO <- build_ROC(network_tfap2aMO_high)
ROC_high_dat_klf17 <- build_ROC(network_klf17_high)
ROC_high_dat_cebpb <- build_ROC(network_cebpb_high)
ROC_high_dat_gata3 <- build_ROC(network_gata3_high)
ROC_high_dat_grhl3 <- build_ROC(network_grhl3_high)
ROC_high_dat_znf750 <- build_ROC(network_znf750_high)


# no prior
### No gene symbols from this network that are in the target genes list ###
#ROC_no_dat <- build_ROC(network_irf6_no)
#ROC_no_dat_tfap2a <- build_ROC(network_tfap2a_no)
#ROC_mod_dat_tfap2a2c <- build_ROC(network_tfap2a2c_no)
#ROC_no_tafp2aMO <- build_ROC(network_tfap2aMO_no)
#ROC_mod_dat_klf17 <- build_ROC(network_klf17_no)
#ROC_mod_dat_cebpb <- build_ROC(network_cebpb_no)
#ROC_mod_dat_gata3 <- build_ROC(network_gata3_no)
#ROC_mod_dat_grhl3 <- build_ROC(network_grhl3_no)
#ROC_mod_dat_znf750 <- build_ROC(network_znf750_no)
