library(data.table)
library(dplyr)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
source("2.1.0.distribute_enrich_fun.r") #load build_enrich_matrix function

# prepare matrix
load("./cgisland_feature_anno_2.Rdata") # we use anno_2, because it is more reliable
load("../00.probe_map_filtered.RData")
anno<-anno[anno$probeID %in% probe_map$id,] # we set all considered probe as background    
cancers <- gsub("_mdnause.txt$", "", list.files("../data/",pattern = "_mdnause.txt$"))  
types <- c("cis","trans")
sig_samcol <- "snps"
bg <- anno
bg_samcol <- "probeID"
bg_regcol <- "CGIposition"
region <- c("island","shelf","shore","opensea")
# build results
cis_df <- list()
trans_df <- list()
for (cancer in cancers) {
    for (type in types) {
        sig_file <- paste0("../", type, "_sig/", cancer, "_", type, "filter_noex.txt")
        if (file.exists(sig_file)) {
        message("Processing ", type, " for ", cancer)
        res <- build_enrich_matrix(
            cancer = cancer,
            type = type,
            sig_file = sig_file,
            sig_samcol = sig_samcol,
            bg = bg,
            bg_samcol = bg_samcol,
            bg_regcol = bg_regcol,
            region = region
        )
        if (type == "cis") {
            cis_df[[cancer]] <- res
        } else {
            trans_df[[cancer]] <- res
        }
        } else {
            message("File not found: ", sig_file)
        }
    }
}
save(cis_df, trans_df, file = "cgisland_feature_data.Rdata")