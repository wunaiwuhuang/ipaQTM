library(data.table)
library(dplyr)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
source("2.1.0.distribute_enrich_fun.r") #load build_enrich_matrix function

# prepare matrix
load("./genomic_feature_anno.Rdata") #load anno
cancers <- gsub("_mdnause.txt$", "", list.files("../data/",pattern = "_mdnause.txt$"))  
types <- c("cis","trans")
sig_samcol <- "snps"
bg <- anno
bg_samcol <- "id"
bg_regcol <- "feature"
region <- c("3'UTR", "5'UTR", "IGR", "TSS1000", "TSS2000", "TSS3000", "1stExon", "OtherExon", "1stIntron", "OtherIntron")
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
save(cis_df, trans_df, file = "genomic_feature_data.Rdata")