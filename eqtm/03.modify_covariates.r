library(data.table)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/eqtm/data/")
files <- list.files(pattern = "covariates.txt$")
cancers <- gsub("_covariates.txt$", "", files)

for(cancer in cancers){
    data <- read.table(paste0(cancer, "_covariates.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    data<-data[!rownames(data) %in% c("gender","age_at_initial_pathologic_diagnosis"),]
    fwrite(data, paste0(cancer, "_covariates_peer.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
}