library(data.table)
library(dplyr)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/eqtm")
cisfiles <- list.files("./cis_sig", pattern = ".txt", full.names = F)
transfiles <- list.files("./trans_sig", pattern = ".txt", full.names = F)

# modify qtm files
load("./00.id2symbol_anno.RData")
colnames(anno) <- c("id", "symbol")
for(i in 1:length(cisfiles)){
    cat("Processing file: ", cisfiles[i] , transfiles[i],"\n")
    cis <- fread(paste0("./cis_sig/", cisfiles[i])) %>% as.data.frame()
    trans <- fread(paste0("./trans_sig/",transfiles[i])) %>% as.data.frame()
    cis$id <- gsub("\\..*", "", cis$gene)
    trans$id <- gsub("\\..*", "", trans$gene)

    cis <- left_join(cis,anno,by = "id")
    trans <- left_join(trans,anno,by = "id")
    # remove symbol is na
    cis <- na.omit(cis)
    trans <- na.omit(trans)
    
    fwrite(cis, file = paste0("./cis_sig/", cisfiles[i]), sep = "\t", row.names = F, col.names = T, quote = F)
    fwrite(trans, file = paste0("./trans_sig/", transfiles[i]), sep = "\t", row.names = F, col.names = T, quote = F)
}