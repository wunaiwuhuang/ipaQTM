library(data.table)
library(dplyr)
library(tidyr)

sp_path<-"/data1/wuguojia/data/IPA_QTM_tcga/spqtm/data/"
setwd(sp_path)
mdna_path<-"/data1/wuguojia/data/IPA_QTM_tcga/data/"

file<-list.files(path=mdna_path,pattern = "_mdnause\\.txt$")
cancers <- gsub("_mdnause\\.txt$", "", file)  # remove suffix
for(cancer in cancers){
    cat("正在处理癌症类型: ", cancer, "\n")
    #load files
    mdna_file <- fread(paste0(mdna_path,cancer, "_mdnause.txt"))%>%as.data.frame()
    sp_file <- fread(paste0(sp_path,cancer, "_spuse.txt"))%>%as.data.frame()
    
    # first overlap with mdna
    sample <- intersect(colnames(mdna_file)[-1],colnames(sp_file)[-1])
    mdna_file <- mdna_file[, c(colnames(mdna_file)[1], sample)]
    sp_file <- sp_file[, c(colnames(sp_file)[1], sample)]
    
    # remove row missing rate > 0.1
    missing_rate_row <- apply(sp_file[, -1], 1, function(x) mean(is.na(x)))
    sp_file <- sp_file[missing_rate_row <= 0.1, ]  # filter rows based on missing rate

    # remove standard deviation < 0.01
    sp_file <- sp_file %>% mutate(across(-id, as.numeric)) # convert to numeric
    sd_values <- apply(sp_file[, -1], 1, sd, na.rm = TRUE)
    sp_file <- sp_file[sd_values >= 0.01, ]  # filter rows based on standard deviation

    # i download the Minimum Average Expression Percentage <0.1 data, and have removed On sex ,chrM chromosomes in 01.ger_location.r, so there is no need to do it again

    # remove col missing rate >0.5
    missing_rate_col <- apply(sp_file[, -1], 2, function(x) mean(is.na(x)))
    sp_file <- sp_file[, c(TRUE, missing_rate_col <= 0.5)]  # keep the first column and filter others

    # second overlap with mdna
    sample <- intersect(colnames(mdna_file)[-1],colnames(sp_file)[-1])
    mdna_file <- mdna_file[, c(colnames(mdna_file)[1], sample)]
    sp_file <- sp_file[, c(colnames(sp_file)[1], sample)]

    # save files
    # mdna location
    mdna_loc <- fread(paste0(mdna_path,cancer, "_mdnaloc.txt"))%>%as.data.frame()
    fwrite(mdna_loc, file=paste0(sp_path, cancer, "_mdnaloc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    # mdna usage file
    fwrite(mdna_file, file=paste0(sp_path, cancer, "_mdnause.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    # sp location
    sp_loc <- fread(paste0(sp_path,cancer, "_sploc.txt"))%>%as.data.frame()
    sp_loc <- sp_loc[sp_loc$id %in% sp_file$id,]  
    sp_loc <- sp_loc[match(sp_file$id, sp_loc$id),]# keep in same order as sp_file
    fwrite(sp_loc, file=paste0(sp_path, cancer, "_sploc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    # sp usage file
    fwrite(sp_file, file=paste0(sp_path, cancer, "_spuse.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}