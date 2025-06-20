library(data.table)

setwd("/data1/wuguojia/data/IPA_QTM_tcga")
#get cancer types
files <- list.files(path="./data/",pattern = "_ipause\\.txt$")
cancers <- gsub("_ipause\\.txt$", "", files)  # remove suffix

for (cancer in cancers) {
    cat("正在处理癌症类型: ", cancer, "\n")
    #load files
    mDNA_file <- paste0("./data/",cancer, "_mdnause.txt")
    IPA_file <- paste0("./data/",cancer, "_ipause.txt")
    if (!file.exists(mDNA_file) | !file.exists(IPA_file)) {
        cat("缺少必要文件，跳过: ", cancer, "\n")
        next
    }
    mDNA_data <- fread(mDNA_file)
    IPA_data <- fread(IPA_file)    
    #deal with ipafile,remove non-cancer samples and repeated samples
    colnames_to_keep <- colnames(IPA_data)[-1]
    colnames_to_keep <- colnames_to_keep[sapply(colnames_to_keep, function(colname) {
        num <- as.numeric(substr(colname, 14, 15))#eg."TCGA.BH.A0DZ.11A.22R.A089.07"
        return(num < 10)
    })]
    colnames_to_keep <- c(colnames(IPA_data)[1],colnames_to_keep)
    IPA_data <- IPA_data[, colnames_to_keep, with = FALSE]
    colnames(IPA_data)[-1] <- substr(colnames(IPA_data)[-1], 1,12)#eg."TCGA.BH.A0DZ"
    cols_to_keep <- colnames(IPA_data)[!duplicated(colnames(IPA_data))]
    IPA_data <- IPA_data[, ..cols_to_keep] #remove duplicated columns
    colnames(IPA_data)<-gsub("\\.","-",colnames(IPA_data))#change "." to "-"
    #overlap
    samples <- intersect(colnames(mDNA_data)[-1], colnames(IPA_data)[-1])
    mDNA_data <- mDNA_data[, c(colnames(mDNA_data)[1], samples), with = FALSE]
    IPA_data <- IPA_data[, c(colnames(IPA_data)[1], samples), with = FALSE]
    #save
    fwrite(mDNA_data, mDNA_file, sep = "\t", row.names = FALSE, quote = FALSE)
    fwrite(IPA_data, IPA_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("所有文件处理完成。\n")
