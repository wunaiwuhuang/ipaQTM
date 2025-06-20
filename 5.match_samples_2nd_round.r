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
    #overlap
    samples <- intersect(colnames(mDNA_data)[-1], colnames(IPA_data)[-1])
    mDNA_data <- mDNA_data[, c(colnames(mDNA_data)[1], samples), with = FALSE]
    IPA_data <- IPA_data[, c(colnames(IPA_data)[1], samples), with = FALSE]
    #save
    fwrite(mDNA_data, mDNA_file, sep = "\t", row.names = FALSE, quote = FALSE)
    fwrite(IPA_data, IPA_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("所有文件处理完成。\n")
