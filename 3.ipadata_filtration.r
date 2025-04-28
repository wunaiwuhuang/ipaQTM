# IPA表达文件过滤去除的标准如下：
#Remove ipa events:
#1.With usage missing rate >0.5
#2.On sex chromosomes and chrM chromosomes
#3.With standard deviation <0.01
#Remove samples:
#1.With usage missing rate >0.5

library(data.table)
setwd("/data1/wuguojia/data/IPA_QTM_tcga")
files <- list.files(path="./data/",pattern = "ipause.txt$")

for (file in files) {
    cat("Processing file:", file, "\n")
    data <- fread(paste0("./data/",file))

    # calculate missing rate for rows (ignore the first column)
    missing_rate_row <- apply(data[, -1, with = FALSE], 1, function(x) mean(is.na(x)))
    data <- data[missing_rate_row <= 0.5, ]  # filter rows based on missing rate

    # remove sex chromosomes and chrM chromosomes
    colnames(data)[1] <- "id"  # rename the first column to "id"
    data <- data[!grepl("^chr[XYM]:", data$id), ] 

    # calculate standard deviation for rows (ignore the first column)
    sd_values <- apply(data[, -1, with = FALSE], 1, sd, na.rm = TRUE)
    data <- data[sd_values >= 0.01, ]  # filter rows based on standard deviation

    # calculate missing rate for cols (ignore the first column)
    missing_rate_col <- apply(data[, -1, with = FALSE], 2, function(x) mean(is.na(x)))
    data <- data[, c(TRUE, missing_rate_col <= 0.5), with = FALSE]  # keep the first column and filter others

    # save the filtered data
    fwrite(data, file=paste0("./data/",file), sep = "\t", row.names = FALSE, quote = FALSE)
}