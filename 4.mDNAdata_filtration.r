# 探针下载地址 https://gdc-hub.s3.us-east-1.amazonaws.com/download/HM450.hg38.manifest.gencode.v36.probeMap
# SNP覆盖参考文献：DOI: 10.1038/nature15393 ；数据下载地址：https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/ ; 实际使用从/data1/wangwenhui/pu-DNAm/1000_Genomes_phase3路径获得数据，
# step1 探针文件过滤去除的标准如下：
# Remove cg sites:
# 1.With beta value missing rate >0.05
# 2.mapping to multiple genomic locations  
# 3.On sex ,chrM chromosomes
# 4.containing known SNPs (from the 1000 Genomes phase 3 with MAF > 0.01)
# 5.with a standard deviation <0.05 
# Remove samples:
# 1.With beta value missing rate >0.05

library(data.table)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/")

#generate good probes
probe_map <- fread("./00.HM450.hg38.manifest.gencode.v36.probeMap")
probe_map <- na.omit(probe_map)
probe_map <- probe_map[, -ncol(probe_map), with = FALSE]
colnames(probe_map) <- c("id","gene","chr","start","end")
# remove mapping to multiple genomic locations
probe_map <- unique(probe_map, by = c("id", "start")) #data.table unique function
# remove sex,chrM chromosomes
probe_map <- probe_map[!(chr %in% c("chrX", "chrY", "chrM"))]
# remove known SNPs
snp_data<-list()
for (chr_num in 1:22) {
    snp_file1 <- sprintf("./1000_Genomes_phase3/SNP.info.chr%d.txt", chr_num)
    snp_file2 <- sprintf("./1000_Genomes_phase3/1000GP.hg38.chr%d.frq", chr_num)
    snp_data1 <- fread(snp_file1)
    snp_data2 <- fread(snp_file2)

    setnames(snp_data1, c("chr", "pos", "rsid", "ref", "alt"))
    setnames(snp_data2, c("chr", "rsid", "A1", "A2", "maf", "other"))
    snp_data1 <- snp_data1[grepl("^[A-Za-z]", snp_data1$rsid), ]
    snp_data2 <- snp_data2[grepl("^[A-Za-z]", snp_data2$rsid), ]
    temp <- cbind(snp_data1[,c("chr", "pos", "rsid", "ref", "alt")], snp_data2[,c("A1", "A2", "maf")]) #all rows are corresponding

    temp$maf <- ifelse(temp$ref == temp$A2, temp$maf, 1 - temp$maf)
    temp <- temp[maf > 0.01]
    snp_data[[chr_num]] <- unique(temp,by = "pos")
}
do.call(rbind, snp_data) -> snp_data
snp_data <- snp_data[!is.na(pos)]
probe_map <- probe_map[!snp_data, on = .(chr, start < pos, end >= pos)]
save(probe_map, file = "./00.probe_map_filtered.RData")

#start filtration
load("./00.probe_map_filtered.RData")
files <- list.files(path="./data/",pattern = "mdnause.txt$")
for (file in files) {
    cat("Processing file:", file, "\n")
    data <- fread(paste0("./data/",file))

    # calculate missing rate for rows (ignore the first column)
    missing_rate_row <- apply(data[, -1, with = FALSE], 1, function(x) mean(is.na(x)))
    data <- data[missing_rate_row <= 0.05, ]  # filter rows based on missing rate

    #included in the probe_map
    colnames(data)[1] <- "id"  # rename the first column to "id"
    data <- data[id %in% probe_map$id]

    #calcaulate standard deviation for rows (ignore the first column)
    sd_values <- apply(data[, -1, with = FALSE], 1, sd, na.rm = TRUE)
    data <- data[sd_values >= 0.05, ]  # filter rows based on standard deviation

    # calculate missing rate for cols (ignore the first column)
    missing_rate_col <- apply(data[, -1, with = FALSE], 2, function(x) mean(is.na(x)))
    data <- data[, c(TRUE, missing_rate_col <= 0.05), with = FALSE]  # keep the first column and filter others
    
    # save the filtered data
    fwrite(data, file=paste0("./data/",file), sep = "\t", row.names = FALSE, quote = FALSE)
}