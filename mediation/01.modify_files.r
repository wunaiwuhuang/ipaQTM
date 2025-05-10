# cp ../../data/*mdnause.txt ./
# cp ../../data/*ipause.txt ./
# cp ../../eqtm/data/*tpmuse.txt ./
# cp ../../trans_sig/*_transfilter_noex.txt ./
# cp ../../eqtm/cis_sig/*_cisfilter_noex.txt ./

library(data.table)
library(dplyr)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/mediation/data")
cancers <- list.files(pattern = ".*_ipause.txt",full.names = F)
cancers <- gsub("_ipause.txt", "", cancers)

#match gene id for tpmuse files
load("/data1/wuguojia/data/IPA_QTM_tcga/eqtm/id2symbol_anno.RData")
anno$symbol[trimws(anno$symbol) == ""] <- NA

for(cancer in cancers){
    cat("dealing with", cancer, "\n")
    # deal with ipaqtm and eqtm files
    ipaqtm <- fread(paste0(cancer, "_transfilter_noex.txt"))
    eqtm <- fread(paste0(cancer, "_cisfilter_noex.txt"))
    eqtm <- eqtm[,c("snps","symbol","statistic","pvalue","FDR","beta","r")]
    colnames(eqtm)[2] <- "gene"
    eqtm$gene[trimws(eqtm$gene) == ""] <- NA
    eqtm <- na.omit(eqtm)
    # get shared cg
    cg <- intersect(unique(ipaqtm$snps), unique(eqtm$snps))  
    ipaqtm <- ipaqtm[which(ipaqtm$snps %in% cg),]
    eqtm <- eqtm[which(eqtm$snps %in% cg),]
    # save files
    fwrite(ipaqtm, paste0(cancer, "_trans_ipaqtm.txt"), sep = "\t", row.names = F, quote = F)
    file.remove(paste0(cancer, "_transfilter_noex.txt"))
    fwrite(eqtm, paste0(cancer, "_cis_eqtm.txt"), sep = "\t", row.names = F, quote = F)
    file.remove(paste0(cancer, "_cisfilter_noex.txt"))
    
    # get selected cg,ipa and gene
    cg <- intersect(unique(ipaqtm$snps), unique(eqtm$snps))
    ipa <- unique(ipaqtm$gene)
    reg <- unique(eqtm$gene)

    # deal with ipause, tpmuse and mdnause files
    ipause <- fread(paste0(cancer, "_ipause.txt"))
    tpmuse <- fread(paste0(cancer, "_tpmuse.txt"))
    mdnause <- fread(paste0(cancer, "_mdnause.txt"))
    colnames(ipause)[1] <- "gene"
    colnames(mdnause)[1] <- "snps"
    colnames(tpmuse)[1] <- "gene"
    # match gene id
    tpm1 <- tpmuse[,1]
    tpm2 <- tpmuse[,2:ncol(tpmuse)]
    tpm1$ensembl <- gsub("\\..*", "", tpm1$gene)
    tpm1 <- left_join(tpm1, anno[,c("ensembl", "symbol")], by = "ensembl")
    tpmuse <- cbind(tpm1, tpm2)
    tpmuse <- tpmuse[!is.na(tpmuse$symbol),-c(1,2)]
    colnames(tpmuse)[1] <- "gene"
    # remove duplicated genes
    ipause <- ipause[!duplicated(ipause$gene),]
    mdnause <- mdnause[!duplicated(mdnause$snps),]
    tpmuse <- tpmuse[!duplicated(tpmuse$gene),]
    # select cg,ipa and gene
    ipause <- ipause[which(ipause$gene %in% ipa),]
    tpmuse <- tpmuse[which(tpmuse$gene %in% reg),]
    mdnause <- mdnause[which(mdnause$snps %in% cg),]
    # get samples
    sample <- intersect(intersect(colnames(ipause[,-1]), colnames(tpmuse[,-1])),colnames(mdnause[,-1]))
    ipause <- ipause[,c(colnames(ipause)[1], sample), with = F] %>% as.data.frame()
    tpmuse <- tpmuse[,c(colnames(tpmuse)[1], sample), with = F] %>% as.data.frame()
    mdnause <- mdnause[,c(colnames(mdnause)[1], sample), with = F] %>% as.data.frame()
    rownames(ipause) <- as.character(ipause[[1]])
    ipause <- ipause[,-1] # remove first column
    rownames(tpmuse) <- as.character(tpmuse[[1]])
    tpmuse <- tpmuse[,-1] # remove first column
    rownames(mdnause) <- as.character(mdnause[[1]])
    mdnause <- mdnause[,-1] # remove first column
    # transpose
    ipause <- as.matrix(ipause) %>% t() %>% as.data.frame()
    ipause <- cbind(sample = rownames(ipause), ipause)
    tpmuse <- as.matrix(tpmuse) %>% t() %>% as.data.frame()
    tpmuse <- cbind(sample = rownames(tpmuse), tpmuse)
    mdnause <- as.matrix(mdnause) %>% t() %>% as.data.frame()
    mdnause <- cbind(sample = rownames(mdnause), mdnause)
    # save
    fwrite(ipause, paste0(cancer, "_ipause.txt"), sep = "\t", row.names = F, quote = F)
    fwrite(tpmuse, paste0(cancer, "_tpmuse.txt"), sep = "\t", row.names = F, quote = F)
    fwrite(mdnause, paste0(cancer, "_mdnause.txt"), sep = "\t", row.names = F, quote = F)
}

cat("all files have been processed.\n")