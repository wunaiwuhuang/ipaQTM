# goto conda enviroment Rnew
library(data.table)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/eqtm")
#prepare to convert ensembl id to gene symbol
    mart <- useMart("ensembl")
    dataset <- as.data.frame(listDatasets(mart))
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    input_list <- as.data.frame(listFilters(mart))
    output_list <- as.data.frame(listAttributes(mart))
    input <- "ensembl_gene_id"
    output <- "external_gene_name"
# list all genes
    load("00.proteincoder.RData")
    gene_id <- unique(geneanno$id_nv)
# convert the ensembl id to gene symbol
    gene_symbol00 <- getBM(attributes = c(input, output), filters = input, values = gene_id, mart = mart)
    gene_symbol01 <- mapIds(org.Hs.eg.db,keys = gene_id,column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")
    gene_symbol02 <- geneanno[,c("id_nv","name")]
    # create annotation file
    anno <- data.frame(ensembl = gene_id, symbol = NA, stringsAsFactors = FALSE)
    # Step 1: 尝试从 gene_symbol02 映射
    map2 <- setNames(trimws(gene_symbol02$name), gene_symbol02$id_nv)
    anno$symbol <- map2[anno$ensembl]
    # Step 2: 对仍为 NA 的尝试从 gene_symbol01
    na_idx <- which(is.na(anno$symbol) | anno$symbol == "")
    map1 <- setNames(trimws(gene_symbol01), names(gene_symbol01))
    anno$symbol[na_idx] <- map1[anno$ensembl[na_idx]]
    # Step 3: 对仍为 NA 的尝试从 gene_symbol00
    na_idx <- which(is.na(anno$symbol) | anno$symbol == "")
    map0 <- setNames(trimws(gene_symbol00$external_gene_name), gene_symbol00$ensembl_gene_id)
    anno$symbol[na_idx] <- map0[anno$ensembl[na_idx]]
    # Step 4: 再次清洗空字符串为 NA
    anno$symbol <- ifelse(is.na(anno$symbol) | trimws(anno$symbol) == "", NA, anno$symbol)

    save(anno, file = "./00.id2symbol_anno.RData")
