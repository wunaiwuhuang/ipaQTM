# goto conda enviroment Rnew
library(data.table)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/eqtm")
cisfiles <- list.files("./cis_sig/", pattern = "_CisSig_noex\\.txt$")
transfiles <- list.files("./trans_sig/", pattern = "_TransSig_noex\\.txt$")

#prepare to convert ensembl id to gene symbol
    mart <- useMart("ensembl")
    dataset <- as.data.frame(listDatasets(mart))
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    input_list <- as.data.frame(listFilters(mart))
    output_list <- as.data.frame(listAttributes(mart))
    input <- "ensembl_gene_id"
    output <- "external_gene_name"
# list all genes
    gene_id=character(0)
    for(i in 1:length(cisfiles)){
        cat("Processing file: ", cisfiles[i] , transfiles[i],"\n")
        cis <- fread(paste0("./cis_sig/", cisfiles[i])) %>% as.data.frame()
        trans <- fread(paste0("./trans_sig/",transfiles[i])) %>% as.data.frame()
        x <- as.character(unique(c(cis$gene, trans$gene)))
        gene_id <- c(gene_id,x)
    }
    gene_id <- unique(sub("\\..*", "", gene_id))
# convert the ensembl id to gene symbol
    gene_symbol00 <- getBM(attributes = c(input, output), filters = input, values = gene_id, mart = mart)
    gene_symbol01 <- mapIds(org.Hs.eg.db,keys = gene_id,column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")
    # create annotation file
    gene_symbol00_unique <- gene_symbol00[!duplicated(gene_symbol00$ensembl_gene_id), ]
    gene_symbol01_unique <- gene_symbol01[!duplicated(names(gene_symbol01))]
    anno <- data.frame(ensembl = gene_id, symbol = NA)
    anno$symbol <- gene_symbol01_unique[anno$ensembl]
    na_index <- which(is.na(anno$symbol))
    if (length(na_index) > 0) {
        map00 <- setNames(gene_symbol00_unique$external_gene_name, gene_symbol00_unique$ensembl_gene_id)
        anno$symbol[na_index] <- map00[anno$ensembl[na_index]]
    }
    save(anno, file = "./id2symbol_anno.RData")
# modify qtm files
    load("./id2symbol_anno.RData")
    anno$symbol[trimws(anno$symbol) == ""] <- NA
    colnames(anno) <- c("id", "symbol")
    for(i in 1:length(cisfiles)){
        cat("Processing file: ", cisfiles[i] , transfiles[i],"\n")
        cis <- fread(paste0("./cis_sig/", cisfiles[i])) %>% as.data.frame()
        trans <- fread(paste0("./trans_sig/",transfiles[i])) %>% as.data.frame()
        cis$id <- gsub("\\..*", "", cis$gene)
        trans$id <- gsub("\\..*", "", trans$gene)

        cis <- left_join(cis,anno,by = "id")
        trans <- left_join(trans,anno,by = "id")

        fwrite(cis, file = paste0("./cis_sig/", cisfiles[i]), sep = "\t", row.names = F, col.names = T, quote = F)
        fwrite(trans, file = paste0("./trans_sig/", transfiles[i]), sep = "\t", row.names = F, col.names = T, quote = F)
    }