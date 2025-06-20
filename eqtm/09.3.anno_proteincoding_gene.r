library(dplyr)
library(data.table)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/eqtm")
load("./00.id2symbol_anno.RData")
#------------------- add protein keywords
    # load function annotation from uniprot
    kw_lines <- readLines("./uniprot_sprot_final.txt")
    ids_raw <- kw_lines[seq(1, length(kw_lines), by = 2)]
    ids_raw <- sub("^>", "", ids_raw)
    kws <- kw_lines[seq(2, length(kw_lines), by = 2)]
    id_list <- c()
    kw_list <- c()
    for (i in seq_along(ids_raw)) {
        id_split <- unlist(strsplit(ids_raw[i], ";"))
        kw_value <- kws[i]
        id_list <- c(id_list, id_split)
        kw_list <- c(kw_list, rep(kw_value, length(id_split)))
    }
    kw_df <- data.frame(ensembl = id_list, kw = kw_list, stringsAsFactors = FALSE)
    kw_df <- unique(kw_df,by="ensembl")
    # add function to anno
    anno <- merge(anno, kw_df, by = "ensembl", all.x = TRUE)
#------------------- filter potential reg
    load("../apa_reg/apa_regulators.RData")
    anno <- as.data.table(anno)
    kws<- unique(trimws(unique(unlist(strsplit(anno[symbol %in% tranddb, kw], ";")))))#use tranddb is better
    kw_all<-fread("./uniprot_keywords.tsv")
    kws<-kws[kws %in% unique(kw_all[kw_all$Category%in%c("Biological process","Molecular function"),Name])]
    kws<-c("mRNA processing","mRNA splicing","RNA-binding","Transcription","Transcription regulation","DNA-binding","Transcription termination","Chromatin regulator")# i choose these items
    escape_regex <- function(strings) {str_replace_all(strings, "([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1")}
    pattern <- paste0("\\b(", paste0(escape_regex(kws), collapse = "|"), ")\\b")
    anno$potential <- mapply(function(kw, sym) {
        kw_list <- trimws(unlist(strsplit(as.character(kw), ";")))
        any(kw_list %in% kws) || (sym %in% apa_reg)
    }, anno$kw, anno$symbol)
    anno[, matched_kw := sapply(strsplit(kw, ";"), function(x) paste(intersect(trimws(x), kws), collapse = ";"))]
save(anno,file="./00.symbol2funtion_anno.RData")