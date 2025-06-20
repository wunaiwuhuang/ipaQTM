# ----------------------- get protein coding genes
    library(data.table)
    library(dplyr)
    setwd("/data1/wuguojia/data/IPA_QTM_tcga/eqtm/data")
    geneanno<- fread("../../00.gencode.v48.annotation.foruse.txt",header=F,sep="\t")%>%as.data.frame()
    geneanno<-geneanno[,c(1,6,2,4,5,8,10,12)]
    colnames(geneanno)<-c("chr","chain","source","start","end","id","type","name")
    geneanno[]<-lapply(geneanno,function(col){
        col <- as.character(col)
        col <-gsub(";","",col)
        col <-gsub('"',"",col)
        return(col)
    })
    geneanno <- subset(geneanno,geneanno$type=="protein_coding"&geneanno$chr %in% paste0("chr", 1:22))
    #remove duplicate id
    geneanno$id_nv<-gsub("\\..*","",geneanno$id)
    geneanno <- geneanno[!duplicated(geneanno$id_nv), ] # in fact no duplicate id, annotation is good
    geneanno$name[grep("^ENSG",geneanno$name)]<-NA
    save(geneanno,file="../00.proteincoder.RData")
# ------------------------------------------------------------------------
library(data.table)
library(dplyr)

mdna_path<-"/data1/wuguojia/data/IPA_QTM_tcga/data/"
setwd("/data1/wuguojia/data/IPA_QTM_tcga/eqtm/data")
files<-list.files(path=".",pattern = "_tpmuse\\.txt$")
cancers <- gsub("_tpmuse\\.txt$", "", files)  # remove suffix

load("../00.proteincoder.RData")
for(cancer in cancers){
    cat("processing cancer: ", cancer, "\n")
#mdnause files and mdnaloc files are ripe, thereis no use to filter
#------------------- tpmuse & mdnause
    tpmuse <- fread(paste0(cancer,"_tpmuse.txt"),header=T,sep="\t",check.names=F) %>% as.data.frame()
    mdnause <- fread(paste0(mdna_path,cancer,"_mdnause.txt"),header=T,sep="\t",check.names=F) %>% as.data.frame()
    # first overlap
    sample <- intersect(colnames(mdnause)[-1],colnames(tpmuse)[-1])
    mdnause <- mdnause[, c(colnames(mdnause)[1], sample)]
    tpmuse <- tpmuse[, c(colnames(tpmuse)[1], sample)]
    # id in geneanno files
    tpmuse$id <-gsub("\\..*","",tpmuse$id)
    tpmuse <- tpmuse[tpmuse$id %in% geneanno$id_nv,]
    tpmuse <- tpmuse[!duplicated(tpmuse$id), ]
    # remove row missing rate > 0.1
    missing_rate_row <- apply(tpmuse[, -1], 1, function(x) mean(is.na(x)))
    tpmuse <- tpmuse[missing_rate_row <= 0.1, ]
    # remove standard deviation < 0.1
    tpmuse <- tpmuse %>% mutate(across(-id, as.numeric)) # convert to numeric
    sd_values <- apply(tpmuse[, -1], 1, sd, na.rm = TRUE)
    tpmuse <- tpmuse[sd_values >= 0.1, ]
    # remove over 80% of samples are less than 0.1
    tpmuse[, -1] <- lapply(tpmuse[, -1], as.numeric)
    expr_count <- apply(tpmuse[, -1], 1, function(x) sum(x >= 0.1, na.rm = TRUE))
    threshold <- 0.2 * (ncol(tpmuse) - 1)
    tpmuse <- tpmuse[expr_count >= threshold, ]
    # remove col missing rate >0.5
    missing_rate_col <- apply(tpmuse[, -1], 2, function(x) mean(is.na(x)))
    tpmuse <- tpmuse[, c(TRUE, missing_rate_col <= 0.5)]
    # second overlap
    sample <- intersect(colnames(mdnause)[-1],colnames(tpmuse)[-1])
    mdnause <- mdnause[, c(colnames(mdnause)[1], sample)]
    tpmuse <- tpmuse[, c(colnames(tpmuse)[1], sample)]
#------------------- tpmloc & mdnaloc
    mdnaloc <- fread(paste0(mdna_path,cancer,"_mdnaloc.txt"),header=T,sep="\t",check.names=F) %>% as.data.frame()
    tpmloc <- geneanno[geneanno$id_nv %in% tpmuse$id,]
    tpmloc <- tpmloc[,c("id_nv","chr","start","end")]
    tpmloc$start <- as.numeric(tpmloc$start)
    tpmloc$end <- as.numeric(tpmloc$end)
    # in same order as tpmuse
    tpmloc <- tpmloc[match(tpmuse$id, tpmloc$id_nv),]
#------------------- save data
    fwrite(tpmuse, file=paste0(cancer, "_tpmuse.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    fwrite(mdnause, file=paste0(cancer, "_mdnause.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    fwrite(tpmloc, file=paste0(cancer, "_tpmloc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    fwrite(mdnaloc, file=paste0(cancer, "_mdnaloc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}