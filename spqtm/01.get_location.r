library(data.table)
library(dplyr)
library(tidyr)

sp_path<-"/data1/wuguojia/data/IPA_QTM_tcga/spqtm/data/"
setwd(sp_path)

file<-list.files(path=sp_path,pattern = "^PSI")
cancers <- sub("PSI_download_(.*)\\.txt", "\\1", file)

loc <- fread("../src/00.TCGA_SpliceSeq_Gene_Structure.txt") %>% distinct(Symbol, Exon, .keep_all = TRUE)
colnames(loc)<-c("gene","chr","chain","exon","start","end")
for(cancer in cancers){
    cat("dealing with file: ", cancer, "\n")
    sp_file <- fread(paste0(sp_path,"PSI_download_",cancer, ".txt"))%>%as.data.frame()
    sp_file <- sp_file[,-((ncol(sp_file)-1):ncol(sp_file))] # remove the last two column
    sp_file[sp_file == "null"] <- NA # set null to na
    # remove from_exon or to_exon is na
    sp_file <- sp_file[!is.na(sp_file$from_exon) & !is.na(sp_file$to_exon), ]
    sp_a <- sp_file[,1:10]
    sp_b <- sp_file[,11:ncol(sp_file)]
    #if type is ES/AA/RI/ME/AT, than use to_exon, else use from_exon
    sp_a <- sp_a %>% mutate(exon = ifelse(splice_type %in% c("ES", "AA", "RI", "ME", "AT"), to_exon, from_exon))
    sp_a <- sp_a[,c("symbol","splice_type","exon","from_exon","to_exon")]
    colnames(sp_a)<-c("gene","type","exon","from_exon","to_exon")
    sp_a$exon <- as.numeric(sp_a$exon)
    #left join  sp_a and loc
    sp_a <- sp_a %>% left_join(loc, by = c("gene"="gene", "exon"="exon"))
    #if type is ES/AA/RI/ME/AT, than use start, else use end
    sp_a <- sp_a %>% mutate(pos = ifelse(type %in% c("ES", "AA", "RI", "ME", "AT"), start, end))
    sp_a$chr <- paste0("chr", sp_a$chr)
    sp_a <- sp_a %>% unite("id", c(chr,gene,type,from_exon,to_exon,chain), sep = ":", remove = FALSE) %>% select(id, everything()) # create id column
    sp_a <- sp_a[,c("id","chr","pos")]
    sp_file <- cbind(sp_a, sp_b)
    
    # remove duplicated id
    sp_file <- sp_file[!duplicated(sp_file$id), ]

    #remove chrX chrY and chrM
    sp_file <- sp_file[!grepl("chrX|chrY|chrM", sp_file$chr), ]
    
    #save location file
    loc_file <- sp_file[,c("id","chr","pos")]
    loc_file$start <- loc_file$pos
    loc_file$end <- loc_file$pos+1
    loc_file <- loc_file[,c("id","chr","start","end")]
    fwrite(loc_file, paste0(sp_path, cancer, "_sploc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    #save sp_file
    sp_file <- sp_file[,-(2:3)]
    colnames(sp_file)[-1] <- gsub("\\_","-",colnames(sp_file)[-1]) # change "-" to "-"
    fwrite(sp_file, paste0(sp_path, cancer, "_spuse.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}