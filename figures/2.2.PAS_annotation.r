library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
mDNA_files <- list.files("../data/",pattern = "_mdnause.txt$")
cancer <- gsub("_mdnause.txt$", "", mDNA_files)  

ipaloc<-list()
for(c in cancer){ipaloc[[c]] <- read.table(paste0("../data/",c,"_ipaloc.txt"),header = T,sep = "\t",quote = "")}
ipaloc <- do.call(rbind,ipaloc)
ipaloc <- unique(ipaloc,by = c("id","chr"))
ipaloc$chain <- str_extract(ipaloc$id, ":[+-]:")
ipaloc$chain <- gsub(":", "", ipaloc$chain)
#if chain is +,then use end,else use start
ipaloc$pos <- ifelse(ipaloc$chain == "+",ipaloc$end,ipaloc$start)
length = 50
#focus on upstream, + chain use minus,- chain use plus
ipaloc$start <- ifelse(ipaloc$chain == "+",ipaloc$pos-length,ipaloc$pos-1) # bed左开右闭
ipaloc$end <- ifelse(ipaloc$chain == "+",ipaloc$pos,ipaloc$pos+length-1)
ipaloc$score = 1
ipaloc<-ipaloc[,c("chr","start","end","id","score","chain")]
#save to bed file
write.table(ipaloc,paste0("../00.ipaloc_upstream50.bed"),sep = "\t",quote = F,row.names = F,col.names = F)

# $ cd /data1/wuguojia/data/IPA_QTM_tcga
# $ bedtools getfasta -fi 00.hg38.fa -bed 00.ipaloc_upstream50.bed -fo 00.ipaloc_sequence.bed -s -name -tab

seq <- read.table("../00.ipaloc_sequence.bed",header = F,sep = "\t",quote = "")
ipaloc$seq <- toupper(seq$V2)
pasmotif <- c("AATAAA","ATTAAA","TATAAA","AATATA","AGTAAA","AAGAAA","AATAAT","AATACA","CATAAA","AATGAA","AACAAA","GATAAA","ACTAAA","ATTATA","AATAGA","AATAAG","ATTACA","AACAAG","AAAAAA","TTTAAA")
# find if seq contain pasmotif
ipaloc$pas <- sapply(ipaloc$seq, function(x) {
  hits <- pasmotif[sapply(pasmotif, function(motif) grepl(motif, x, fixed = TRUE))]
  if (length(hits) > 0) {
    return(paste(hits, collapse = ";"))
  } else {
    return(NA)  # 没有匹配的可以返回 NA 或 ""
  }
})
# add ipaloc$pos
ipaloc$pos <- ifelse(ipaloc$chain == "+",ipaloc$end,ipaloc$start+1)
save(ipaloc,file = "./pas_motif.RData")

load("./pas_motif.RData")
load("../00.probe_map_filtered.RData")
ipaloc <- na.omit(ipaloc)
paslist <- list()
for(c in cancer){
  print(paste0("Processing ",c))
  cis <- read.table(paste0("../cis_sig/", c, "_cisfilter_noex.txt"),header = T,sep = "\t",quote = "")
  cis <- cis[which(cis$gene %in% ipaloc$id),]
  # add ipaloc$seq ,start and end info
  cis <- merge(cis,ipaloc[,c("id","pos","pas")],by.x = "gene",by.y = "id",all.x = T)
  cis <- cis[,c("gene","pos","pas","snps")]
  threshold <- 100 #100bp
  cis$up <- cis$pos - threshold
  cis$down <- cis$pos + threshold
  # add probe_map start and end info
  cis <- merge(cis,probe_map[,c("id","start","end")],by.x = "snps",by.y = "id",all.x = T)
  # if up <= start and down >= end, then keep this row
  cis <- cis[which(cis$up <= cis$start & cis$down >= cis$end),]
  # if no entry, then fill in NA
  if(nrow(cis) == 0){
    cis <- data.frame(gene = c, pos = NA, pas = NA, snps = NA, up = NA, down = NA, start = NA, end = NA)
  }
  cis$cancer <- c
  paslist[[c]] <- cis
}
paslist <- do.call(rbind,paslist)
paslist <- paslist[,c("cancer","pas")]
# remove ; in pas, expand list
paslist <- tidyr::separate_rows(paslist, pas, sep = ";")
paslist <- na.omit(paslist)
source("./0.color.r")
pas_colors <- generate_cancer_colors(unique(paslist$pas))
pas_count <- paslist %>%group_by(cancer, pas) %>% summarise(count = n(), .groups = "drop")
#reorder
cancer_order <- pas_count %>% group_by(cancer) %>% summarise(total = sum(count)) %>% arrange(desc(total)) %>% pull(cancer)
pas_count$cancer <- factor(pas_count$cancer, levels = cancer_order)
# draw
p <- ggplot(pas_count, aes(x = cancer, y = count, fill = pas)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = pas_colors) +
      labs(x = "Cancer Type", y = "PAS Motif Count", fill = "PAS Motif") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
ggsave(filename="./2.2.pas_annotation.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)