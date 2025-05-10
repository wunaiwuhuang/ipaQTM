# 00.human.eCLIP.RBPs.POSTAR3.txt data were downloaded from POSTAR3 database
# 00.human.RBPs.ENCODE.tsv were downloaded from ENCODE database
library(data.table)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
med_path <- "/data1/wuguojia/data/IPA_QTM_tcga/mediation/"
load(paste0(med_path,"./all_sig_result.RData"))

rbp1<-fread("../00.human.eCLIP.RBPs.txt",header = F,sep = "\t")
colnames(rbp1)<-c("chr","start","end","peak_id","chain","gene","source","sample","ref","score")
rbp2<-fread("../00.human.RBPs.ENCODE.tsv",header = T,sep = "\t")

rbpname<-unique(rbp$gene)
med <- unique(result$gene)