setwd("/data1/wuguojia/data/IPA_QTM_tcga/mediation")
load("./all_sig_result.RData")

brca<- subset(result, cancer_type == "BRCA")
gene <-as.data.frame(table(brca$gene))