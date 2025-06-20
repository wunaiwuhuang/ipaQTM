library(dplyr)
library(data.table)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/mediation")
files <- list.files(path="./result/",pattern = ".*_mediate_all.RData",full.names = F)
load("../eqtm/00.symbol2funtion_anno.RData")
putates<-unique(na.omit(anno[anno$potential,symbol]))#potential apa_regulators
result <- list()
for(file in files){
    cancer <- gsub("_mediate_all.RData", "", file)
    load(paste0("./result/",file))
    rtmp <- rtmp[rtmp$gene%in%putates, ] 
    # p.adjust
    rtmp$FDR_ACME<-p.adjust(rtmp$ACME_pval,method="BH")
    rtmp$FDR_TE<-p.adjust(rtmp$total_effect_pvalr,method="BH")
    result[[cancer]] <- rtmp
}
result <- do.call(rbind, result)
# add gene location and cpg location
    # gene location
    load("/data1/wuguojia/data/IPA_QTM_tcga/eqtm/00.id2symbol_anno.RData")
    load("/data1/wuguojia/data/IPA_QTM_tcga/eqtm/00.proteincoder.RData")
    geneanno<-left_join(geneanno,anno,by=c("id_nv" = "ensembl"))%>%na.omit()
    geneanno<-geneanno[,c("symbol","chr","start","end")]%>%.[!duplicated(geneanno$symbol),]
    colnames(geneanno)<-c("gene","chr_gene","start_gene","end_gene")
    result<-left_join(result,geneanno,by=c("gene"="gene"))
    # cpg location
    load("../00.probe_map_filtered.RData")
    cg_loc<-probe_map[,c("id","chr","start","end")]
    colnames(cg_loc)<-c("cpg","chr_cpg","start_cpg","end_cpg")
    cg_loc<-cg_loc[!duplicated(cg_loc$cpg), ]
    result<-left_join(result, cg_loc, by = "cpg")
    #left chr_gene == chr_cpg
    result <- result[result$chr_gene == result$chr_cpg, ]
#-------------- significance filtration
    ## step1: remove NA 
    result <- result[!is.na(result$ACME_Estimate) & result$ACME_Estimate != "", ]
    ## step2: TE FDR<=0.05
    result <- result[result$FDR_TE <= 0.05, ]
    result <- result[sign(result$TE_CI_lower) == sign(result$TE_CI_Upper), ]
    ## step3: ACME FDR<=0.05
    result <- result[result$FDR_ACME <= 0.05, ]
    result <- result[sign(result$ACME_CI_lower) == sign(result$ACME_CI_Upper), ]
    ## step4: remove harmful ACME
    result <- result[sign(result$ACME_Estimate) == sign(result$total_effect_Estimate) & (sign(result$ACME_Estimate) == sign(result$ADE_Estimate) | result$ADE_pval > 0.05),] # ACME and TE must follow same direction,but ADE won't obey the rule if it is not significant.
    ## step5: prop.mediate should be reasonable
    result <- result[result$Prop.Mediated_Estimate > 0 & result$Prop.Mediated_Estimate < 1, ]
    result <- result[result$Prop.Mediated_pval <= 0.05, ]
    result <- result[sign(result$Prop_CI_lower) == sign(result$Prop_CI_Upper), ]
#---------------- threshold filtration
    threshold<-seq(0, 0.8, by = 0.01)
    df<-data.frame(matrix(0,nrow=length(threshold),ncol=5,dimnames=list(NULL,c("thresh","all_gene","all_reg","brca_gene","brca_reg"))))
    df$thresh=threshold
    load("../apa_reg/apa_regulators.RData")#known apa regulator
    for(i in 1:dim(df)[1]){
        thresh<-df$thresh[i]
        # prop >= thresh, and gene should be more than once
        tpm_all<-result[result$Prop.Mediated_Estimate>=thresh,]
        gene_counts <- table(tpm_all$gene)
        genes_to_keep <- names(gene_counts[gene_counts > 1])
        tpm_all <- tpm_all[tpm_all$gene %in% genes_to_keep, ]
        # prop >= thresh
        tpm_brca<-tpm_all[tpm_all$cancer_type=="BRCA",]

        df$all_gene[i]=length(unique(tpm_all$gene))
        df$all_reg[i]=length(intersect(unique(tpm_all$gene),apa_reg))
        df$brca_gene[i]=length(unique(tpm_brca$gene))
        df$brca_reg[i]=length(intersect(unique(tpm_brca$gene),apa_reg))
    }
    df$all_per=df$all_reg/df$all_gene *100
    df$brca_per=df$brca_reg/df$brca_gene *100
    # choose thresh=0.05 as a good filtration
    prop_thresh=0.1
    pair_thresh=1
    result<-result[result$Prop.Mediated_Estimate>=prop_thresh,]
    gene_counts <- table(result$gene)
    genes_to_keep <- names(gene_counts[gene_counts > pair_thresh])
    result <- result[result$gene %in% genes_to_keep, ]
# check overlap with apa_reg
    load("../apa_reg/apa_regulators.RData")
    gene<-unique(result$gene)
    #brca<-result[result$cancer_type=="BRCA",]
    #gene<-unique(brca$gene)
    gene_sets <- list(
        Gene = gene,
        TranDB = tranddb,
        Harmonizome = harmonizome)
    upset(
        fromList(gene_sets),
        order.by = "freq",
        sets = c("Gene", "TranDB", "Harmonizome"),
        mainbar.y.label = "Number of Shared Genes",
        sets.x.label = "Total Genes per Set",
        text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2))
# save data
save(result, file = "./all_sig_result.RData")