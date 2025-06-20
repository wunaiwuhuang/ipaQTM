library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(scales)
library(dplyr)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
load("../mediation/importance_coef_sep.RData")
load("../mediation/all_sig_result.RData")
load("../apa_reg/apa_regulators.RData")
#----------------------- select coefficient
    # test_path select top coefficients
    test_path<-test_path[test_path$overlap_50==max(test_path$overlap_50),]
    test_path<-test_path[test_path$overlap_100==max(test_path$overlap_100),]
    test_path<-test_path[test_path$overlap_200==max(test_path$overlap_200),] 
    test_path<-test_path[test_path$w3==max(test_path$w3),] # assign prop more weight

    # test_str select top coefficients
    test_str<-test_str[test_str$overlap_50==max(test_str$overlap_50),]
    test_str<-test_str[test_str$overlap_100==max(test_str$overlap_100),]
    test_str<-test_str[test_str$w4 == max(test_str$w4),]# assign cancer type more weight
#-------------------------- filter gene
    result <- result %>%
        group_by(gene) %>%
        mutate(cancer_count = n_distinct(cancer_type)) %>%
        ungroup() %>%
        filter(cancer_count > 1 | gene %in% apa_reg) %>%
        dplyr::select(-cancer_count) #left known apa regulators and cancer type>1 gene
#----------------------- pre load data
    imp<-result[,c("cpg","gene","ipa","cancer_type","ACME_Estimate","ACME_CI_lower","ACME_CI_Upper","Prop.Mediated_Estimate","Prop_CI_lower","Prop_CI_Upper")]# ACME大小、稳定性；中介比例大小、稳定性
    # 1. ACME normalization
    imp$ACME_Estimate_norm <- rescale(log1p(abs(imp$ACME_Estimate)), to = c(0, 1))
    # 2. ACME CI normalization
    imp$ACME_CI_width <- imp$ACME_CI_Upper - imp$ACME_CI_lower
    imp$ACME_CI_width <- log1p( 1 / (imp$ACME_CI_width + 1e-6))
    imp$ACME_CI_norm <- rescale(abs(imp$ACME_CI_width), to = c(0, 1))
    # 3. Prop normalization
    imp$Prop_norm <- rescale(log1p(abs(imp$Prop.Mediated_Estimate)), to = c(0, 1))
    # 4. Prop CI normalization
    imp$Prop_CI_width <- imp$Prop_CI_Upper - imp$Prop_CI_lower
    imp$Prop_CI_width <- log1p( 1 / (imp$Prop_CI_width + 1e-6))
    imp$Prop_CI_norm <- rescale(abs(imp$Prop_CI_width), to = c(0, 1))
    setDT(imp)

    df<- data.frame(matrix(0,ncol=4,nrow=length(unique(result$gene)),dimnames=list(unique(result$gene),c("pairs","ipas","cgs","cancers"))))
    for(i in unique(result$gene)){
        temp <- result[result$gene==i,]
        df[i,1] <- length(temp$ipa)
        df[i,2] <- length(unique(temp$ipa))
        df[i,3] <- length(unique(temp$cpg))
        df[i,4] <- length(unique(temp$cancer_type))
    }
    df <- as.data.frame(df)
    df$gene <- rownames(df)
    df[, c("pairs", "ipas", "cgs","cancers")] <- lapply(df[, c("pairs", "ipas", "cgs","cancers")], as.integer)# 中介路径数、IPA事件种类、受调甲基位点种类、涵盖癌症种类
    df[1:4]<-lapply(df[1:4],function(x)rescale(log1p(x),to=c(0,1)))
    setDT(df)
#----------------------- start calculation
    # pathscore importance
    w=test_path[1,]
    imp_copy <- copy(imp)
    imp_copy[, score := w$w1 * ACME_Estimate_norm +
                    w$w2 * ACME_CI_norm +
                    w$w3 * Prop_norm +
                    w$w4 * Prop_CI_norm]
    imp2 <- imp_copy[, .(
        mean_score = mean(score),
        median_score = median(score),
        max_score = max(score),
        sum_score = sum(score)
    ), by = gene]  #gene score, use median_score  
    # structurescore importance
    w=test_str[1,]
    df_copy <- copy(df)
    df2 <- df_copy[, str := w$w1 * pairs + w$w2 * ipas + w$w3 * cgs + w$w4 * cancers]
    res <- merge(df2,imp2,by="gene")
    res <- res[,c("gene","str","median_score")]
    colnames(res) <- c("gene","str_score","path_score")
save(res,file="./importance_mediation_gene2.RData")