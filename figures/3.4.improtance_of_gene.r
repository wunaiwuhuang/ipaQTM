library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(scales)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
load("../mediation/importance_coef_com.RData")
load("../mediation/all_sig_result.RData")
load("../apa_reg/apa_regulators.RData")
# select coefficient
    test_coef<-test_coef[test_coef$overlap_50==max(test_coef$overlap_50),]
    test_coef<-test_coef[test_coef$overlap_100==max(test_coef$overlap_100),]
    weight_grid<-test_coef[1:8]
    weight_grid<-as.data.table(weight_grid)
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
    df[1:4]<-lapply(df[1:4],function(x)rescale(log1p(x),to=c(1,3)))
    setDT(df)
#----------------------- start calculation
    # use 10 cores to calculate
    n_cores <- 10
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    test_coef <- foreach(i = 1:nrow(weight_grid), .combine = rbind, .packages = c("data.table", "dplyr", "scales") ) %dopar% {
        w <- weight_grid[i]
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
        df_copy <- copy(df)
        df2 <- df_copy[, str := w$w5 * pairs + w$w6 * ipas + w$w7 * cgs + w$w8 * cancers]
        res <- merge(df2,imp2,by="gene")
        res <- res[gene%in%unique(imp_copy[imp_copy$cancer_type=="BRCA",]$gene),]#select brca
        res <-res[,score:=median_score * str][order(-score)]
        genes_ranked <- res$gene
        top_genes <- intersect(genes_ranked[1:200], apa_reg)
        top_genes_str <- paste(top_genes, collapse = ";")
        data.frame(
            w1 = w$w1,
            w2 = w$w2,
            w3 = w$w3,
            w4 = w$w4,
            w5 = w$w5,
            w6 = w$w6,
            w7 = w$w7,
            w8 = w$w8,            
            overlap_50 = length(intersect(genes_ranked[1:50], apa_reg)),
            overlap_100 = length(intersect(genes_ranked[1:100], apa_reg)),
            overlap_150 = length(intersect(genes_ranked[1:150], apa_reg)),
            overlap_200 = length(top_genes),
            top_200gene = top_genes_str
        ) # number of overlapped known-apa regulators
    }
    # close c1
    stopCluster(cl)
# re-select coefficient
    test_coef<-test_coef[test_coef$overlap_100==max(test_coef$overlap_100),]
    test_coef<-test_coef[test_coef$overlap_200==max(test_coef$overlap_200),]
    test_coef<-test_coef[test_coef$overlap_50==max(test_coef$overlap_50),]
    # select final weight
    test_coef<-test_coef[test_coef$w2==test_coef$w6,]
    weight <- test_coef[1:8]%>%as.data.table()%>%.[3]
#----------------------- calculate final results
#> weight
#     w1  w2   w3   w4   w5  w6  w7   w8
#1: 0.14 0.5 0.24 0.12 0.22 0.5 0.1 0.18
    w <- weight
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
    df_copy <- copy(df)
    df2 <- df_copy[, str := w$w5 * pairs + w$w6 * ipas + w$w7 * cgs + w$w8 * cancers]
    res <- merge(df2,imp2,by="gene")
    res <-res[,score:=median_score * str][order(-score)]
save(res,file="./importance_mediation_gene.RData")