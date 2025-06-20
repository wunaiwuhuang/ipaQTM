library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(scales)
library(dplyr)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/mediation")
load("./importance_coef_sep.RData")
load("./all_sig_result.RData")
load("../apa_reg/apa_regulators.RData")
# test_path select top coefficients
test_path<-test_path[test_path$overlap_50==max(test_path$overlap_50),]
test_path<-test_path[test_path$overlap_100>median(test_path$overlap_100),]
test_path<-test_path[test_path$overlap_200>=median(test_path$overlap_200),] 

# test_str select top coefficients
test_str<-test_str[test_str$overlap_50==max(test_str$overlap_50),]
test_str<-test_str[test_str$overlap_100==max(test_str$overlap_100),]
test_str<-test_str[test_str$w4 >= 0.28,]# assign cancer type more weight

# select final coefficients combination
path<-test_path[,1:4]
str<-test_str[,1:4]
steps <- seq(0.2, 0.8, by = 0.01)
coef <- as.data.table(expand.grid(a = steps, b = steps))
coef <- coef[abs(a + b - 1) < 1e-6]
weight_grid <- merge(path, str, by = NULL)
weight_grid <- merge(weight_grid, coef, by = NULL)
colnames(weight_grid)[1:8] <- paste0("w", 1:8)
weight_grid<-as.data.table(weight_grid)
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
    # use 10 cores to calculate
    n_cores <- 8
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
        res <-res[,score:=w$a*median_score + w$b*str][order(-score)]
        genes_ranked <- res$gene
        data.frame(
            a = w$a,
            b = w$b,
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
            overlap_200 = length(intersect(genes_ranked[1:200], apa_reg))
        ) # number of overlapped known-apa regulators
    }
    # close c1
    stopCluster(cl)
    save(test_coef,file="./importance_coef_com.RData")