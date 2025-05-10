setwd("/data1/wuguojia/data/IPA_QTM_tcga/mediation")
files <- list.files(path="./result/",pattern = ".*_mediate_all.RData",full.names = F)

result <- list()
for(file in files){
    cancer <- gsub("_mediate_all.RData", "", file)
    load(paste0("./result/",file))
    rtmp$FDR_ACME<-p.adjust(rtmp$ACME_pval,method="BH")
    result[[cancer]] <- rtmp
}
result <- do.call(rbind, result)
## 步骤1：去除缺失项和空值
result <- result[!is.na(result$ACME_Estimate) & result$ACME_Estimate != "", ]
## 步骤2：过滤假阳性
result <- result[result$FDR_ACME <= 0.05, ]
#result <- result[result$ACME_pval <= 0.05, ]
## 步骤3：ACME效应量下限=0.1
result <- result[abs(result$ACME_Estimate) >= 0.1, ]
## 步骤4: 去除ACME与ADE反向的路径
result <- result[sign(result$ACME_Estimate) == sign(result$ADE_Estimate), ]
## 步骤5: 保留Prop.Mediated_Estimate 0.1-1
result <- result[result$Prop.Mediated_Estimate >= 0 & result$Prop.Mediated_Estimate <= 1, ]
## 步骤6: Prop.Mediated_pval <= 0.05
result <- result[result$Prop.Mediated_pval <= 0.05, ]

save(result, file = "./all_sig_result.RData")