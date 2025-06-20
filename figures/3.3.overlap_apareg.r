library(data.table)
library(UpSetR)
library(ggplot2)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

load("../apa_reg/apa_regulators.RData")
med_path <- "/data1/wuguojia/data/IPA_QTM_tcga/mediation/"
load(paste0(med_path,"./all_sig_result.RData"))
# all cancer
    gene <- unique(result$gene)
    gene_sets <- list(
        Gene = gene,
        TranDB = tranddb,
        Harmonizome = harmonizome
    )
    pdf("./3.3.upset_allcancer.pdf", width = 8, height = 6)
    upset(
        fromList(gene_sets),  # 将列表转换为适合UpSetR的格式
        order.by = "freq",    # 按交集频率排序
        sets = c("Gene", "TranDB", "Harmonizome"),  # 指定集合顺序
        mainbar.y.label = "Number of Shared Genes",  # Y轴标签
        sets.x.label = "Total Genes per Set",        # X轴标签
        text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2)  # 调整文本大小
    )
    dev.off()
# BRCA
    brca <- subset(result,result$cancer_type=="BRCA")
    gene <- unique(brca$gene)
    gene_sets <- list(
        Gene = gene,
        TranDB = tranddb,
        Harmonizome = harmonizome
    )
    pdf("./3.3.upset_brca.pdf", width = 8, height = 6)
    upset(
        fromList(gene_sets),  # 将列表转换为适合UpSetR的格式
        order.by = "freq",    # 按交集频率排序
        sets = c("Gene", "TranDB", "Harmonizome"),  # 指定集合顺序
        mainbar.y.label = "Number of Shared Genes",  # Y轴标签
        sets.x.label = "Total Genes per Set",        # X轴标签
        text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2)  # 调整文本大小
    )
    dev.off()