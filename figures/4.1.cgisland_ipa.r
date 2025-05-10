library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
load("./cgisland_feature_anno.Rdata")
files<-list.files(path="../cis_sig/",pattern="*cisfilter_noex.txt",full.names=F)
cancers<-gsub("_cisfilter_noex.txt","",files)
df<-data.frame(matrix(0,nrow=length(cancers),ncol=4,dimnames=list(cancers,c("island","shore","shelf","opensea"))))
anno<-as.data.frame(anno)
for(file in files){
    cancer<-gsub("_cisfilter_noex.txt","",file)
    temp<-fread(paste0("../cis_sig/",file),header = T,sep = "\t")
    temp<-temp[,c("snps","gene")]%>%as.data.frame()%>%unique()
    temp$cancer<-cancer
    temp<-left_join(temp,anno,by=c("snps"="probeID"))
    df[cancer,"island"]<-dim(temp[temp$CGIposition=="Island",])[1]
    df[cancer,"shore"]<-dim(temp[temp$CGIposition=="Shore",])[1]
    df[cancer,"shelf"]<-dim(temp[temp$CGIposition=="Shelf",])[1]
    df[cancer,"opensea"]<-dim(temp[temp$CGIposition=="Opensea",])[1]
}
# stack plot
    df_long <- df %>%
        rownames_to_column("cancer") %>%
        pivot_longer(cols = -cancer, names_to = "region", values_to = "value") %>%
        group_by(cancer) %>%
        mutate(prop = value / sum(value)) %>%
        ungroup()
    df_long$region <- factor(df_long$region, levels = c("island", "shore", "shelf", "opensea"))#reorder the levels
    p<-ggplot(df_long, aes(x = cancer, y = prop, fill = region)) +
        geom_bar(stat = "identity") +
        labs(x = "Cancer Type", y = "Proportion", fill = "Region") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
    ggsave(filename="./4.1.cgisland_ipa_stack.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)
# box plot
    df_prop <- df / rowSums(df)  # 每行归一化，得到概率（占比）
    df_long <- pivot_longer(as.data.frame(df_prop), cols = everything(), names_to = "Region", values_to = "Proportion")
    df_long$Cancer <- rep(rownames(df), each = ncol(df))
    df_long$Region <- factor(df_long$Region, levels = c("island", "shore", "shelf", "opensea"))
    df_long$Proportion <- as.numeric(df_long$Proportion)
    my_comparisons <- combn(levels(df_long$Region), 2, simplify = FALSE)
    # 加载颜色
    source("./0.color.r")
    cancer_colors <- generate_cancer_colors(unique(df_long$Cancer))
    # 绘图
    p <- ggplot(df_long, aes(x = Region, y = Proportion)) +
    geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8), outlier.shape = NA) + 
    geom_point(aes(color = Cancer, group = Region),  # 按 Region 分组，确保每类单独排布
                size = 2,
                position = position_dodge(width = 0.8)) +  # 确保散点在 box 里排成一列
    stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",  
                     label = "p.signif") +
    scale_color_manual(values = cancer_colors, name = "Cancer Type") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    ) +
    labs(
        title = "Distribution of Genomic Regions across Cancers",
        x = "Genomic Region",
        y = "Proportion"
    ) +
    coord_cartesian(ylim = c(0, 0.8))
    ggsave(filename="./4.1.cgisland_ipa_box.pdf",plot=p,device="pdf",width=6,height=6,units="in",dpi=300)
