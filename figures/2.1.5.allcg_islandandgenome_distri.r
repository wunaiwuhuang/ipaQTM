library(ggplot2)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
source("./2.1.0.distribute_enrich_fun.r") #load compute_enrichment_results function
################################ cgisland feature draw ################################
#------------------------------- calculate enrichment coefficients
    load("./cgisland_feature_data.Rdata")
    cis_result <- compute_enrichment_results(cis_df,"cis")
    cis_result <- subset(cis_result, P_value_Chi <= 0.05)
    trans_result <- compute_enrichment_results(trans_df,"trans")
    trans_result <- subset(trans_result, P_value_Chi <= 0.05)
#-------------------------------- draw box plot
    # load data
    total_results <- rbind(cis_result, trans_result)
    region_order <- c("island", "shore", "shelf", "opensea")
    total_results$Region <- factor(total_results$Region, levels = region_order)
    # calculate point sizes
    total_results$logP <- -log10(total_results$P_value_Chi)
    total_results$size <- cut(total_results$logP, 
                            breaks = c(-Inf, 2, 10, 100, Inf), 
                            labels = c(1, 2, 3, 4),
                            right = TRUE) # must be TRUE,otherwise inf won't be included
    total_results$size <- as.numeric(as.character(total_results$size))    
    # load color
    source("./0.color.r")
    cancer_colors <- generate_cancer_colors(unique(total_results$Cancer))
    # draw points and boxes
    p <- ggplot(total_results, aes(x = Region, y = Odds_Ratio, fill = Type)) +
        geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +  # 让 box 之间有间距
        geom_point(data = total_results,
                aes(color = Cancer, size = size,group = Type), # 强制按 Type 分组散点,以便确定散点能稳定排成一列
                position = position_dodge(width = 0.8)) +  # 确保散点在 box 里排成一列
        scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
        scale_color_manual(values = cancer_colors,name = "Cancer Type") +
        scale_size_continuous(name = "-log10(P-value)",breaks = c(1, 2, 3, 4),labels = c("<2", "2~10", "10~100", ">100"),range = c(1, 6)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
        labs(title = "Genomic Feature", x = "Region", y = "Fold Enrichment")+
        coord_cartesian(ylim = c(0, 2)) # only show 0-2
    ggsave(filename="./2.1.5.cgisland_feature_anno.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)

################################ geneomic feature draw ################################
# -------------------------------- calculate enrichment coefficients
    load("./genomic_feature_data.Rdata")
    cis_result <- compute_enrichment_results(cis_df,"cis")
    cis_result <- subset(cis_result, P_value_Chi <= 0.05)
    trans_result <- compute_enrichment_results(trans_df,"trans")
    trans_result <- subset(trans_result, P_value_Chi <= 0.05)
# --------------------------------------- draw box plot
    # load data
    total_results <- rbind(cis_result, trans_result)
    region_order <- c("1stExon", "OtherExon", "1stIntron", "OtherIntron","3'UTR","5'UTR","IGR","TSS1000","TSS2000","TSS3000")
    total_results$Region <- factor(total_results$Region, levels = region_order)
    # calculate point sizes
    total_results$logP <- -log10(total_results$P_value_Chi)
    total_results$size <- cut(total_results$logP, 
                            breaks = c(-Inf, 2, 10, 100, Inf), 
                            labels = c(1, 2, 3, 4),
                            right = TRUE) # must be TRUE,otherwise inf won't be included
    total_results$size <- as.numeric(as.character(total_results$size))    
    # load color
    source("./0.color.r")
    cancer_colors <- generate_cancer_colors(unique(total_results$Cancer))
    # draw points and boxes
    p <- ggplot(total_results, aes(x = Region, y = Odds_Ratio, fill = Type)) +
        geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +  # 让 box 之间有间距
        geom_point(data = total_results,
                aes(color = Cancer, size = size,group = Type), # 强制按 Type 分组散点,以便确定散点能稳定排成一列
                position = position_dodge(width = 0.8)) +  # 确保散点在 box 里排成一列
        scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
        scale_color_manual(values = cancer_colors,name = "Cancer Type") +
        scale_size_continuous(name = "-log10(P-value)",breaks = c(1, 2, 3, 4),labels = c("<2", "2~10", "10~100", ">100"),range = c(1, 6)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
        labs(title = "Genomic Feature", x = "Region", y = "Fold Enrichment")+
        coord_cartesian(ylim = c(0, 4)) # only show 0-2
    ggsave(filename="./2.1.5.genomic_feature_anno.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)
