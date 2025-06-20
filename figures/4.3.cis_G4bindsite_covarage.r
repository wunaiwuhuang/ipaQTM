library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
# ---------------------- visualisation ----------------------
    load("../G4site/g4_enrichment_results.RData")
    g4_coverage<-rbindlist(results_list)
    g4_sig <- g4_coverage[g4_coverage$p_value <= 0.05, ]
    g4_sig$logP <- -log10(g4_sig$p_value)
    g4_sig$size <- cut(g4_sig$logP, 
                        breaks = c(-Inf, 2, 4, 8, Inf), 
                        labels = c(1, 2, 3, 4),
                        right = TRUE) # must be TRUE,otherwise inf won't be included
    g4_sig$size <- as.numeric(as.character(g4_sig$size))
    cis_data <- g4_sig %>% filter(type == "cis")
    trans_data <- g4_sig %>% filter(type == "trans")
    order <- paste0("chr",1:22)
    # plot function
    plot_chr_enrichment <- function(data, order, type_label,ylim_range = c(0, 5)) {
        data$chr <- factor(data$chr, levels = order)  # order x axis
        # draw box and scatter plot
        ggplot(data, aes(x = chr, y = OR, fill = type)) +
            geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +
            geom_point(data = data,
                    aes(color = cancer, size = size,group = type), 
                    position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
            scale_color_manual(values = cancer_colors,name = "Cancer Type") +
            scale_size_continuous(name = "-log10(P-value)",
                                breaks = c(1, 2, 3, 4),
                                labels = c("<2", "2~4", "4~8", ">8"),
                                range = c(1, 6)) +
            coord_cartesian(ylim = ylim_range) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
            labs(title = paste0("Top chrs (", type_label, ")"),
                x = "Transcription Factor", y = "Odds Ratio")
    }
    # 加载颜色定义
    source("./0.color.r")
    cancer_colors <- generate_cancer_colors(unique(g4_sig$cancer))
    p_cis <- plot_chr_enrichment(cis_data, order, "cis",ylim_range = c(0,2))
    p_trans <- plot_chr_enrichment(trans_data, order, "trans",ylim_range = c(0,4))
    ggsave(filename="./4.3.g4_bindsite_cis_coverage.pdf",plot=p_cis,device="pdf",width=16,height=6,units="in",dpi=300)
    ggsave(filename="./4.3.g4_bindsite_trans_coverage.pdf",plot=p_trans,device="pdf",width=16,height=6,units="in",dpi=300)