setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

med_path <- "/data1/wuguojia/data/IPA_QTM_tcga/mediation/"
load(paste0(med_path,"./all_sig_result.RData"))

cancer <- unique(result$cancer_type)
df <- data.frame(matrix(0,nrow=length(cancer), ncol=4, dimnames=list(cancer, c("pair", "cg","gene", "ipa"))))
for(i in cancer){
    tmp <- subset(result,result$cancer_type==i)
    df[i, "pair"] <- dim(tmp)[1]
    df[i, "cg"] <- length(unique(tmp$cpg))
    df[i, "gene"] <- length(unique(tmp$gene))
    df[i, "ipa"] <- length(unique(tmp$ipa))
}

draw_df <- df
data_types <- colnames(draw_df)
draw_df$cancer <- rownames(draw_df)
draw_list <- list()
for (data_type in data_types) {
    plot_df <- draw_df %>%
    dplyr::select(cancer, !!sym(data_type))
    plot_df$cancer <- factor(plot_df$cancer, levels = plot_df$cancer)
    p <- ggplot(plot_df, aes(x = cancer, y = !!sym(data_type))) +
    geom_bar(stat = "identity", fill = "#1f77b4") +
    theme_minimal(base_size = 14) +
    labs(title = paste(toupper(data_type), "pair count by cancer"),
        x = "Cancer type",
        y = "Pair count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #save the plot
    draw_list[[data_type]] <- p
}
# Combine the plots into a single figure
combined_plot <- plot_grid(
    draw_list$pair, 
    draw_list$cg, 
    draw_list$gene, 
    draw_list$ipa, 
    ncol = 2,  
    rel_widths = c(1, 1),  
    align = "h"
)
# Save the combined plot
output_file <- paste0("./3.1.mediation_statistic.pdf")
ggsave(output_file, combined_plot, width = 12, height = 8)