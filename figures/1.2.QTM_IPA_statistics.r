library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

files <- list.files("../data/", pattern = "_mdnause.txt$")
cancer <- gsub("_mdnause.txt$", "", files)  

################################ statistics
  result_df <- data.frame(
    cis_qtm = numeric(length(cancer)),
    cis_ipa = numeric(length(cancer)),
    cis_pair = numeric(length(cancer)),
    trans_qtm = numeric(length(cancer)),
    trans_ipa = numeric(length(cancer)),
    trans_pair = numeric(length(cancer)),
    row.names = cancer)

  for (c in cancer) {
    cis_file <- paste0("../cis_sig/", c, "_cisfilter_noex.txt")
    if (file.exists(cis_file)) {
      cis_data <- read.table(cis_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      result_df[c, "cis_qtm"] <- length(unique(cis_data$snps))
      result_df[c, "cis_ipa"] <- length(unique(cis_data$gene))
      result_df[c, "cis_pair"] <- dim(cis_data)[1]
    } #cis files
    trans_file <- paste0("../trans_sig/", c, "_transfilter_noex.txt")
    if (file.exists(trans_file)) {
      trans_data <- read.table(trans_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      result_df[c, "trans_qtm"] <- length(unique(trans_data$snps))
      result_df[c, "trans_ipa"] <- length(unique(trans_data$gene))
      result_df[c,"trans_pair"] <- dim(trans_data)[1]
    } #trans files
  }
################################ draw figure: QTM and IPA unique numbers
  draw_df <- result_df[,c("cis_qtm","cis_ipa","trans_qtm","trans_ipa")]
  data_types <- c("cis", "trans")
  for (data_type in data_types) {
    # QTM data
    qtm_col <- paste0(data_type, "_qtm")
    qtm_df <- draw_df %>%
    dplyr::select(!!sym(qtm_col)) %>%
    mutate(cancer = rownames(draw_df)) %>%
    arrange(!!sym(qtm_col))  # reorder
    # reorder the data frame by QTM values
    qtm_df$cancer <- factor(qtm_df$cancer, levels = qtm_df$cancer[order(-qtm_df[[qtm_col]])])

    qtm_plot <- ggplot(qtm_df, aes(x = !!sym(qtm_col) / 10000, y = cancer)) +
    geom_bar(stat = "identity", fill = "#B22222") +
    scale_x_reverse(name = "QTM (10^4)") +  # X 轴反向
    labs(y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10), axis.ticks.y = element_blank())

    # ipa data
    ipa_col <- paste0(data_type, "_ipa")
    ipa_df <- draw_df %>%
    dplyr::select(!!sym(ipa_col),!!sym(qtm_col)) %>%
    mutate(cancer = rownames(draw_df)) %>%
    arrange(!!sym(qtm_col))  # ensure the order of cancer types is the same as in qtm_df
    # reorder the data frame by QTM values
    ipa_df$cancer <- factor(ipa_df$cancer, levels = qtm_df$cancer[order(-qtm_df[[qtm_col]])])

    ipa_plot <- ggplot(ipa_df, aes(x = !!sym(ipa_col) / 1000, y = cancer)) +
    geom_bar(stat = "identity", fill = "#003366") +
    scale_x_continuous(name = "IPA (10^3)") +
    labs(y = NULL) +
    theme_minimal()

    # combined plot
    final_plot <- plot_grid(
    qtm_plot, 
    ipa_plot, 
    ncol = 2,  
    rel_widths = c(1, 1),  
    align = "h"
    )
    # save the plot
    output_file <- paste0("1.2.",data_type, "_unique.pdf")
    ggsave(output_file, final_plot, width = 6, height = 6)
  }
################################ draw figure: QTM and IPA pairs
  draw_df <- result_df[,c("cis_pair","trans_pair")]
  data_types <- c("cis", "trans")
  draw_df$cancer <- rownames(draw_df)
  for (data_type in data_types) {
    qtm_col <- paste0(data_type, "_pair")
    plot_df <- draw_df %>%
      dplyr::select(cancer, !!sym(qtm_col)) %>%
      arrange(desc(!!sym(qtm_col)))
    plot_df$cancer <- factor(plot_df$cancer, levels = plot_df$cancer)
    
    p <- ggplot(plot_df, aes(x = cancer, y = !!sym(qtm_col))) +
      geom_bar(stat = "identity", fill = ifelse(data_type == "cis", "#1f77b4", "#ff7f0e")) +
      theme_minimal(base_size = 14) +
      labs(title = paste(toupper(data_type), "pair count by cancer"),
          x = "Cancer type",
          y = "Pair count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #save the plot
    output_file <- paste0("1.2.",data_type, "_pair.pdf")
    ggsave(output_file, p, width = 6, height = 6)
  }