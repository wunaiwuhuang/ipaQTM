# 加载必要的库
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(data.table)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

files <- list.files("../data/", pattern = "_mdnause.txt")
cancer <- gsub("_mdnause.txt", "", files)  
################################ statistics
    result_df <- data.frame(
        cpg = numeric(length(cancer)),
        ipa = numeric(length(cancer)),
        sample = numeric(length(cancer)),
        row.names = cancer
    )
    for (c in cancer) {
        cpg_file <- paste0("../data/", c, "_mdnause.txt")
        if (file.exists(cpg_file)) {
            print(paste0("dealing with file: ",cpg_file)) 
            cpg_data <- fread(cpg_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
            result_df[c, "cpg"] <- length(unique(cpg_data$id))
        }
        ipa_file <- paste0("../data/", c, "_ipause.txt")
        if (file.exists(ipa_file)) {
            print(paste0("dealing with file: ",ipa_file))
            ipa_data <- fread(ipa_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
            result_df[c, "ipa"] <- length(unique(ipa_data$id))
            result_df[c, "sample"] <- length(colnames(ipa_data)) - 1  # 减去第一列的 id 列
        }
    }
################################ sample and IPA ################################
    draw_df <- result_df
    # 处理 sample 数据（横轴朝左）
    sample_df <- draw_df %>%
    select(sample) %>%
    mutate(cancer = rownames(draw_df)) %>%
    arrange(-sample)  # 按 sample 数值从大到小排序
    # 更新因子水平，使其按 sample 排序
    sample_df$cancer <- factor(sample_df$cancer, levels = sample_df$cancer)
    sample_plot <- ggplot(sample_df, aes(x = sample / 100, y = cancer)) +
        geom_bar(stat = "identity", fill = "#B22222") +
        scale_x_reverse(name = "sample (10^2)") +  # X 轴反向
        labs(y = NULL) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 10), axis.ticks.y = element_blank())
    # 处理 IPA 数据（横轴朝右）
    ipa_df <- draw_df %>%
    select(ipa, sample) %>%
    mutate(cancer = rownames(draw_df)) %>%
    arrange(-sample)  # 确保 IPA 也按照 sample 的顺序排列
    # 更新因子水平
    ipa_df$cancer <- factor(ipa_df$cancer, levels = sample_df$cancer)
    ipa_plot <- ggplot(ipa_df, aes(x = ipa / 1000, y = cancer)) +
        geom_bar(stat = "identity", fill = "#003366") +
        scale_x_continuous(name = "IPA (10^3)") +
        labs(y = NULL) +
        theme_minimal()
    # 使用 cowplot::plot_grid() 拼接两张图
    final_plot <- plot_grid(
        sample_plot, 
        ipa_plot, 
        ncol = 2,  # 2 列布局
        rel_widths = c(1, 1),  # 两张图等宽
        align = "h"
    )
    # 保存最终图到 PDF 文件
    ggsave("1.1.sample_IPA_statistics.pdf", final_plot, width = 6, height = 6)
################################ sample and CPG ################################
    draw_df <- result_df
    # 处理 sample 数据（横轴朝左）
    sample_df <- draw_df %>%
    select(sample) %>%
    mutate(cancer = rownames(draw_df)) %>%
    arrange(-sample)  # 按 sample 数值从大到小排序
    # 更新因子水平，使其按 sample 排序
    sample_df$cancer <- factor(sample_df$cancer, levels = sample_df$cancer)
    sample_plot <- ggplot(sample_df, aes(x = sample / 100, y = cancer)) +
        geom_bar(stat = "identity", fill = "#B22222") +
        scale_x_reverse(name = "sample (10^2)") +  # X 轴反向
        labs(y = NULL) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 10), axis.ticks.y = element_blank())
    # 处理 cpg 数据（横轴朝右）
    cpg_df <- draw_df %>%
    select(sample, cpg) %>%
    mutate(cancer = rownames(draw_df)) %>%
    arrange(-sample)  # 确保 cpg 也按照 CPG 的顺序排列
    # 更新因子水平
    cpg_df$cancer <- factor(cpg_df$cancer, levels = cpg_df$cancer)
    cpg_plot <- ggplot(cpg_df, aes(x = cpg / 10000, y = cancer)) +
        geom_bar(stat = "identity", fill = "#003366") +
        scale_x_continuous(name = "cpg (10^4)") +
        labs(y = NULL) +
        theme_minimal()
    # 使用 cowplot::plot_grid() 拼接两张图
    final_plot <- plot_grid(
        sample_plot, 
        cpg_plot, 
        ncol = 2,  # 2 列布局
        rel_widths = c(1, 1),  # 两张图等宽
        align = "h"
    )
    # 保存最终图到 PDF 文件
    ggsave("1.1.sample_CPG_statistics.pdf", final_plot, width = 6, height = 6)
