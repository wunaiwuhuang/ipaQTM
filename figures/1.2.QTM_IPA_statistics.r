# 加载必要的库
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# 设置工作目录
setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

# 获取符合模式的文件列表，并提取 cancer 名称
files <- list.files("../data/", pattern = "_mdnause.txt$")
cancer <- gsub("_mdnause.txt$", "", files)  

################################ statistics
  # 创建一个空数据框，行名为cancer，列名为 cis_qtm, cis_ipa, trans_qtm, trans_ipa
  result_df <- data.frame(
    cis_qtm = numeric(length(cancer)),
    cis_ipa = numeric(length(cancer)),
    trans_qtm = numeric(length(cancer)),
    trans_ipa = numeric(length(cancer)),
    row.names = cancer
  )
  # 遍历每种cancer类型
  for (c in cancer) {
    # 读取 cis 数据
    cis_file <- paste0("../cis_sig/", c, "_cisfilter_noex.txt")
    if (file.exists(cis_file)) {
      cis_data <- read.table(cis_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      result_df[c, "cis_qtm"] <- length(unique(cis_data$snps))
      result_df[c, "cis_ipa"] <- length(unique(cis_data$gene))
    }
    # 读取 trans 数据
    trans_file <- paste0("../trans_sig/", c, "_transfilter_noex.txt")
    if (file.exists(trans_file)) {
      trans_data <- read.table(trans_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      result_df[c, "trans_qtm"] <- length(unique(trans_data$snps))
      result_df[c, "trans_ipa"] <- length(unique(trans_data$gene))
    }
  }
################################ draw figure
  #draw_df<-result_df
  draw_df <- result_df

  # 创建一个包含数据类型的列表
  data_types <- c("cis", "trans")
  # 循环处理 cis 和 trans 数据
  for (data_type in data_types) {
    # 处理 QTM 数据（横轴朝左）
    qtm_col <- paste0(data_type, "_qtm")
    qtm_df <- draw_df %>%
      select(!!sym(qtm_col)) %>%
      mutate(cancer = rownames(draw_df)) %>%
      arrange(!!sym(qtm_col))  # 按 QTM 数值从小到大排序

    # 更新因子水平，使其按 qtm 排序
    #qtm_df$cancer <- factor(qtm_df$cancer, levels = qtm_df$cancer)
    qtm_df$cancer <- factor(qtm_df$cancer, levels = qtm_df$cancer[order(-qtm_df[[qtm_col]])])

    qtm_plot <- ggplot(qtm_df, aes(x = !!sym(qtm_col) / 10000, y = cancer)) +
      geom_bar(stat = "identity", fill = "#B22222") +
      scale_x_reverse(name = "QTM (10^4)") +  # X 轴反向
      labs(y = NULL) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10), axis.ticks.y = element_blank())

    # 处理 IPA 数据（横轴朝右）
    ipa_col <- paste0(data_type, "_ipa")
    ipa_df <- draw_df %>%
      select(!!sym(ipa_col),!!sym(qtm_col)) %>%
      mutate(cancer = rownames(draw_df)) %>%
      arrange(!!sym(qtm_col))  # 确保 IPA 也按照 QTM 的顺序排列

    # 更新因子水平
    #ipa_df$cancer <- factor(ipa_df$cancer, levels = qtm_df$cancer)
    ipa_df$cancer <- factor(ipa_df$cancer, levels = qtm_df$cancer[order(-qtm_df[[qtm_col]])])


    ipa_plot <- ggplot(ipa_df, aes(x = !!sym(ipa_col) / 1000, y = cancer)) +
      geom_bar(stat = "identity", fill = "#003366") +
      scale_x_continuous(name = "IPA (10^3)") +
      labs(y = NULL) +
      theme_minimal()

    # 使用 cowplot::plot_grid() 拼接两张图
    final_plot <- plot_grid(
      qtm_plot, 
      ipa_plot, 
      ncol = 2,  # 2 列布局
      rel_widths = c(1, 1),  # 两张图等宽
      align = "h"
    )

    # 保存最终图到 PDF 文件
    output_file <- paste0("1.2.",data_type, "_plot.pdf")
    ggsave(output_file, final_plot, width = 6, height = 6)
  }

################################