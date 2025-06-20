library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

# 获取符合模式的文件列表，并提取 cancer 名称
mDNA_files <- list.files("../data/", pattern = "_mdnause.txt$")
cancer <- gsub("_mdnause.txt$", "", mDNA_files)  

# 初始化数据框
all_data <- data.frame()

# 循环读取所有癌症类型的cis和trans数据
for (c in cancer) {
  # 读取cis数据
  cis_file <- paste0("../cis_sig/", c, "_cisfilter.txt")
  if (file.exists(cis_file)) {
    cis_data <- fread(cis_file, header = TRUE, stringsAsFactors = FALSE)
    cis_data$type <- "cis"
    cis_data$c <- c
    all_data <- rbind(all_data, cis_data)
  }
  
  # 读取trans数据
  trans_file <- paste0("../trans_sig/", c, "_transfilter.txt")
  if (file.exists(trans_file)) {
    trans_data <- fread(trans_file, header = TRUE, stringsAsFactors = FALSE)
    trans_data$type <- "trans"
    trans_data$c <- c
    all_data <- rbind(all_data, trans_data)
  }
}

# 创建绘图数据
plot_data <- all_data %>%
  filter(!is.na(Fraction_genotype))

# 绘制直方图
hist_plot <- ggplot(plot_data, aes(x = Fraction_genotype, fill = type)) +
  geom_histogram(aes(y = ..count..), binwidth = 0.001, position = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
  labs(x = NULL, y = "Number of APA events", 
       title = "Distribution of APA variation explained by CpG sites") +
  theme_minimal() +
  theme(legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  xlim(0, 0.5)  # 限制横轴范围

# 箱线图
box_plot <- ggplot(plot_data, aes(x = Fraction_genotype, y = type, fill = type)) +
  geom_boxplot(alpha = 0.7,outlier.shape = NA) +
  scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
  labs(x = "Fraction of genotype", y = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  xlim(0, 0.5)  # 保持横轴范围一致

# 组合图
combined_plot <- plot_grid(hist_plot, box_plot, ncol = 1, align = "v", rel_heights = c(3.5, 1))

# 保存图形
ggsave("1.4.ciaandtrans_explain_proportion.pdf", combined_plot, width = 10, height = 8, dpi = 300)
