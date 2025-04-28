library(data.table)
library(ggplot2)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

# 获取符合模式的文件列表，并提取 cancer 名称
files <- list.files("../data/", pattern = "_mdnause.txt$")
cancer <- gsub("_mdnause.txt$", "", files)  

# 1. 读取所有 ../cis_sig/cancer_cisfilter_noex.txt 文件并建立数据框 df
df_list <- lapply(cancer, function(c) {
  file_path <- paste0("../cis_sig/", c, "_cisfilter_noex.txt")
  if (file.exists(file_path)) {
    data <- fread(file_path, select = c("snps", "gene", "pvalue"))
    data[, qtm := snps]  # 赋值 qtm 列
    return(data[, .(qtm, gene, pvalue)])
  }
  return(NULL)
})
df <- rbindlist(df_list, use.names = TRUE, fill = TRUE)

# 2. 读取所有 ../cancer_IPAloc.txt 文件并建立数据框 df_ipasite
df_ipasite_list <- lapply(cancer, function(c) {
  file_path <- paste0("../data/", c, "_ipaloc.txt")
  if (file.exists(file_path)) {
    data <- fread(file_path, select = c("id", "start","end"))
    colnames(data)[1] <- "gene"
    return(unique(data))
  }
  return(NULL)
})
df_ipasite <- unique(rbindlist(df_ipasite_list, use.names = TRUE, fill = TRUE))
df_ipasite$chain <- str_extract(df_ipasite$gene, ":[+-]:")
df_ipasite$chain <- gsub(":", "", df_ipasite$chain)
#if chain is +,then use end,else use start
df_ipasite$pos <- ifelse(df_ipasite$chain == "+",df_ipasite$end,df_ipasite$start)
df_ipasite<-subset(df_ipasite, select = -c(start, end))

# 3. 读取 ../mDNAloc.txt 文件并建立数据框 df_qtmsite
if (file.exists("../00.probe_map_filtered.RData")) {
  load("../00.probe_map_filtered.RData")
  df_qtmsite <- probe_map[, c("id", "start","end")]
  df_qtmsite$pos <- (df_qtmsite$start+df_qtmsite$end)/2
  setnames(df_qtmsite, c("id", "pos"), c("qtm", "qtmsite"))
} else {
  df_qtmsite <- data.table(qtm = character(), qtmsite = numeric())
}
df_qtmsite <- subset(df_qtmsite, select = -c(start, end))

# 4. 将 df 的 qtmsite, ipasite 列对应填入数据
df <- merge(df, df_qtmsite, by = "qtm", all.x = TRUE)
df <- merge(df, df_ipasite, by = "gene", all.x = TRUE)

# 5. calculate istance and draw figure
# if +, than pos-qtmsite, if -, than qtmsite-pos
df$distance <- ifelse(df$chain == "+", df$pos - df$qtmsite, df$qtmsite - df$pos)
df$log<- -log10(df$pvalue)

p <- ggplot(df, aes(x = distance, y = log)) +
  geom_point(size = 0.01) +  # 设置点的大小
  labs(x = "Distance", y = "-(Log10 pvalue)", title = "Scatter Plot of Distance vs Log") +
  theme_minimal()
ggsave("1.3.qtm_ipa_sitedistance.pdf", plot = p, width = 4, height = 4, dpi = 300)