# to statistic correlation
library(data.table)
library(ggplot2)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
files <- list.files("../data/", pattern = "_mdnause.txt")
cancers <- gsub("_mdnause.txt", "", files)

result <- data.frame(matrix(0,nrow=length(cancers),ncol=3))
colnames(result) <- c("pos","med","neg")
rownames(result) <- cancers

# 2. 读取所有 ../cancers_IPAloc.txt 文件并建立数据框 df_ipasite
df_ipasite_list <- lapply(cancers, function(c) {
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


################################ cis
    for(cancer in cancers){
        cat(paste0("dealing with file: ",cancer,"_cisfilter_noex.txt\n"))
        cis <- fread(paste0("../cis_sig/",cancer,"_cisfilter_noex.txt"))
        pos <- nrow(subset(cis,cis$r >= 0.1))
        med <- nrow(subset(cis,cis$r < 0.1 & cis$r > -0.1))
        neg <- nrow(subset(cis,cis$r <= -0.1))
        result[cancer,"pos"] <- pos
        result[cancer,"med"] <- med
        result[cancer,"neg"] <- neg
    }
    result_cis <- result
    result_cis$per <- result_cis$pos/result_cis$neg
    #draw
    df <- result_cis
    df$Cancer <- rownames(df)
    df_long <- df %>%
    dplyr::select(Cancer, pos, neg) %>%
    pivot_longer(cols = c(pos, neg), names_to = "Type", values_to = "Count")
    scale_factor <- max(df$pos + df$neg) / max(df$per)
    #start drawing
    p <-ggplot() +
        geom_bar(data = df_long,aes(x = Cancer, y = Count, fill = Type),stat = "identity") +
        scale_fill_manual(values = c(pos = "red", neg = "blue"),name = "Type") +
        #line
        geom_line(data = df,aes(x = Cancer, y = per * scale_factor, group = 1),color = "black", size = 0.8) +
        geom_point(data = df,aes(x = Cancer, y = per * scale_factor),color = "black", size = 2) +
        # add horizontal line at per = 1 (scaled)
        geom_hline(yintercept = 1 * scale_factor, color = "darkgreen", linetype = "dashed", size = 1) +
        #axis
        scale_y_continuous(name = "Count (pos + neg)",sec.axis = sec_axis(~ . / scale_factor, name = "pos/neg")) +
        labs(title = "Cis Associations by Cancer",x = "Cancer Type") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y.right = element_text(color = "black"),
            axis.title.y.left = element_text(color = "black")
        )
    ggsave(filename="./2.4.correlation_all_cis.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)

################################ trans
    for(cancer in cancers){
        cat(paste0("dealing with file: ",cancer,"_transfilter_noex.txt\n"))
        trans <- fread(paste0("../trans_sig/",cancer,"_transfilter_noex.txt"))
        pos <- nrow(subset(trans,trans$r >= 0.1))
        med <- nrow(subset(trans,trans$r < 0.1 & trans$r > -0.1))
        neg <- nrow(subset(trans,trans$r <= -0.1))
        result[cancer,"pos"] <- pos
        result[cancer,"med"] <- med
        result[cancer,"neg"] <- neg
    }
    result_trans <- result
    result_trans$per <- result_trans$pos/result_trans$neg
    #draw
    df <- result_trans
    df$Cancer <- rownames(df)
    df_long <- df %>%
    dplyr::select(Cancer, pos, neg) %>%
    pivot_longer(cols = c(pos, neg), names_to = "Type", values_to = "Count")
    scale_factor <- max(df$pos + df$neg) / max(df$per)
    #start drawing
    p <-ggplot() +
        geom_bar(data = df_long,aes(x = Cancer, y = Count, fill = Type),stat = "identity") +
        scale_fill_manual(values = c(pos = "red", neg = "blue"),name = "Type") +
        #line
        geom_line(data = df,aes(x = Cancer, y = per * scale_factor, group = 1),color = "black", size = 0.8) +
        geom_point(data = df,aes(x = Cancer, y = per * scale_factor),color = "black", size = 2) +
        # add horizontal line at per = 1 (scaled)
        geom_hline(yintercept = 1 * scale_factor, color = "darkgreen", linetype = "dashed", size = 1) +
        #axis
        scale_y_continuous(name = "Count (pos + neg)",sec.axis = sec_axis(~ . / scale_factor, name = "pos/neg")) +
        labs(title = "Cis Associations by Cancer",x = "Cancer Type") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y.right = element_text(color = "black"),
            axis.title.y.left = element_text(color = "black")
        )
    ggsave(filename="./2.4.correlation_all_trans.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)

################################ cis-100bp,1000bp,10000bp,100000bp
    thresholds <- c(100,1000,10000,100000)
    for(threshold in thresholds){
        for(cancer in cancers){
            cat(paste0("dealing with file: ",cancer,"_cisfilter_noex.txt\n"))
            cis <- fread(paste0("../cis_sig/",cancer,"_cisfilter_noex.txt"))
            cis <- left_join(cis,df_ipasite,by="gene")
            cis <- left_join(cis,df_qtmsite,by=c("snps"="qtm"))
            cis$dist <- abs(cis$pos - cis$qtmsite)
            cis <- subset(cis, cis$dist <= threshold)
            pos <- nrow(subset(cis,cis$r >= 0.1))
            med <- nrow(subset(cis,cis$r < 0.1 & cis$r > -0.1))
            neg <- nrow(subset(cis,cis$r <= -0.1))
            result[cancer,"pos"] <- pos
            result[cancer,"med"] <- med
            result[cancer,"neg"] <- neg
        }
        result_cis <- result
        result_cis$per <- result_cis$pos/result_cis$neg
        #draw
        df <- result_cis
        df$Cancer <- rownames(df)
        df_long <- df %>%
        dplyr::select(Cancer, pos, neg) %>%
        pivot_longer(cols = c(pos, neg), names_to = "Type", values_to = "Count")
        scale_factor <- max(df$pos + df$neg) / max(df$per)
        #start drawing
        p <-ggplot() +
            geom_bar(data = df_long,aes(x = Cancer, y = Count, fill = Type),stat = "identity") +
            scale_fill_manual(values = c(pos = "red", neg = "blue"),name = "Type") +
            #line
            geom_line(data = df,aes(x = Cancer, y = per * scale_factor, group = 1),color = "black", size = 0.8) +
            geom_point(data = df,aes(x = Cancer, y = per * scale_factor),color = "black", size = 2) +
            # add horizontal line at per = 1 (scaled)
            geom_hline(yintercept = 1 * scale_factor, color = "darkgreen", linetype = "dashed", size = 1) +
            #axis
            scale_y_continuous(name = "Count (pos + neg)",sec.axis = sec_axis(~ . / scale_factor, name = "pos/neg")) +
            labs(title = "Cis Associations by Cancer",x = "Cancer Type") +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title.y.right = element_text(color = "black"),
                axis.title.y.left = element_text(color = "black")
            )
        ggsave(filename=paste0("./2.4.correlation_all_cis_",threshold,"bp.pdf"),plot=p,device="pdf",width=10,height=6,units="in",dpi=300)        
    }
###################################################################################