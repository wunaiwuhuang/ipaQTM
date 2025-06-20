# to statistic correlation
library(data.table)
library(ggplot2)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
files <- list.files("../data/", pattern = "_mdnause.txt")
cancers <- gsub("_mdnause.txt", "", files)

result <- data.frame(matrix(0,nrow=length(cancers),ncol=3))
colnames(result) <- c("pos","med","neg")
rownames(result) <- cancers
result_up <- result
result_down <- result

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
        cis <- left_join(cis,df_ipasite,by="gene")
        cis <- left_join(cis,df_qtmsite,by=c("snps"="qtm"))
        cis$dist <- ifelse(cis$chain == "+", cis$pos - cis$qtmsite, cis$qtmsite - cis$pos)
        
        #positive means cg site in upstream of ipa site
        cis_up <- subset(cis, cis$dist >= 0)
        up_pos <- nrow(subset(cis_up,cis_up$r >= 0.1))
        up_med <- nrow(subset(cis_up,cis_up$r < 0.1 & cis_up$r > -0.1))
        up_neg <- nrow(subset(cis_up,cis_up$r <= -0.1))        
        #negative means cg site in downstream of ipa site
        cis_down <- subset(cis, cis$dist < 0)
        down_pos <- nrow(subset(cis_down,cis_down$r >= 0.1))
        down_med <- nrow(subset(cis_down,cis_down$r < 0.1 & cis_down$r > -0.1))
        down_neg <- nrow(subset(cis_down,cis_down$r <= -0.1))        
        
        result_up[cancer,"pos"] <- up_pos
        result_up[cancer,"med"] <- up_med
        result_up[cancer,"neg"] <- up_neg

        result_down[cancer,"pos"] <- down_pos
        result_down[cancer,"med"] <- down_med
        result_down[cancer,"neg"] <- down_neg        
    }
    result_cis_up <- result_up
    result_cis_up$per <- result_cis_up$pos/result_cis_up$neg
    result_cis_down <- result_down
    result_cis_down$per <- result_cis_down$pos/result_cis_down$neg

    #combine up and down
    draw_df <- data.frame(
        Cancer = rownames(result),
        up_num = result_cis_up$pos + result_cis_up$neg,
        down_num = result_cis_down$pos + result_cis_down$neg,
        up_per = result_cis_up$per,
        down_per = result_cis_down$per
    )
    # convert to long data
    bar_data <- draw_df %>%
    dplyr::select(Cancer, up_num, down_num) %>%
    pivot_longer(cols = c(up_num, down_num),
                names_to = "Direction",
                values_to = "Count") %>%
    mutate(Direction = factor(Direction, levels = c("down_num", "up_num")))
    line_data <- draw_df %>%
    dplyr::select(Cancer, up_per, down_per)
    #set scale
    scale_factor <- max(draw_df$up_num + draw_df$down_num) / max(draw_df$up_per,draw_df$down_per)
    #start drawing
    p <- ggplot() +
        geom_bar(data = bar_data, aes(x = Cancer, y = Count, fill = Direction),
                stat = "identity") +
        geom_line(data = line_data, aes(x = Cancer, y = up_per * scale_factor, group = 1),
                color = "red", size = 1) +
        geom_point(data = line_data, aes(x = Cancer, y = up_per * scale_factor),
                color = "red", size = 1.5) +
        geom_line(data = line_data, aes(x = Cancer, y = down_per * scale_factor, group = 1),
                color = "blue", size = 1) +
        geom_point(data = line_data, aes(x = Cancer, y = down_per * scale_factor),
                color = "blue", size = 1.5) +
        #set bar color
        scale_fill_manual(values = c("down_num" = "lightblue", "up_num" = "lightcoral")) +
        # add horizontal line at per = 1 (scaled)
        geom_hline(yintercept = 1 * scale_factor, color = "darkgreen", linetype = "dashed", size = 1) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y.right = element_text(color = "black"),
            axis.title.y.left = element_text(color = "black")
        )+
        labs(x = "Cancer", y = "Number of Events", fill = "Direction") +
        scale_y_continuous(
        name = "Number of Events",
        sec.axis = sec_axis(~ . / scale_factor, name = "pos/neg")
        )
    ggsave(filename="./2.4.correlation_updown_all.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)

################################ cis-100bp,1000bp,10000bp,100000bp
    thresholds <- c(100,1000,10000,100000)
    for(threshold in thresholds){
        for(cancer in cancers){
            cat(paste0("dealing with file: ",cancer,"_cisfilter_noex.txt\n"))
            cis <- fread(paste0("../cis_sig/",cancer,"_cisfilter_noex.txt"))
            cis <- left_join(cis,df_ipasite,by="gene")
            cis <- left_join(cis,df_qtmsite,by=c("snps"="qtm"))
            cis$dist <- ifelse(cis$chain == "+", cis$pos - cis$qtmsite, cis$qtmsite - cis$pos)
            cis <- subset(cis, abs(cis$dist) <= threshold)
            #positive means cg site in upstream of ipa site
            cis_up <- subset(cis, cis$dist >= 0)
            up_pos <- nrow(subset(cis_up,cis_up$r >= 0.1))
            up_med <- nrow(subset(cis_up,cis_up$r < 0.1 & cis_up$r > -0.1))
            up_neg <- nrow(subset(cis_up,cis_up$r <= -0.1))        
            #negative means cg site in downstream of ipa site
            cis_down <- subset(cis, cis$dist < 0)
            down_pos <- nrow(subset(cis_down,cis_down$r >= 0.1))
            down_med <- nrow(subset(cis_down,cis_down$r < 0.1 & cis_down$r > -0.1))
            down_neg <- nrow(subset(cis_down,cis_down$r <= -0.1))        
            
            result_up[cancer,"pos"] <- up_pos
            result_up[cancer,"med"] <- up_med
            result_up[cancer,"neg"] <- up_neg

            result_down[cancer,"pos"] <- down_pos
            result_down[cancer,"med"] <- down_med
            result_down[cancer,"neg"] <- down_neg        
        }
        result_cis_up <- result_up
        result_cis_up$per <- result_cis_up$pos/result_cis_up$neg
        result_cis_down <- result_down
        result_cis_down$per <- result_cis_down$pos/result_cis_down$neg

        #combine up and down
        draw_df <- data.frame(
            Cancer = rownames(result),
            up_num = result_cis_up$pos + result_cis_up$neg,
            down_num = result_cis_down$pos + result_cis_down$neg,
            up_per = result_cis_up$per,
            down_per = result_cis_down$per
        )
        # replace NaN and Inf
        draw_df$up_per[is.nan(draw_df$up_per)] <- 0
        if (any(is.infinite(draw_df$up_per))) {
        max_up_per <- max(draw_df$up_per[is.finite(draw_df$up_per)], na.rm = TRUE)
        draw_df$up_per[is.infinite(draw_df$up_per)] <- max_up_per
        }
        draw_df$down_per[is.nan(draw_df$down_per)] <- 0
        if (any(is.infinite(draw_df$down_per))) {
        max_down_per <- max(draw_df$down_per[is.finite(draw_df$down_per)], na.rm = TRUE)
        draw_df$down_per[is.infinite(draw_df$down_per)] <- max_down_per
        }
        # convert to long data
        bar_data <- draw_df %>%
        dplyr::select(Cancer, up_num, down_num) %>%
        pivot_longer(cols = c(up_num, down_num),
                    names_to = "Direction",
                    values_to = "Count") %>%
        mutate(Direction = factor(Direction, levels = c("down_num", "up_num")))
        line_data <- draw_df %>%
        dplyr::select(Cancer, up_per, down_per)
        #set scale
        scale_factor <- max(draw_df$up_num + draw_df$down_num) / max(draw_df$up_per,draw_df$down_per)
        #start drawing
        p <- ggplot() +
            geom_bar(data = bar_data, aes(x = Cancer, y = Count, fill = Direction),
                    stat = "identity") +
            geom_line(data = line_data, aes(x = Cancer, y = up_per * scale_factor, group = 1),
                    color = "red", size = 1) +
            geom_point(data = line_data, aes(x = Cancer, y = up_per * scale_factor),
                    color = "red", size = 1.5) +
            geom_line(data = line_data, aes(x = Cancer, y = down_per * scale_factor, group = 1),
                    color = "blue", size = 1) +
            geom_point(data = line_data, aes(x = Cancer, y = down_per * scale_factor),
                    color = "blue", size = 1.5) +
            #set bar color
            scale_fill_manual(values = c("down_num" = "lightblue", "up_num" = "lightcoral")) +
            # add horizontal line at per = 1 (scaled)
            geom_hline(yintercept = 1 * scale_factor, color = "darkgreen", linetype = "dashed", size = 1) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title.y.right = element_text(color = "black"),
                axis.title.y.left = element_text(color = "black")
            )+
            labs(x = "Cancer", y = "Number of Events", fill = "Direction") +
            scale_y_continuous(
            name = "Number of Events",
            sec.axis = sec_axis(~ . / scale_factor, name = "pos/neg")
            )
        ggsave(filename=paste0("./2.4.correlation_updown_",threshold,"bp.pdf"),plot=p,device="pdf",width=10,height=6,units="in",dpi=300)    
    }
###################################################################################