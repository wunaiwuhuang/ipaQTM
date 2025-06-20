setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

med_path <- "/data1/wuguojia/data/IPA_QTM_tcga/mediation/"
load(paste0(med_path,"./all_sig_result.RData"))
load("../apa_reg/apa_regulators.RData")
#-------------------------------------------  original ---------------------
#------------------------ by cancer_type
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
        labs(title = paste(toupper(data_type), " count by cancer"),
            x = "Cancer type",
            y = "count") +
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
#------------------------ by cancer_number
    df<- data.frame(matrix(0,ncol=4,nrow=length(unique(result$gene)),dimnames=list(unique(result$gene),c("pairs","ipas","cgs","cancers"))))
    for(i in unique(result$gene)){
        temp <- result[result$gene==i,]
        df[i,1] <- length(temp$ipa)
        df[i,2] <- length(unique(temp$ipa))
        df[i,3] <- length(unique(temp$cpg))
        df[i,4] <- length(unique(temp$cancer_type))
    }
    df <- as.data.frame(df)
    df$gene <- rownames(df)
    df[, c("pairs", "ipas", "cgs","cancers")] <- lapply(df[, c("pairs", "ipas", "cgs","cancers")], as.integer)
    # all gene
    df_pie<-df %>% count(cancers) %>% arrange(desc(cancers)) %>% mutate(perc=n/sum(n)*100,label=paste0(cancers," (",round(perc,1),"%)"),label = factor(label, levels = label))
    p1<- ggplot(df_pie, aes(x = "", y = n, fill = label)) +
        geom_col(width = 1, color = "white") +
        coord_polar(theta = "y") +
        theme_void() +
        labs(fill = "Cancer number") +
        ggtitle(paste0("Distribution of all gene (",sum(df_pie$n),")")) +
        theme(plot.title = element_text(hjust = 0.5),legend.position = "right",) +
        scale_fill_brewer(palette = "Set3")
    ggsave(filename = "./3.1.cancer_num_allgene.pdf", plot = p1, width = 6, height = 6)
    # known regulators
    df_pie<-df[df$gene %in% apa_reg,] %>% count(cancers) %>% arrange(desc(cancers)) %>% mutate(perc=n/sum(n)*100,label=paste0(cancers," (",round(perc,1),"%)"),label = factor(label, levels = label))
    p2<- ggplot(df_pie, aes(x = "", y = n, fill = label)) +
        geom_col(width = 1, color = "white") +
        coord_polar(theta = "y") +
        theme_void() +
        labs(fill = "Cancer number") +
        ggtitle(paste0("Distribution of known reg (",sum(df_pie$n),")")) +
        theme(plot.title = element_text(hjust = 0.5),legend.position = "right",) +
        scale_fill_brewer(palette = "Set3")
    ggsave(filename = "./3.1.cancer_num_knownreg.pdf", plot = p2, width = 6, height = 6)
#------------------------------------------- after adjustment --------------
#-------------------------- filter gene
    load(paste0(med_path,"./all_sig_result.RData"))
    load("../apa_reg/apa_regulators.RData")
    result <- result %>%
        group_by(gene) %>%
        mutate(cancer_count = n_distinct(cancer_type)) %>%
        ungroup() %>%
        filter(cancer_count > 1 | gene %in% apa_reg) %>%
        dplyr::select(-cancer_count) #left known apa regulators and cancer type>1 gene
    result <- result[!grepl("^ZNF", result$gene), ] # remove ZNF family genes
#------------------------ by cancer_type
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
        labs(title = paste(toupper(data_type), " count by cancer"),
            x = "Cancer type",
            y = "count") +
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
    output_file <- paste0("./3.1.mediation_statistic_afteradj.pdf")
    ggsave(output_file, combined_plot, width = 12, height = 8)
#------------------------ by cancer_number
    df<- data.frame(matrix(0,ncol=4,nrow=length(unique(result$gene)),dimnames=list(unique(result$gene),c("pairs","ipas","cgs","cancers"))))
    for(i in unique(result$gene)){
        temp <- result[result$gene==i,]
        df[i,1] <- length(temp$ipa)
        df[i,2] <- length(unique(temp$ipa))
        df[i,3] <- length(unique(temp$cpg))
        df[i,4] <- length(unique(temp$cancer_type))
    }
    df <- as.data.frame(df)
    df$gene <- rownames(df)
    df[, c("pairs", "ipas", "cgs","cancers")] <- lapply(df[, c("pairs", "ipas", "cgs","cancers")], as.integer)
    # all gene
    df_pie<-df %>% count(cancers) %>% arrange(desc(cancers)) %>% mutate(perc=n/sum(n)*100,label=paste0(cancers," (",round(perc,1),"%)"),label = factor(label, levels = label))
    p1<- ggplot(df_pie, aes(x = "", y = n, fill = label)) +
        geom_col(width = 1, color = "white") +
        coord_polar(theta = "y") +
        theme_void() +
        labs(fill = "Cancer number") +
        ggtitle(paste0("Distribution of all gene (",sum(df_pie$n),")")) +
        theme(plot.title = element_text(hjust = 0.5),legend.position = "right",) +
        scale_fill_brewer(palette = "Set3")
    ggsave(filename = "./3.1.cancer_num_allgene_afteradj.pdf", plot = p1, width = 6, height = 6)
    # known regulators
    df_pie<-df[df$gene %in% apa_reg,] %>% count(cancers) %>% arrange(desc(cancers)) %>% mutate(perc=n/sum(n)*100,label=paste0(cancers," (",round(perc,1),"%)"),label = factor(label, levels = label))
    p2<- ggplot(df_pie, aes(x = "", y = n, fill = label)) +
        geom_col(width = 1, color = "white") +
        coord_polar(theta = "y") +
        theme_void() +
        labs(fill = "Cancer number") +
        ggtitle(paste0("Distribution of known reg (",sum(df_pie$n),")")) +
        theme(plot.title = element_text(hjust = 0.5),legend.position = "right",) +
        scale_fill_brewer(palette = "Set3")
    ggsave(filename = "./3.1.cancer_num_knownreg_afteradj.pdf", plot = p2, width = 6, height = 6)
