setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
source("./2.1.0.distribute_enrich_fun.r") #load compute_enrichment_results function,build_enrich_matrix function

load("/data1/wuguojia/data/IPA_QTM_tcga/mediation/all_sig_result.RData")
cancer <- unique(result$cancer_type)
####################################  cgisland location
# --------------------------------- calculate enrichment coefficients
    # prepare matrix
    load("./cgisland_feature_anno_2.Rdata") # we use anno_2, because it is more reliable
    load("../00.probe_map_filtered.RData")
    anno<-anno[anno$probeID %in% probe_map$id,] # we set all considered probe as background    
    cancers <- gsub("_mdnause.txt$", "", list.files("../data/",pattern = "_mdnause.txt$"))  
    types <- c("trans")
    sig_samcol <- "snps"
    bg <- anno
    bg_samcol <- "probeID"
    bg_regcol <- "CGIposition"
    region <- c("island","shelf","shore","opensea")
    # build results
    cis_df <- list()
    trans_df <- list()
    for (cancer in cancers) {
        for (type in types) {
            sig_file <- paste0("../", type, "_sig/", cancer, "_", type, "filter_noex.txt")
            if (file.exists(sig_file)) {
            message("Processing ", type, " for ", cancer)
            res <- build_enrich_matrix(
                cancer = cancer,
                type = type,
                sig_file = sig_file,
                sig_samcol = sig_samcol,
                bg = bg,
                bg_samcol = bg_samcol,
                bg_regcol = bg_regcol,
                region = region
            )
            if (type == "cis") {
                cis_df[[cancer]] <- res
            } else {
                trans_df[[cancer]] <- res
            }
            } else {
                message("File not found: ", sig_file)
            }
        }
    }
    trans_result <- compute_enrichment_results(trans_df,"trans")
    trans_result <- subset(trans_result, P_value_Chi <= 0.05)
#-------------------------------- visualization
    #draw box plot
    total_results <- trans_result
    region_order <- region
    total_results$Region <- factor(total_results$Region, levels = region_order)
    # calculate point sizes
    total_results$logP <- -log10(total_results$P_value_Chi)
    total_results$size <- cut(total_results$logP, 
                            breaks = c(-Inf, 1, 2, 5, Inf), 
                            labels = c(1, 2, 3, 4),
                            right = TRUE) # must be TRUE,otherwise inf won't be included
    total_results$size <- as.numeric(as.character(total_results$size))    
    # load color
    source("./0.color.r")
    cancer_colors <- generate_cancer_colors(unique(total_results$Cancer))
    # draw points and boxes
    p <- ggplot(total_results, aes(x = Region, y = Odds_Ratio, fill = Type)) +
        geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +  # 让 box 之间有间距
        geom_point(data = total_results,
                aes(color = Cancer, size = size,group = Type), # 强制按 Type 分组散点,以便确定散点能稳定排成一列
                position = position_dodge(width = 0.8)) +  # 确保散点在 box 里排成一列
        scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
        scale_color_manual(values = cancer_colors,name = "Cancer Type") +
        scale_size_continuous(name = "-log10(P-value)",breaks = c(1, 2, 3, 4),labels = c("<1", "1~2", "2~5", ">5"),range = c(1, 6)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
        labs(title = "Genomic Feature", x = "Region", y = "Fold Enrichment")+
        coord_cartesian(ylim = c(0, 3)) # only show 0-3
    ggsave(filename="./3.2.cgisland_cgloc.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)
####################################  genomic location
# --------------------------------- calculate enrichment coefficients
    # prepare matrix
    load("./genomic_feature_anno.Rdata") #load anno
    cancers <- gsub("_mdnause.txt$", "", list.files("../data/",pattern = "_mdnause.txt$"))  
    types <- c("trans")
    sig_samcol <- "snps"
    bg <- anno
    bg_samcol <- "id"
    bg_regcol <- "feature"
    region <- c("3'UTR", "5'UTR", "IGR", "TSS1000", "TSS2000", "TSS3000", "1stExon", "OtherExon", "1stIntron", "OtherIntron")
    # build results
    cis_df <- list()
    trans_df <- list()
    for (cancer in cancers) {
        for (type in types) {
            sig_file <- paste0("../", type, "_sig/", cancer, "_", type, "filter_noex.txt")
            if (file.exists(sig_file)) {
            message("Processing ", type, " for ", cancer)
            res <- build_enrich_matrix(
                cancer = cancer,
                type = type,
                sig_file = sig_file,
                sig_samcol = sig_samcol,
                bg = bg,
                bg_samcol = bg_samcol,
                bg_regcol = bg_regcol,
                region = region
            )
            if (type == "cis") {
                cis_df[[cancer]] <- res
            } else {
                trans_df[[cancer]] <- res
            }
            } else {
                message("File not found: ", sig_file)
            }
        }
    }
    trans_result <- compute_enrichment_results(trans_df,"trans")
    trans_result <- subset(trans_result, P_value_Chi <= 0.05)
#-------------------------------- visualization
    # draw box plot
    total_results <- trans_result
    region_order <- region
    total_results$Region <- factor(total_results$Region, levels = region_order)
    # calculate point sizes
    total_results$logP <- -log10(total_results$P_value_Chi)
    total_results$size <- cut(total_results$logP, 
                            breaks = c(-Inf, 1, 2, 5, Inf), 
                            labels = c(1, 2, 3, 4),
                            right = TRUE) # must be TRUE,otherwise inf won't be included
    total_results$size <- as.numeric(as.character(total_results$size))    
    # load color
    source("./0.color.r")
    cancer_colors <- generate_cancer_colors(unique(total_results$Cancer))
    # draw points and boxes
    p <- ggplot(total_results, aes(x = Region, y = Odds_Ratio, fill = Type)) +
        geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +  # 让 box 之间有间距
        geom_point(data = total_results,
                aes(color = Cancer, size = size,group = Type), # 强制按 Type 分组散点,以便确定散点能稳定排成一列
                position = position_dodge(width = 0.8)) +  # 确保散点在 box 里排成一列
        scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
        scale_color_manual(values = cancer_colors,name = "Cancer Type") +
        scale_size_continuous(name = "-log10(P-value)",breaks = c(1, 2, 3, 4),labels = c("<1", "1~2", "2~5", ">5"),range = c(1, 6)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
        labs(title = "Genomic Feature", x = "Region", y = "Fold Enrichment")+
        coord_cartesian(ylim = c(0, 5)) # only show 0-3
    ggsave(filename="./3.2.genomic_cgloc.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)

