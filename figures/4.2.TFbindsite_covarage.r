library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
# ---------------------- calculate TF enrichment ----------------------
    load("../TFsite/TF_cg_coverage_matrix.RData")
    calc_enrichment <- function(df, mode = c("cis", "trans")) {
        mode <- match.arg(mode)
        sig_col <- paste0(mode, "_sig")
        nosig_col <- paste0(mode, "_nosig")
        # construct 2x2 matrix
        results <- lapply(1:nrow(df), function(i) {
            a <- df[[sig_col]][i]
            b <- df[[nosig_col]][i]
            c <- sum(df[[sig_col]]) - a
            d <- sum(df[[nosig_col]]) - b
            mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
            # fisher test
            fisher <- fisher.test(mat)
            data.frame(
                TF = rownames(df)[i],
                OR = fisher$estimate,
                p_value = fisher$p.value,
                sig_count = a,
                nosig_count = b
            )
        })
        result_df <- do.call(rbind, results)
        result_df$p_adj <- p.adjust(result_df$p_value, method = "BH")
        return(result_df)
    }

    all_enrich_results <- lapply(names(df_list), function(cancer) {
        df <- df_list[[cancer]]
        
        cis_result <- calc_enrichment(df, mode = "cis")
        cis_result$cancer <- cancer
        cis_result$type <- "cis"
        
        trans_result <- calc_enrichment(df, mode = "trans")
        trans_result$cancer <- cancer
        trans_result$type <- "trans"
        
        rbind(cis_result, trans_result)
    })

    final_results <- do.call(rbind, all_enrich_results)
    tf_coverage <- final_results
    save(tf_coverage, file = "./TF_enrichment_results.RData")
# ---------------------- visualisation ----------------------
    load("./TF_enrichment_results.RData")
    tf_sig <- tf_coverage[tf_coverage$p_adj <= 0.05, ]
    tf_sig$logP <- -log10(tf_sig$p_adj)
    tf_sig$size <- cut(tf_sig$logP, 
                        breaks = c(-Inf, 2, 10, 100, Inf), 
                        labels = c(1, 2, 3, 4),
                        right = TRUE) # must be TRUE,otherwise inf won't be included
    tf_sig$size <- as.numeric(as.character(tf_sig$size))
    # calculate top TFs
    threshold <- 30
    cis_tf_median <- tf_sig %>% filter(type == "cis") %>% group_by(TF) %>% summarize(median_OR = median(OR, na.rm = TRUE), .groups = "drop") %>% filter(median_OR >= 1)
    trans_tf_median <- tf_sig %>% filter(type == "trans") %>% group_by(TF) %>% summarize(median_OR = median(OR, na.rm = TRUE), .groups = "drop") %>% filter(median_OR >= 1)# first get median OR >= 1 TFs
    top_cis <- tf_sig %>% filter(type == "cis",TF %in% cis_tf_median$TF) %>% count(TF, sort = TRUE) %>% slice_max(n, n = threshold)
    top_trans <- tf_sig %>% filter(type == "trans",TF %in% trans_tf_median$TF) %>% count(TF, sort = TRUE) %>% slice_max(n, n = threshold)
    cis_data <- tf_sig %>% filter(type == "cis", TF %in% top_cis$TF)
    trans_data <- tf_sig %>% filter(type == "trans", TF %in% top_trans$TF)
    # plot function
    plot_tf_enrichment <- function(data, tf_set, type_label,ylim_range = c(1, 3)) {
        data$TF <- factor(data$TF, levels = tf_set$TF)  # order x axis
        # draw box and scatter plot
        ggplot(data, aes(x = TF, y = OR, fill = type)) +
            geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +
            geom_point(data = data,
                    aes(color = cancer, size = size,group = type), 
                    position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
            scale_color_manual(values = cancer_colors,name = "Cancer Type") +
            scale_size_continuous(name = "-log10(P-value)",
                                breaks = c(1, 2, 3, 4),
                                labels = c("<2", "2~10", "10~100", ">100"),
                                range = c(1, 8)) +
            coord_cartesian(ylim = ylim_range) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
            labs(title = paste0("Top TFs (", type_label, ")"),
                x = "Transcription Factor", y = "Odds Ratio")
    }
    # 加载颜色定义
    source("./0.color.r")
    cancer_colors <- generate_cancer_colors(unique(tf_sig$cancer))
    p_cis <- plot_tf_enrichment(cis_data, top_cis, "cis",ylim_range = c(1, 2.5))
    p_trans <- plot_tf_enrichment(trans_data, top_trans, "trans",ylim_range = c(1, 3))
    ggsave(filename="./4.2.tf_bindsite_cis_coverage.pdf",plot=p_cis,device="pdf",width=16,height=6,units="in",dpi=300)
    ggsave(filename="./4.2.tf_bindsite_trans_coverage.pdf",plot=p_trans,device="pdf",width=16,height=6,units="in",dpi=300)