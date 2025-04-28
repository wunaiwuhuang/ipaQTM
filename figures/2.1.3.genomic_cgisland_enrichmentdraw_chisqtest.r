################################ geneomic feature draw ################################
    library(ggplot2)
    setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
    load("./genomic_feature_data.Rdata")
################################ calculate enrichment coefficients
    compute_results <- function(df,label) {
        results <- list()
        for (cancer in rownames(df)) {
            cancer_results <- list()
            for (region in colnames(df)[1:(ncol(df)/2)]) { # 前一半是 sig_xxx
                sig_count <- df[cancer, region]
                nosig_count <- df[cancer, gsub("^sig_", "nosig_", region)]
                total_sig <- sum(df[cancer, 1:(ncol(df)/2)])
                total_nosig <- sum(df[cancer, ((ncol(df)/2)+1):ncol(df)])
                # 构建 2x2 表 并计算
                contingency_table <- matrix(c(sig_count, nosig_count, total_sig - sig_count, total_nosig - nosig_count), nrow = 2, byrow = TRUE)
                test_result <- fisher.test(contingency_table, alternative = "greater")
                chi_result <- chisq.test(contingency_table)
                g_test_stat <- function(table) {
                    observed <- as.vector(table)
                    expected <- rowSums(table) %o% colSums(table) / sum(table)
                    G <- 2 * sum(observed * log(observed / expected))
                    p_value <- pchisq(G, df = 1, lower.tail = FALSE)
                    return(p_value)
                }
                g_test_p_value <- g_test_stat(contingency_table)
                # construct result table
                cancer_results[[region]] <- c(
                    "p_value_f" = test_result$p.value,
                    "p_value_chi" = chi_result$p.value,
                    "p_value_g" = g_test_p_value,
                    "odds_ratio" = (sig_count / (total_sig - sig_count)) / (nosig_count / (total_nosig - nosig_count))
                )
            }
            results[[cancer]] <- cancer_results
        }
        # combine all results
        df_results <- do.call(rbind, lapply(names(results), function(cancer) {
            data.frame(Cancer = cancer, 
                        Region = names(results[[cancer]]),
                        P_value_F = sapply(results[[cancer]], function(x) x["p_value_f"]),
                        P_value_Chi = sapply(results[[cancer]], function(x) x["p_value_chi"]),
                        P_value_G = sapply(results[[cancer]], function(x) x["p_value_g"]),
                        Odds_Ratio = sapply(results[[cancer]], function(x) x["odds_ratio"]),
                        Type = label)
        }))
        #remove unsuitable odds ratio
        df_results <- df_results[!is.na(df_results$Odds_Ratio) & is.finite(df_results$Odds_Ratio), ]
        return(df_results)
    }
    cis_result <- compute_results(cis_df,"cis")
    trans_result <- compute_results(trans_df,"trans")
################################ draw box plot
    # load data
    total_results <- rbind(cis_result, trans_result)
    total_results$Region <- gsub("sig_", "", total_results$Region)
    region_order <- c("1stExon", "OtherExon", "1stIntron", "OtherIntron","3'UTR","5'UTR","IGR","TSS1000","TSS2000","TSS3000")
    total_results$Region <- factor(total_results$Region, levels = region_order)
    # calculate point sizes
    total_results$logP <- -log10(total_results$P_value_Chi)
    total_results$size <- cut(total_results$logP, 
                            breaks = c(-Inf, 2, 10, 100, Inf), 
                            labels = c(1, 2, 3, 4),
                            right = FALSE)
    total_results$size <- as.numeric(as.character(total_results$size))    
    # load color
    source("./0.color.r")
    cancer_colors <- generate_cancer_colors(unique(total_results$Cancer))
    # draw points and boxes
    p <- ggplot(total_results, aes(x = Region, y = Odds_Ratio, fill = Type)) +
        geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +  # 让 box 之间有间距
        geom_point(data = total_results[total_results$logP > 2, ],
                aes(color = Cancer, size = size,group = Type), # 强制按 Type 分组散点,以便确定散点能稳定排成一列
                position = position_dodge(width = 0.8)) +  # 确保散点在 box 里排成一列
        scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
        scale_color_manual(values = cancer_colors,name = "Cancer Type") +
        scale_size_continuous(name = "-log10(P-value)",breaks = c(1, 2, 3, 4),labels = c("<2", "2~10", "10~100", ">100"),range = c(1, 6)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
        labs(title = "Genomic Feature", x = "Region", y = "Fold Enrichment")+
        coord_cartesian(ylim = c(0, 2)) # only show 0-2
    ggsave(filename="./2.1.1.genomic_feature_anno.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)
################################ cgisland feature draw ################################
################################ calculate enrichment coefficients
    library(ggplot2)
    setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")
    load("./cgisland_feature_data.Rdata")
    cis_result <- compute_results(cis_df,"cis")
    trans_result <- compute_results(trans_df,"trans")
################################ draw box plot
    # load data
    total_results <- rbind(cis_result, trans_result)
    total_results$Region <- gsub("sig_", "", total_results$Region)
    region_order <- c("Island", "Shore", "Shelf", "Opensea")
    total_results$Region <- factor(total_results$Region, levels = region_order)
    # calculate point sizes
    total_results$logP <- -log10(total_results$P_value_Chi)
    total_results$size <- cut(total_results$logP, 
                            breaks = c(-Inf, 2, 10, 100, Inf), 
                            labels = c(1, 2, 3, 4),
                            right = FALSE)
    total_results$size <- as.numeric(as.character(total_results$size))    
    # load color
    source("./0.color.r")
    cancer_colors <- generate_cancer_colors(unique(total_results$Cancer))
    # draw points and boxes
    p <- ggplot(total_results, aes(x = Region, y = Odds_Ratio, fill = Type)) +
        geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +  # 让 box 之间有间距
        geom_point(data = total_results[total_results$logP > 2, ],
                aes(color = Cancer, size = size,group = Type), # 强制按 Type 分组散点,以便确定散点能稳定排成一列
                position = position_dodge(width = 0.8)) +  # 确保散点在 box 里排成一列
        scale_fill_manual(values = c("cis" = "#FFB6C1", "trans" = "#ADD8E6")) +
        scale_color_manual(values = cancer_colors,name = "Cancer Type") +
        scale_size_continuous(name = "-log10(P-value)",breaks = c(1, 2, 3, 4),labels = c("<2", "2~10", "10~100", ">100"),range = c(1, 6)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
        labs(title = "Genomic Feature", x = "Region", y = "Fold Enrichment")+
        coord_cartesian(ylim = c(0, 2)) # only show 0-2
    ggsave(filename="./2.1.2.cgisland_feature_anno.pdf",plot=p,device="pdf",width=10,height=6,units="in",dpi=300)
