library(data.table)
library(UpSetR)
library(ggplot2)
library(dplyr)
library(scales)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

load("../mediation/all_sig_result.RData")
load("./importance_mediation_gene2.RData")
load("../apa_reg/apa_regulators.RData")
#--------------------------- plot functions
    #select top gene
    select_topgene <- function(data, column_name,threshold) {
        topgene <- data %>%
            arrange(desc(.data[[column_name]])) %>%
            head(threshold)
        genes <- unique(topgene$gene)
        return(genes)
    }
    # plot topgene gene in ring
    plot_topgene_ring <- function(data, column_name, title_suffix = "",threshold) {
        topgene <- data %>%
            arrange(desc(.data[[column_name]])) %>%
            head(threshold) %>%
            mutate(
                is_apa = gene %in% apa_reg, #overlap with known apa-reg
                scaled_value = ( .data[[column_name]] - min(.data[[column_name]])) / (max(.data[[column_name]]) - min(.data[[column_name]])) * 0.5 + 1,
                id = row_number(),
                label = paste0(sprintf("%.2f", .data[[column_name]]), "  ", gene),
                angle = 90 - 360 * (id - 0.5) / threshold,  # 固定100个基因的角度计算
                hjust = ifelse(angle < -90, 1, 0),
                angle = ifelse(angle < -90, angle + 180, angle),
                label_y = scaled_value + 0.08  # 每个标签略高于柱顶
            )
        p<-ggplot(topgene, aes(x = factor(id), y = scaled_value ,fill = is_apa)) +
            geom_col(color = "black", width = 1, linewidth = 0.1) +  # 使用geom_col替代geom_bar
            geom_text(
            aes(y = label_y, label = label, angle = angle, hjust = hjust),
            size = 3.5, color = "black"
            ) +
            coord_polar(start = 0) +
            theme_void() + 
            scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray70")) +
            theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
            ) +
            labs(title = paste("Top ",threshold," Genes for", column_name, title_suffix))
        return(p)
    }
#-------------------------------------------  original ---------------------
#--------------------------- all cancer: str_score
    # ring plot
    threshold <- 100
    p_allring <- plot_topgene_ring(res, "str_score", " all cancer",threshold)
    ggsave("./3.5.topgene_ring_all_str.pdf", p_allring, width = 8, height = 8, dpi = 300)
    # get top 100 genes intersect with apa_regulators
    threshold <- 100
    top_all <-select_topgene(res,"str_score",threshold)
    # upset
    gene_sets <- list(
        TREND_DB=tranddb,
        HARMONIZOME=harmonizome,
        `topgene_cg-GENE-ipa`=top_all
    )
    pdf("./3.5.topgene_upset_all_str.pdf", width = 8, height = 6)
    upset(
        fromList(gene_sets),  # 将列表转换为适合UpSetR的格式
        order.by = "freq",    # 按交集频率排序
        sets = c("TREND_DB","HARMONIZOME","topgene_cg-GENE-ipa"),  # 指定集合顺序
        mainbar.y.label = "Number of Shared Genes",  # Y轴标签
        sets.x.label = "Total Genes per Set",        # X轴标签
        text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2)  # 调整文本大小
    )
    dev.off()
#--------------------------- all cancer: path_score
    # ring plot
    threshold <- 100
    p_allring <- plot_topgene_ring(res, "path_score", " all cancer",threshold)
    ggsave("./3.5.topgene_ring_all_path.pdf", p_allring, width = 8, height = 8, dpi = 300)
    # get top 100 genes intersect with apa_regulators
    threshold <- 100
    top_all <-select_topgene(res,"path_score",threshold)
    # upset
    gene_sets <- list(
        TREND_DB=tranddb,
        HARMONIZOME=harmonizome,
        `topgene_cg-GENE-ipa`=top_all
    )
    pdf("./3.5.topgene_upset_all_path.pdf", width = 8, height = 6)
    upset(
        fromList(gene_sets),  # 将列表转换为适合UpSetR的格式
        order.by = "freq",    # 按交集频率排序
        sets = c("TREND_DB","HARMONIZOME","topgene_cg-GENE-ipa"),  # 指定集合顺序
        mainbar.y.label = "Number of Shared Genes",  # Y轴标签
        sets.x.label = "Total Genes per Set",        # X轴标签
        text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2)  # 调整文本大小
    )
    dev.off()
#--------------------------- all cancer: overlap str-path
    threshold <- 100
    path_top <- select_topgene(res,"path_score",threshold)
    str_top <-select_topgene(res,"str_score",threshold)
    # draw a venn diagram
    library(VennDiagram)
    venn.plot <- draw.triple.venn(
        area1 = length(path_top),
        area2 = length(str_top),
        area3 = length(apa_reg),
        n12 = length(intersect(path_top, str_top)),
        n23 = length(intersect(str_top, apa_reg)),
        n13 = length(intersect(path_top, apa_reg)),
        n123 = length(Reduce(intersect, list(path_top, str_top, apa_reg))),
        category = c("Path Score", "Str Score", "APA Reg"),
        fill = c("skyblue", "pink", "lightgreen"),
        lty = "blank",
        cex = 1.5,
        cat.cex = 1.2
    )
    grid.draw(venn.plot)
    intersect(path_top,str_top)
    intersect(path_top,apa_reg)
    intersect(str_top,apa_reg)
    intersect(intersect(path_top,str_top),apa_reg)
#------------------------------------------- after adjustment --------------
    load("../mediation/all_sig_result.RData")
    load("./importance_mediation_gene2.RData")
    load("../apa_reg/apa_regulators.RData")
    result <- result %>%
        group_by(gene) %>%
        mutate(cancer_count = n_distinct(cancer_type)) %>%
        ungroup() %>%
        filter(cancer_count > 1 | gene %in% apa_reg) %>%
        dplyr::select(-cancer_count) #left known apa regulators and cancer type>1 gene
    result <- result[!grepl("^ZNF", result$gene), ] # remove ZNF family genes
    res <- res[!grepl("^ZNF", res$gene), ] # remove ZNF family genes
#--------------------------- all cancer: str_score
    # ring plot
    threshold <- 100
    p_allring <- plot_topgene_ring(res, "str_score", " all cancer",threshold)
    ggsave("./3.5.topgene_ring_all_str_afteradj.pdf", p_allring, width = 8, height = 8, dpi = 300)
    # get top 100 genes intersect with apa_regulators
    threshold <- 100
    top_all <-select_topgene(res,"str_score",threshold)
    # upset
    gene_sets <- list(
        TREND_DB=tranddb,
        HARMONIZOME=harmonizome,
        `topgene_cg-GENE-ipa`=top_all
    )
    pdf("./3.5.topgene_upset_all_str_afteradj.pdf", width = 8, height = 6)
    upset(
        fromList(gene_sets),  # 将列表转换为适合UpSetR的格式
        order.by = "freq",    # 按交集频率排序
        sets = c("TREND_DB","HARMONIZOME","topgene_cg-GENE-ipa"),  # 指定集合顺序
        mainbar.y.label = "Number of Shared Genes",  # Y轴标签
        sets.x.label = "Total Genes per Set",        # X轴标签
        text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2)  # 调整文本大小
    )
    dev.off()
#--------------------------- all cancer: path_score
    # ring plot
    threshold <- 100
    p_allring <- plot_topgene_ring(res, "path_score", " all cancer",threshold)
    ggsave("./3.5.topgene_ring_all_path_afteradj.pdf", p_allring, width = 8, height = 8, dpi = 300)
    # get top 100 genes intersect with apa_regulators
    threshold <- 100
    top_all <-select_topgene(res,"path_score",threshold)
    # upset
    gene_sets <- list(
        TREND_DB=tranddb,
        HARMONIZOME=harmonizome,
        `topgene_cg-GENE-ipa`=top_all
    )
    pdf("./3.5.topgene_upset_all_path_afteradj.pdf", width = 8, height = 6)
    upset(
        fromList(gene_sets),  # 将列表转换为适合UpSetR的格式
        order.by = "freq",    # 按交集频率排序
        sets = c("TREND_DB","HARMONIZOME","topgene_cg-GENE-ipa"),  # 指定集合顺序
        mainbar.y.label = "Number of Shared Genes",  # Y轴标签
        sets.x.label = "Total Genes per Set",        # X轴标签
        text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2)  # 调整文本大小
    )
    dev.off()
#--------------------------- all cancer: overlap str-path
    threshold <- 100
    path_top <- select_topgene(res,"path_score",threshold)
    str_top <-select_topgene(res,"str_score",threshold)
    # draw a venn diagram
    library(VennDiagram)
    venn.plot <- draw.triple.venn(
        area1 = length(path_top),
        area2 = length(str_top),
        area3 = length(apa_reg),
        n12 = length(intersect(path_top, str_top)),
        n23 = length(intersect(str_top, apa_reg)),
        n13 = length(intersect(path_top, apa_reg)),
        n123 = length(Reduce(intersect, list(path_top, str_top, apa_reg))),
        category = c("Path Score", "Str Score", "APA Reg"),
        fill = c("skyblue", "pink", "lightgreen"),
        lty = "blank",
        cex = 1.5,
        cat.cex = 1.2
    )
    grid.draw(venn.plot)
    intersect(path_top,str_top)
    intersect(path_top,apa_reg)
    intersect(str_top,apa_reg)
    intersect(intersect(path_top,str_top),apa_reg)