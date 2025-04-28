
################################ cgisland feature data
    #download annotation files from https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz
    library(data.table)
    library(dplyr)
    library(rlang)

    setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

    # 获取符合模式的文件列表，并提取 cancer 名称
    mDNA_files <- list.files("../data/",pattern = "_mdnause.txt$")
    cancer <- gsub("_mdnause.txt$", "", mDNA_files)  

    # 定义基因组区域
    regions <- c("sig_Opensea","sig_Shelf","sig_Shore","sig_Island","nosig_Opensea","nosig_Shelf","nosig_Shore","nosig_Island")

    # 创建数据框
    cis_df <- data.frame(matrix(0, nrow = length(cancer), ncol = length(regions)))
    trans_df <- data.frame(matrix(0, nrow = length(cancer), ncol = length(regions)))
    rownames(cis_df) <- cancer
    rownames(trans_df) <- cancer
    colnames(cis_df) <- regions
    colnames(trans_df) <- regions

    #load annotation file
    anno<-fread("../00.HM450.hg38.manifest.gencode.v36.tsv")
    anno<-unique(anno[,c("probeID","CGIposition")])
    anno <- anno %>%
        mutate(CGIposition = case_when(
            CGIposition %in% c("N_Shore", "S_Shore") ~ "Shore",
            CGIposition %in% c("N_Shelf", "S_Shelf") ~ "Shelf",
            CGIposition == "Island" ~ "Island",
            is.na(CGIposition) ~ "Opensea",
            TRUE ~ CGIposition ))
    save(anno, file = "./cgisland_feature_anno.Rdata")

    load("./cgisland_feature_anno.Rdata")
    # 定义一个通用函数来处理 cis 和 trans
    process_data <- function(cancer_type, cg_path,mdna_path, df) {
        # 读取数据
        tem_cg <- fread(cg_path) %>% as.data.frame()
        tem_cg <- unique(tem_cg$snps)
        #match to annotation
        tem_cg <- anno[anno$probeID %in% tem_cg,]

        tem_mdna <- fread(mdna_path) %>% as.data.frame()
        tem_mdna <- unique(tem_mdna$id)
        #match to annotation
        tem_mdna <- anno[anno$probeID %in% tem_mdna,]

        # 定义 CGI 位置类别
        cgi_regions <- c("Island", "Shore", "Shelf", "Opensea")
        # 计算 CGI 位置的 sig 和 nosig 数量
        cgi_sig_count <- table(tem_cg$CGIposition)
        cgi_nosig_count <- table(tem_mdna$CGIposition)
        # 遍历 cgi_regions 计算 sig 和 nosig 数量
        for (region in cgi_regions) {
            sig_count <- cgi_sig_count[region] %||% 0
            nosig_count <- (cgi_nosig_count[region] %||% 0) - sig_count
            # 更新 df
            df[cancer_type, paste0("sig_", region)] <- sig_count
            df[cancer_type, paste0("nosig_", region)] <- nosig_count
        }
        return(df)
    }

    # 计算 cis 数据
    for (c in cancer) {
        cg_path <- paste0("../cis_sig/", c, "_cisfilter_noex.txt")
        mdna_path <- paste0("../data/", c, "_mdnause.txt")
        if (file.exists(cg_path)) {
            print(paste("Processing CIS file:", cg_path))
            cis_df <- process_data(c, cg_path,mdna_path, cis_df)
        }
    }

    # 计算 trans 数据
    for (c in cancer) {
        cg_path <- paste0("../trans_sig/", c, "_transfilter_noex.txt")
        mdna_path <- paste0("../data/", c, "_mdnause.txt")
        if (file.exists(cg_path)) {
            print(paste("Processing TRANS file:", cg_path))
            trans_df <- process_data(c, cg_path,mdna_path, trans_df)
        }
    }

    save(cis_df, trans_df, file = "cgisland_feature_data.Rdata")
