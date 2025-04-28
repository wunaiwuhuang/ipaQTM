################################ genomic feature data
    # this part i choose to use chipseeker to complete
    library(data.table)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(dplyr)
    library(rlang)
    library(stringr)

    setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

    # 获取符合模式的文件列表，并提取 cancer 名称
    mDNA_files <- list.files("../data/",pattern = "_mdnause.txt$")
    cancer <- gsub("_mdnause.txt$", "", mDNA_files)  

    # 定义基因组区域
    regions <- c("sig_1stExon", "sig_1stIntron", "sig_OtherExon", "sig_OtherIntron",
    "sig_3'UTR", "sig_5'UTR", "sig_IGR", "sig_TSS1000", "sig_TSS2000", "sig_TSS3000",
    "nosig_1stExon", "nosig_1stIntron", "nosig_OtherExon", "nosig_OtherIntron",
    "nosig_3'UTR", "nosig_5'UTR", "nosig_IGR", "nosig_TSS1000", "nosig_TSS2000", "nosig_TSS3000")

    # 创建数据框
    cis_df <- data.frame(matrix(0, nrow = length(cancer), ncol = length(regions)))
    trans_df <- data.frame(matrix(0, nrow = length(cancer), ncol = length(regions)))
    rownames(cis_df) <- cancer
    rownames(trans_df) <- cancer
    colnames(cis_df) <- regions
    colnames(trans_df) <- regions

    # 读取mDNAloc.txt
    load("../00.probe_map_filtered.RData")
    mDNAloc <- probe_map[, .(chr, start,end, id)]
    fwrite(mDNAloc, "./temp.bed", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
    # 使用ChIPseeker注释
    peakAnno <- annotatePeak("./temp.bed",tssRegion=c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,annoDb="org.Hs.eg.db")
    # 计算各类注释的数量
    anno<-data.frame(id=peakAnno@anno$V4,feature=peakAnno@anno$annotation)
        anno$feature <- case_when(
        anno$feature == "Distal Intergenic" ~"IGR" ,
        anno$feature == "Downstream (<=300bp)" ~"IGR" ,#只有18个，扔到这里吧
        anno$feature == "3' UTR" ~ "3'UTR",
        anno$feature == "5' UTR" ~"5'UTR" ,
        anno$feature == "Promoter (<=1kb)" ~"TSS1000" ,
        anno$feature == "Promoter (1-2kb)" ~"TSS2000" ,
        anno$feature == "Promoter (2-3kb)" ~"TSS3000" ,
        TRUE ~ anno$feature
        )
        # 处理包含 "exon 1 of" 和 "intron 1 of" 的情况
        anno$feature <- ifelse(str_detect(anno$feature, "exon 1 of"), "1stExon", anno$feature)
        anno$feature <- ifelse(str_detect(anno$feature, "exon") & !str_detect(anno$feature, "exon 1 of"), "OtherExon", anno$feature)
        anno$feature <- ifelse(str_detect(anno$feature, "intron 1 of"), "1stIntron", anno$feature)
        anno$feature <- ifelse(str_detect(anno$feature, "intron") & !str_detect(anno$feature, "intron 1 of"), "OtherIntron", anno$feature)
    save(anno,file="./genomic_feature_anno.Rdata")
    file.remove("./temp.bed")

    load("./genomic_feature_anno.Rdata") #load anno
    # 定义一个通用函数来处理 cis 和 trans
    process_data <- function(cancer_type, cg_path,mdna_path, df) {
        tem_cg <- fread(cg_path) %>% as.data.frame()
        tem_cg <- unique(tem_cg$snps)
        #match to annotation
        tem_cg <- anno[anno$id %in% tem_cg,]

        tem_mdna <- fread(mdna_path) %>% as.data.frame()
        tem_mdna <- unique(tem_mdna$id)
        #match to annotation
        tem_mdna <- anno[anno$id %in% tem_mdna,]

        # 计算 sig 和 nosig 各类别的数量
        regions <- c("3'UTR", "5'UTR", "IGR", "TSS1000", "TSS2000", "TSS3000", 
                    "1stExon", "OtherExon", "1stIntron", "OtherIntron")
        # 遍历 regions 计算 sig 和 nosig 的数量
        for (region in regions) {
            # 计算 sig 数目
            sig_count <- sum(tem_cg$feature == region, na.rm = TRUE)
            # 计算 nosig 数目
            nosig_count <- sum(tem_mdna$feature == region, na.rm = TRUE) - sig_count
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

    # 保存结果
    save(cis_df, trans_df, file = "genomic_feature_data.Rdata")