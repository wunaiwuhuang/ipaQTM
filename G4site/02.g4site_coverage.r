library(data.table)
library(dplyr)
library(GenomicRanges) # thounsand times faster than R
setwd("/data1/wuguojia/data/IPA_QTM_tcga/G4site")
# ------------------ background annotation ---------------------
    load("./00.g4site_foruse.RData")
    load("../00.probe_map_filtered.RData") 
    g4site<-out
    # use genomecranges to overlap
    probe_gr <- GRanges(seqnames = probe_map$chr,
                        ranges = IRanges(start = probe_map$start, end = probe_map$end),
                        id = probe_map$id)
    g4_gr <- GRanges(seqnames = g4site$chr,
                    ranges = IRanges(start = g4site$start, end = g4site$end),
                    id = g4site$id)
    #get annos
    hits <- findOverlaps(probe_gr, g4_gr)
    probe_ids <- mcols(probe_gr)$id[queryHits(hits)]
    overlap_count <- as.data.frame(table(probe_ids))
    colnames(overlap_count) <- c("id", "g4num")
    all <- data.frame(id = probe_map$id)
    anno <- left_join(all, overlap_count, by = "id")
    anno$g4num[is.na(anno$g4num)] <- 0
    anno <- left_join(anno,probe_map,by="id")
    anno <- anno[,c("id","chr","g4num")]
    save(anno, file = "./g4_cg_coverage_allanno.RData")
# ------------------ build matirx ---------------------
    cg_cis <- "/data1/wuguojia/data/IPA_QTM_tcga/cis_sig/"
    files_cis <- list.files(cg_cis, pattern = "cisfilter_noex.txt", full.names = F)
    cg_trans <- "/data1/wuguojia/data/IPA_QTM_tcga/trans_sig/"
    files_trans <- list.files(cg_trans, pattern = "transfilter_noex.txt", full.names = F)
    cancers<- gsub("_cisfilter_noex.txt", "", files_cis)
    load("./g4_cg_coverage_allanno.RData")
    results_list <- list()
    for(cancer in cancers){
        cat("Processing cancer:", cancer, "\n")
        cis <- fread(paste0(cg_cis, cancer, "_cisfilter_noex.txt"))
        trans <- fread(paste0(cg_trans, cancer, "_transfilter_noex.txt"))
        cis_sig_cpg <- unique(cis$snps)
        trans_sig_cpg <- unique(trans$snps)    
        # add sig/not sig status
        anno$cis <- ifelse(anno$id %in% cis_sig_cpg, TRUE, FALSE)
        anno$trans <- ifelse(anno$id %in% trans_sig_cpg, TRUE, FALSE)
        # g4 overlap status
        anno$g4_overlap <- anno$g4num > 0
        # calculate each chr respectively
        calc_chr_enrichment <- function(df, sig_col){
            chrs <- unique(df$chr)
            res <- lapply(chrs, function(chr){
                sub_df <- df[df$chr == chr, ]
                # 2x2 ：  
                #        g4 overlap   no g4 overlap
                # sig        a             b
                # non-sig    c             d
                a <- sum(sub_df[[sig_col]] & sub_df$g4_overlap)
                b <- sum(sub_df[[sig_col]] & (!sub_df$g4_overlap))
                c <- sum((!sub_df[[sig_col]]) & sub_df$g4_overlap)
                d <- sum((!sub_df[[sig_col]]) & (!sub_df$g4_overlap))
                # use fisher test to calculate enrichment
                mat <- matrix(c(a, b, c, d), nrow=2, byrow=TRUE)     
                ft <- fisher.test(mat)
                # results
                data.frame(
                    cancer = cancer,
                    chr = chr,
                    type = sig_col,
                    OR = ft$estimate,
                    p_value = ft$p.value,
                    g4_overlap_sig = a,
                    no_g4_overlap_sig = b,
                    g4_overlap_nosig = c,
                    no_g4_overlap_nosig = d
                )
            })
            do.call(rbind, res)
        }
        # calculate cis and trans
        cis_res <- calc_chr_enrichment(anno, "cis")
        trans_res <- calc_chr_enrichment(anno, "trans")
        results_list[[cancer]] <- rbind(cis_res, trans_res)
    }
    save(results_list,file="./g4_enrichment_results.RData")