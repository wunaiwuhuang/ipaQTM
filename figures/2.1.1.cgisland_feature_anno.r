library(data.table)
library(dplyr)
library(rlang)
library(GenomicRanges)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

#download annotation files from https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz
anno<-fread("../00.HM450.hg38.manifest.gencode.v36.tsv")
anno<-unique(anno[,c("probeID","CGIposition")])
anno <- anno %>%
    mutate(CGIposition = case_when(
        CGIposition %in% c("N_Shore", "S_Shore") ~ "shore",
        CGIposition %in% c("N_Shelf", "S_Shelf") ~ "shelf",
        CGIposition == "Island" ~ "island",
        is.na(CGIposition) ~ "opensea",
        TRUE ~ CGIposition ))
save(anno, file = "./cgisland_feature_anno_1.Rdata")

#wangwenhui download cg island sites from web and use self-defined cg location. i copy cg sites from /data1/wangwenhui/genome_anno/cpg/cpgIslandExt.txt to /data1/wuguojia/data/IPA_QTM_tcga/00.cpgIslandExt.txt
anno<-fread("../00.HM450.hg38.manifest.gencode.v36.tsv")
anno<-na.omit(unique(anno[,c("probeID","CpG_chrm","CpG_beg","CpG_end")]))
cgi<-fread("../00.cpgIslandExt.txt")
cgi <- cgi[, c(2,3,4,5)]
colnames(cgi) <- c("chr", "start", "end", "CGi_ID")
anno_gr <- GRanges(seqnames = anno$CpG_chrm,
                ranges = IRanges(start = anno$CpG_beg, end = anno$CpG_end),
                probeID = anno$probeID)
cgi_gr <- GRanges(seqnames = cgi$chr,
                ranges = IRanges(start = cgi$start, end = cgi$end),
                CGi_ID = cgi$CGi_ID)                    
# shore/shelf (unstream or downstream 2kb / 4kb）
shore_up <- flank(cgi_gr, 2000, start = TRUE, both = TRUE)
shelf_up <- flank(cgi_gr, 4000, start = TRUE, both = TRUE)
shelf_up <- setdiff(shelf_up, shore_up)                    
anno$CGIposition <- "opensea"#initiate as opensea
# island
hits <- findOverlaps(anno_gr, cgi_gr)
anno$CGIposition[queryHits(hits)] <- "island"
# shore
hits <- findOverlaps(anno_gr, shore_up)
anno$CGIposition[queryHits(hits)] <- "shore"
# shelf
hits <- findOverlaps(anno_gr, shelf_up)
anno$CGIposition[queryHits(hits)] <- "shelf"
# island > shore > shelf
anno$CGIposition[queryHits(findOverlaps(anno_gr, cgi_gr))] <- "island"
anno$CGIposition[setdiff(queryHits(findOverlaps(anno_gr, shore_up)), 
                        queryHits(findOverlaps(anno_gr, cgi_gr)))] <- "shore"
anno$CGIposition[setdiff(queryHits(findOverlaps(anno_gr, shelf_up)),
                        union(queryHits(findOverlaps(anno_gr, cgi_gr)),
                            queryHits(findOverlaps(anno_gr, shore_up))))] <- "shelf"
# save as anno_2
save(anno, file = "./cgisland_feature_anno_2.Rdata")
