library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures")

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