library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/figures/")

eqtm_cis <- list.files("../eqtm/cis_sig",pattern = "cisfilter_noex.txt$",full.names = T)
ipaqtm_cis <- list.files("../cis_sig",pattern = "cisfilter_noex.txt$",full.names = T)
eqtm_trans <- list.files("../eqtm/trans_sig",pattern = "transfilter_noex.txt$",full.names = T)
ipaqtm_trans <- list.files("../trans_sig",pattern = "transfilter_noex.txt$",full.names = T)

# create a dataframe to store the results
results <- data.frame(matrix(0, nrow = length(ipaqtm_cis), ncol = 4))
rownames(results) <- gsub("_cisfilter_noex.txt$", "", basename(ipaqtm_cis))
colnames(results) <- c("cisipa","cisoverlap","transipa","transoverlap")

# calculate the overlap for cis
for (i in 1:length(ipaqtm_cis)) {
    # get the cancer type
    cancer_type <- gsub("_cisfilter_noex.txt$", "", basename(ipaqtm_cis[i]))
    cat("Processing cancer: ", cancer_type , "\n")
    # read the data
    ipaqtm <- fread(ipaqtm_cis[i]) %>% as.data.frame()
    eqtm <- fread(eqtm_cis[i]) %>% as.data.frame()
    
    # extract gene symbol from ipaqtm
    ipaqtm$symbol <- sapply(strsplit(ipaqtm$gene, ":"), function(x) x[2])

    #calculate overlap cg (in share genes)
    common_genes <- intersect(ipaqtm$symbol, eqtm$symbol)
    all_overlap_cg <- character(0)
    for (gene in common_genes) {
        ipa_cg_gene <- unique(ipaqtm$snps[ipaqtm$symbol == gene])
        eqtm_cg_gene <- unique(eqtm$snps[eqtm$symbol == gene])
        overlap_cg_gene <- intersect(ipa_cg_gene, eqtm_cg_gene)
        all_overlap_cg <- c(all_overlap_cg, overlap_cg_gene)
    }
    overlap_cg <- unique(all_overlap_cg)
    # save results
    # results[cancer_type, "cisipa"] <- length(unique(ipaqtm$snps))
    results[cancer_type, "cisipa"] <- length(unique(ipaqtm$snps[ipaqtm$symbol %in% common_genes]))
    results[cancer_type, "cisoverlap"] <- length(overlap_cg)
}
# calculate the overlap for trans
for (i in 1:length(ipaqtm_trans)) {
    # get the cancer type
    cancer_type <- gsub("_transfilter_noex.txt$", "", basename(ipaqtm_trans[i]))
    cat("Processing cancer: ", cancer_type , "\n")
    # read the data
    ipaqtm <- fread(ipaqtm_trans[i]) %>% as.data.frame()
    eqtm <- fread(eqtm_trans[i]) %>% as.data.frame()

    # extract gene symbol from ipaqtm
    ipaqtm$symbol <- sapply(strsplit(ipaqtm$gene, ":"), function(x) x[2])

    #calculate overlap cg (in share genes)
    common_genes <- intersect(ipaqtm$symbol, eqtm$symbol)
    all_overlap_cg <- character(0)
    for (gene in common_genes) {
        ipa_cg_gene <- unique(ipaqtm$snps[ipaqtm$symbol == gene])
        eqtm_cg_gene <- unique(eqtm$snps[eqtm$symbol == gene])
        overlap_cg_gene <- intersect(ipa_cg_gene, eqtm_cg_gene)
        all_overlap_cg <- c(all_overlap_cg, overlap_cg_gene)
    }
    overlap_cg <- unique(all_overlap_cg)
    # save results
    # results[cancer_type, "transipa"] <- length(unique(ipaqtm$snps))
    results[cancer_type, "transipa"] <- length(unique(ipaqtm$snps[ipaqtm$symbol %in% common_genes]))
    results[cancer_type, "transoverlap"] <- length(overlap_cg)
}
# save the results
save(results, file = "./eqtm_overlap.RData")

load("./eqtm_overlap.RData")
# draw the figure
results$Cancer <- rownames(results)
# reshape cis 数据（long format）
cis_data <- results %>%
    dplyr::select(Cancer, cisipa, cisoverlap) %>%
    pivot_longer(cols = c(cisipa, cisoverlap), names_to = "Type", values_to = "Count") %>%
    mutate(Color = ifelse(Type == "cisipa", "lightblue", "blue"))
# cis_order <- results %>% arrange(desc(cisipa)) %>% pull(Cancer)
# cis_data$Cancer <- factor(cis_data$Cancer, levels = cis_order)
# reshape trans 数据
trans_data <- results %>%
    dplyr::select(Cancer, transipa, transoverlap) %>%
    pivot_longer(cols = c(transipa, transoverlap), names_to = "Type", values_to = "Count") %>%
    mutate(Color = ifelse(Type == "transipa", "lightcoral", "firebrick"))
# trans_order <- results %>% arrange(desc(transipa)) %>% pull(Cancer)
# trans_data$Cancer <- factor(trans_data$Cancer, levels = trans_order)
# Plot cis
p1 <- ggplot(cis_data, aes(x = Cancer, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_manual(values = c("cisipa" = "lightblue", "cisoverlap" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "cis IPA & Overlap", y = "Count", x = "")

# Plot trans
p2 <- ggplot(trans_data, aes(x = Cancer, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_manual(values = c("transipa" = "lightcoral", "transoverlap" = "firebrick")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "trans IPA & Overlap", y = "Count", x = "")

#save figure
ggsave(filename="./2.3.1.eqtm_overlap_cis.pdf",plot=p1,device="pdf",width=10,height=6,units="in",dpi=300)
ggsave(filename="./2.3.1.eqtm_overlap_trans.pdf",plot=p2,device="pdf",width=10,height=6,units="in",dpi=300)