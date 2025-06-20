
library(data.table)
library(dplyr)
library(ggplot2)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/mediation/")
#------------------------ calculation -------------------
    cancers <- list.files(path="./data/",pattern = ".*_ipause.txt",full.names = F)
    cancers <- gsub("_ipause.txt", "", cancers)
    r_cutoffs <- seq(0.1, 0.3, by = 0.01)
    re_list<-list()
    for(cancer in cancers){
        cat(sprintf("start processing %s...\n", cancer))
        result<-data.frame(matrix(0,nrow=length(r_cutoffs),ncol=3,dimnames=list(NULL,c("r","path","gene"))))
        result$r <- r_cutoffs
        # trans ipaqtm data
        iqtm=read.table(paste0("./data/",cancer,"_trans_ipaqtm.txt"),header=T,sep="\t",check.names=F)
        # cis eqtm data
        eqtm=read.table(paste0("./data/",cancer,"_cis_eqtm.txt"),header=T,sep="\t",check.names=F)
        for(rval in r_cutoffs){
            e <- subset(eqtm,abs(eqtm$r)>=rval)
            i =iqtm[,1:2]%>%unique()
            e =e[,1:2]%>%unique()
            bridge=merge(e,i,by='snps')
            bridge=unique(bridge)
            colnames(bridge)=c('cg','gene','ipa')
            result[result$r==rval,"path"] <- dim(bridge)[1]
            result[result$r==rval,"gene"] <- length(unique(bridge$gene))
            result$cancer <- cancer
        }
        re_list[[cancer]] <- result
    }
    result <- do.call(rbind, re_list)
    result$plog<-log2(result$path)
    result$glog<-log2(result$gene)
    save(result, file = "./eqtm_rvalue_threshold.RData")
#------------------------ visualisation -------------------
    load("./eqtm_rvalue_threshold.RData")
    # order cancer types by the number of paths at r = 0.1
    cancer_order <- result %>%
    filter(r == 0.1) %>%
    arrange(desc(plog)) %>%
    pull(cancer)
    result$cancer <- factor(result$cancer, levels = cancer_order)
    # draw the plot
    p1 <- ggplot(result, aes(x = cancer, y = plog, group = factor(r), color = factor(r))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(title = "log2(Path count) vs. Cancer",
        x = "Cancer Type", y = "log2(Path count)", color = "r cutoff") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.border = element_rect(color = "black", fill = NA))
    p2 <- ggplot(result, aes(x = cancer, y = glog, group = factor(r), color = factor(r))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(title = "log2(Gene count) vs. Cancer",
        x = "Cancer Type", y = "log2(Gene count)", color = "r cutoff") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.border = element_rect(color = "black", fill = NA))

    ggsave(filename = "./02.pathnumber_alter.pdf", plot = p1, width = 10, height = 6)
    ggsave(filename = "./02.genenumber_alter.pdf", plot = p2, width = 10, height = 6)