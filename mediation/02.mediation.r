library(mediation)
library(data.table)
library(parallel)
library(dplyr)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/mediation/")
cancers <- list.files(path="./data/",pattern = ".*_ipause.txt",full.names = F)
cancers <- gsub("_ipause.txt", "", cancers)

process_cancer <- function(cancer) {
    cat(sprintf("start processing %s...\n", cancer))
    # trans ipaqtm data
    iqtm=read.table(paste0("./data/",cancer,"_trans_ipaqtm.txt"),header=T,sep="\t",check.names=F)
    # cis eqtm data
    eqtm=read.table(paste0("./data/",cancer,"_cis_eqtm.txt"),header=T,sep="\t",check.names=F)
    # ipausage data
    ipause=fread(paste0("./data/",cancer,"_ipause.txt"),header=T,sep="\t",check.names=F)
    ipause<- as.data.frame(ipause)
    rownames(ipause) <- ipause[,1]
    ipause=ipause[,-1]
    # mdnause data
    mdnause=fread(paste0("./data/",cancer,"_mdnause.txt"),header=T,sep="\t",check.names=F)
    mdnause<- as.data.frame(mdnause)
    rownames(mdnause) <- mdnause[,1]
    mdnause=mdnause[,-1]
    # expression data
    tpmuse=fread(paste0("./data/",cancer,"_tpmuse.txt"),header=T,sep="\t",check.names=F)
    tpmuse<- as.data.frame(tpmuse)
    rownames(tpmuse) <- tpmuse[,1]
    tpmuse=tpmuse[,-1]

    cocg <- intersect(unique(iqtm$snps),unique(eqtm$snps))
    iqtm <- iqtm[which(iqtm$snps %in% cocg),]
    eqtm <- eqtm[which(eqtm$snps %in% cocg),]
    iqtm=iqtm[,1:2]
    eqtm=eqtm[,1:2]
    bridge=merge(eqtm,iqtm,by='snps')
    bridge=unique(bridge)
    colnames(bridge)=c('cg','gene','ipa')
    # build result datafrome
    rtmp=data.frame(cpg=NA,ipa=NA,gene=NA,ACME_Estimate=NA,
                    ACME_pval=NA,ACME_CI_lower=NA,ACME_CI_Upper=NA,
                    ADE_Estimate=NA,ADE_pval=NA,ADE_CI_lowerr=NA,
                    ADE_CI_Upperr=NA,   total_effect_Estimater=NA,
                    total_effect_pvalr=NA,TE_CI_lower=NA,
                    TE_CI_Upper=NA, Prop.Mediated_Estimate=NA,
                    Prop.Mediated_pval=NA,Prop_CI_lower=NA, 
                    Prop_CI_Upper=NA,cancer_type=NA)    
    # start calculating mediation                    
    for(i in 1:dim(bridge)[1]){
        cg <- bridge[i,1]
        gene <- bridge[i,2]
        ipa <- bridge[i,3]
        tmp = cbind(mdnause[,cg],ipause[,ipa],tpmuse[,gene]) %>% as.data.frame()
        colnames(tmp)=c('Predictor','Outcome','Mediator')
        model_mediator <- lm(Mediator ~ Predictor, data = tmp)
        model_outcome <- lm(Outcome ~ Predictor + Mediator, data = tmp)
        set.seed(123.234)
        mediation_result <- mediate(model_mediator,model_outcome,treat ='Predictor',mediator = 'Mediator',bootstrap = 1000)
        rr1=summary(mediation_result)
        rtmp[i,1]=cg
        rtmp[i,2]=ipa
        rtmp[i,3]=gene
        rtmp[i,4]=rr1$d0
        rtmp[i,5]=rr1$d0.p
        rtmp[i,6]=rr1$d0.ci[1]
        rtmp[i,7]=rr1$d0.ci[2]
        rtmp[i,8]=rr1$z0
        rtmp[i,9]=rr1$z0.p
        rtmp[i,10]=rr1$z0.ci[1]
        rtmp[i,11]=rr1$z0.ci[2]
        rtmp[i,12]=rr1$tau.coef
        rtmp[i,13]=rr1$tau.p
        rtmp[i,14]=rr1$tau.ci[1]
        rtmp[i,15]=rr1$tau.ci[2]
        rtmp[i,16]=rr1$n.avg
        rtmp[i,17]=rr1$n.avg.p
        rtmp[i,18]=rr1$n.avg.ci[1]
        rtmp[i,19]=rr1$n.avg.ci[2]
        rtmp[i,20]=cancer
        # output progress
        if (i %% 1000 == 0) {
            cat(sprintf("Already done %d/%d in cancer %s\n", i, dim(bridge)[1],cancer))
        }
    }
    save(rtmp,file=paste0("./result/",cancer,"_mediate_all.RData"))
    cat(sprintf("All done in cancer %s \n",cancer))
}
# run the function in parallel
dir.create("./result/log", showWarnings = FALSE)  # 确保日志目录存在

wrapped_process_cancer <- function(cancer) {
    log_file <- paste0("./result/log/", cancer, "_log.txt")
    sink(log_file, split = TRUE)  # 重定向cat()输出到log文件，同时也输出到控制台
    on.exit(sink())               # 函数结束时关闭sink
    process_cancer(cancer)
}

mclapply(cancers, wrapped_process_cancer, mc.cores = 12)
