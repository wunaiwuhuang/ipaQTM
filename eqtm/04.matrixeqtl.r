library(MatrixEQTL)
library(dplyr)
library(parallel)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/eqtm/")
files <- list.files("./data",pattern = "tpmuse.txt$")
cancers <- gsub("_tpmuse.txt$", "", files)

#cancers = c("TGCT","PRAD")
process_cancer <- function(cancer) {

  tpmloc=read.table(paste0("./data/",cancer,"_tpmloc.txt"),header=T,sep='\t')
  mdnaloc=read.table(paste0("./data/",cancer,"_mdnaloc.txt"),header=T,sep='\t')
  
  file1=c(paste0("./data/",cancer,"_mdnause.txt"))
  file2=c(paste0("./data/",cancer,"_tpmuse.txt"))
  file3=c(paste0("./data/",cancer,"_covariates_peer.txt"))

  snps = SlicedData$new();
  snps$fileDelimiter = "\t";     
  snps$fileOmitCharacters = "NA"; 
  snps$fileSkipRows = 1;          
  snps$fileSkipColumns = 1;       
  snps$fileSliceSize = 2000;      
  snps$LoadFile(file1)

  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      
  gene$fileOmitCharacters = "NA";
  gene$fileSkipRows = 1;          
  gene$fileSkipColumns = 1;       
  gene$fileSliceSize = 2000;     
  gene$LoadFile(file2)

  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      
  cvrt$fileOmitCharacters = "NA"; 
  cvrt$fileSkipRows = 1;          
  cvrt$fileSkipColumns = 1;       
  cvrt$fileSliceSize = 2000;      
  cvrt$LoadFile(file3)

  output_cis=paste0("./cis_sig/",cancer,'_CisSig_noex.txt')
  output_trans=paste0("./trans_sig/",cancer,'_TransSig_noex.txt')
  
  #a=as.numeric(nrow(snps))
  #b=as.numeric(nrow(gene))
  #p_trans = 0.05/(a*b)
  p_trans = 1e-6
  p_cis = 0.05

  me=Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_trans,
    pvOutputThreshold = p_trans,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = output_cis,
    pvOutputThreshold.cis = p_cis,
    snpspos = mdnaloc,
    genepos = tpmloc,
    cisDist = 1e6,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )

  dfFull = me$param$dfFull
  
  tstat_cis = me$cis$eqtls$statistic;
  r_cis = tstat_cis / sqrt( dfFull + tstat_cis^2 )
  cis=me$cis$eqtls
  cis$r=r_cis
  
  tstat_trans=me$trans$eqtls$statistic
  r_trans=tstat_trans / sqrt( dfFull + tstat_trans^2 )
  trans=me$trans$eqtls
  trans$r=r_trans

  write.table(cis,output_cis,sep='\t',row.names=F,quote=F)
  write.table(trans,output_trans,sep='\t',row.names=F,quote=F)
}

mclapply(cancers, process_cancer, mc.cores = 8)  # parallel processing