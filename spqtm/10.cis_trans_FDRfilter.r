library(data.table)
library(dplyr)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/spqtm")
############################# for explaination files #############################
cisfiles <- list.files("./cis_sig/", pattern = "_CisSig\\.txt$")
transfiles <- list.files("./trans_sig/", pattern = "_TransSig\\.txt$")
#cis
for(file in cisfiles){
    print(paste0("deal with ",file))
    cancer<-basename(file)
    cancer<-gsub("_CisSig\\.txt$","",cancer)
    temp <- fread(paste0("./cis_sig/",cancer,"_CisSig.txt"))
    temp <- temp%>%filter(FDR<0.05)
    fwrite(temp,paste0("./cis_sig/",cancer,"_cisfilter.txt"),sep="\t",quote=FALSE)
}
#trans
for(file in transfiles){
    print(paste0("deal with ",file))
    cancer<-basename(file)
    cancer<-gsub("_TransSig\\.txt$","",cancer)
#    temp_mdna<-fread(paste0(cancer,"_mDNAuse.tsv"))
#    temp_mdna<-nrow(temp_mdna)
#    temp_ipa<-fread(paste0(cancer,"_IPAuse.txt"))
#    temp_ipa<-nrow(temp_ipa)
#    pthreshold<-1/(temp_mdna*temp_ipa)
    temp <- fread(paste0("./trans_sig/",cancer,"_TransSig.txt"))
#    temp <- temp%>%filter(`p-value`<pthreshold)
    temp <- temp%>%filter(FDR<0.05 & `p-value`<1e-8)
    fwrite(temp,paste0("./trans_sig/",cancer,"_transfilter.txt"),sep="\t",quote=FALSE)
}

############################# for no_explaination files ##########################
cisfiles <- list.files("./cis_sig/", pattern = "_CisSig_noex\\.txt$")
transfiles <- list.files("./trans_sig/", pattern = "_TransSig_noex\\.txt$")
#cis
for(file in cisfiles){
    print(paste0("deal with ",file))
    cancer<-basename(file)
    cancer<-gsub("_CisSig_noex\\.txt$","",cancer)
    temp <- fread(paste0("./cis_sig/",cancer,"_CisSig_noex.txt"))
    temp <- temp%>%filter(FDR<0.05)
    fwrite(temp,paste0("./cis_sig/",cancer,"_cisfilter_noex.txt"),sep="\t",quote=FALSE)
}
#trans
for(file in transfiles){
    print(paste0("deal with ",file))
    cancer<-basename(file)
    cancer<-gsub("_TransSig_noex\\.txt$","",cancer)
#    temp_mdna<-fread(paste0(cancer,"_mDNAuse.tsv"))
#    temp_mdna<-nrow(temp_mdna)
#    temp_ipa<-fread(paste0(cancer,"_IPAuse.txt"))
#    temp_ipa<-nrow(temp_ipa)
#    pthreshold<-1/(temp_mdna*temp_ipa)
    temp <- fread(paste0("./trans_sig/",cancer,"_TransSig_noex.txt"))
#    temp <- temp%>%filter(pvalue<pthreshold)
    temp <- temp%>%filter(FDR<0.05 & pvalue<1e-8)
    fwrite(temp,paste0("./trans_sig/",cancer,"_transfilter_noex.txt"),sep="\t",quote=FALSE)
}