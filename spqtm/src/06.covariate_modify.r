library(data.table)
library(dplyr)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/spqtm/data/")
files <- list.files(pattern = "spuse.txt$")
cancers <- gsub("_spuse.txt$", "", files)

clinical <- read.table("/data1/wuguojia/data/IPA_QTM_tcga/00.all.info.txt",header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
clinical$id <- rownames(clinical)
clinical <- clinical[,c("id","Tumour","gender","age")] #male is 1, female is 0

for(cancer in cancers){
    ###################################
    # deal with covariate data
    print(paste0("Processing ", cancer," covariate and clinical data"))
    # Read the clinical data
    cli <- subset(clinical, Tumour == cancer)
    # Read the covariate data
    covariate <- fread(paste0(cancer, "_covariates.txt"))
    covariate <- as.data.frame(t(covariate[,-1]))
    covariate$id <- rownames(covariate)
    data <- merge(cli, covariate, by = "id", all.x = TRUE) #must be 100% matching
    rownames(data) <- data$id
    # reorder the data to match the samples in spuse
    spuse <- fread(paste0(cancer, "_spuse.txt"))
    samples <- colnames(spuse)[-1]  
    data <- data[match(samples, data$id),]
    # Write the merged data to a new file
    data <- data[, !(colnames(data) %in% c("id", "Tumour"))] %>% t()

    # add PC_1 - PC_5 to cater for subsequent analysis, set to same value
    # give artificial values, very small and avoid multicollinearity
    data <- as.data.frame(t(data))
    n <- nrow(data)
    set.seed(123)
    data$PC_1 <- runif(n, 0.003, 0.005)
    data$PC_2 <- seq(0.003, 0.005, length.out = n)
    data$PC_3 <- rnorm(n, mean = 0.003, sd = 0.0005)
    data$PC_4 <- log1p(data$PC_2)
    data$PC_5 <- sin(seq(0, pi/6, length.out = n)) * 0.001
    #reorder the columns
    original_names <- names(data)
    pc_names <- paste0("PC_", 1:5)
    pc_data <- data[, pc_names]
    data <- data[, !colnames(data)%in%pc_names]
    data <- cbind(data[, 1:2], pc_data, data[, 3:ncol(data)])
    data <- as.data.frame(t(data))
    #add peer name
    n <- nrow(data)-7
    rownames(data)[1:2]<-c("gender","age")
    rownames(data)[3:7]<-paste0("PC_",1:5)
    rownames(data)[8:(n+7)]<-paste0("peer_",1:n)
    # save the data
    fwrite(data, paste0(cancer, "_covariates.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
    data_peer <- data[-(1:7),]# only retain peer data
    fwrite(data_peer, paste0(cancer, "_covariates_peer.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
}
