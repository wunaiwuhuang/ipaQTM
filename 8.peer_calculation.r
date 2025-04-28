# enter conda activate peer
# peer factors (N<150 use 15 peers, 150<=N<250 use 30 peers ,N>=250 use 35 peers)
# reference: https://github.com/PMBio/peer/wiki/Tutorial

library(peer)
library(data.table)
library(parallel)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/data/")
files <- list.files(pattern = "ipause.txt$")

process_peer_file <- function(file){
    cat("Processing file:", file, "\n")
    data <- fread(file, header = T, sep = "\t", data.table = F)
    expr <- as.matrix(data[, -1]) # Remove the first column (gene names)
    expr <- t(expr) # N rows and G columns, where N is the number of samples, and G is the number of genes

    # Get the number of samples
    num_samples <- nrow(expr)
    
    # Determine the number of factors based on the number of samples
    if (num_samples < 150) {
        num_factors <- 15
    } else if (num_samples < 250) {
        num_factors <- 30
    } else {
        num_factors <- 35
    }
    
    # Create a PEER object
    model <- PEER()
    
    # set the observed data,
    PEER_setPhenoMean(model,as.matrix(expr)) #output NULL response means no error here

    # Set the number of factors
    PEER_setNk(model, num_factors) 
    
    # Run the PEER model 
    PEER_update(model)
    
    # Get the factor matrix : inferred confounders (NxK matrix)
    factors <- PEER_getX(model)
    
    # add sample names to the row names
    rownames(factors) <- rownames(expr)
    factors <- t(as.data.frame(factors))

    # Save the factors to a file
    cancer <- gsub("_ipause.txt", "", file)
    fwrite(factors, paste0(cancer,"_covariates.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
}
# Process each file in parallel
mclapply(files, process_peer_file, mc.cores = 12)