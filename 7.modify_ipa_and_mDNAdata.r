# for ipa data Data modifying:
# 1.Impute (R package)
# 2.Quantile normalization
# 3.Inverse normal transformation
# for mDNA data Data modifying:
# 1.Impute (R package)

library(data.table)
library(impute)
library(preprocessCore)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/data/")
#################################### deal with ipa data #####################################
    files_ipause <- list.files(pattern = "ipause.txt$")
    for (file in files_ipause) {
        cat("Processing IPAuse file:", file, "\n")
        data <- fread(file)

        # knn function to impute missing values
        imputed_values <- impute.knn(as.matrix(data[, -1, with = FALSE]))$data
        data[, 2:ncol(data)] <- as.data.table(imputed_values)
        
        # Quantile normalization (for all samples)
        mat <- as.matrix(data[, 2:ncol(data), with = FALSE])
        mat <- preprocessCore::normalize.quantiles(mat)

        # Inverse normal transformation (for IPA event)
        mat <- t(apply(mat, 1, function(line){
        qnorm((rank(line) - 0.5) / length(line))
        }))

        # Convert back to data
        data[, 2:ncol(data)] <- as.data.table(mat)

        # Save the modified data
        fwrite(data, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }

################################### deal with mDNA data #####################################
    files_mdna <- list.files(pattern = "mdnause.txt$")
    for (file in files_mdna) {
        cat("Processing mDNAuse file:", file, "\n")
        data <- fread(file)

        # knn function to impute missing values
        imputed_values <- impute.knn(as.matrix(data[, -1, with = FALSE]))$data
        data[, 2:ncol(data)] <- as.data.table(imputed_values)

        # Save the modified data
        fwrite(data, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
