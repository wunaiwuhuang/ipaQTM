# for sp data Data modifying:
# 1.Impute (R package)
# 2.Quantile normalization
# 3.Inverse normal transformation

library(data.table)
library(impute)
library(preprocessCore)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/spqtm/data/")
#################################### deal with ipa data #####################################
    files_spuse <- list.files(pattern = "spuse.txt$")
    for (file in files_spuse) {
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