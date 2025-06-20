library(data.table)
library(dplyr)
library(GenomicRanges)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/TFsite/")
source("../figures/2.1.0.distribute_enrich_fun.r") #load compute_enrichment_results function
# background probes
load("../00.probe_map_filtered.RData") # we set all considered probe as background    
probe_gr <- GRanges(seqnames = probe_map$chr,
                    ranges = IRanges(start = probe_map$start, end = probe_map$end),
                    id = probe_map$id)
# tf list
all_tf_files <- list.files("./all_tf", pattern = "*.txt", full.names = TRUE)
tf_names <- gsub("\\.txt$", "", basename(all_tf_files))
# universe parameters
cancers <- gsub("_cisfilter_noex.txt", "", list.files("/data1/wuguojia/data/IPA_QTM_tcga/cis_sig/", pattern = "cisfilter_noex.txt", full.names = FALSE))
types <- c("cis","trans")
sig_samcol <- "snps"
bg_samcol <- "id"
bg_regcol <- "TF"
# start calculation
cis_df <- list()
trans_df <- list()
for(i in seq_along(all_tf_files)){
  message("Processing TF: ", i, "/", length(all_tf_files))
  #---------- build tf annotation
  tf_file <- all_tf_files[i]
  tf_name <- tf_names[i]
  tf_data <- fread(tf_file, header = FALSE, data.table = FALSE)
  colnames(tf_data) <- c("chr", "start", "end")
  tf_data_gr <- GRanges(seqnames = tf_data$chr,
                ranges = IRanges(start = tf_data$start, end = tf_data$end))    
  anno <- probe_map
  anno[[bg_regcol]] <- NA #initiate as NA
  hits <- findOverlaps(probe_gr, tf_data_gr)
  anno[[bg_regcol]][queryHits(hits)] <- tf_name
  #---------- build tf annotation
  for(cancer in cancers){
    bg <- anno
    region <- tf_name    
    for (type in types) {
      sig_file <- paste0("../", type, "_sig/", cancer, "_", type, "filter_noex.txt")
      if (file.exists(sig_file)) {
        res <- build_enrich_matrix(
            cancer = cancer,
            type = type,
            sig_file = sig_file,
            sig_samcol = sig_samcol,
            bg = bg,
            bg_samcol = bg_samcol,
            bg_regcol = bg_regcol,
            region = region
        )
        if (type == "cis") {
          if (is.null(cis_df[[cancer]])) {
            cis_df[[cancer]] <- res
          } else {
            cis_df[[cancer]] <- rbind(cis_df[[cancer]], res)
          }
        } else {
          if (is.null(trans_df[[cancer]])) {
            trans_df[[cancer]] <- res
          } else {
            trans_df[[cancer]] <- rbind(trans_df[[cancer]], res)
          }
        }
      } else {
          message("File not found: ", sig_file)
      }
    }
  }
}

save(cis_df,trans_df, file = "./TF_cg_coverage_matrix.RData")
