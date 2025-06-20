library(mediation)
library(data.table)
library(parallel)
library(dplyr)

setwd("/data1/wuguojia/data/IPA_QTM_tcga/mediation/")
cancers <- list.files(path = "./data/", pattern = ".*_ipause.txt", full.names = FALSE)
cancers <- gsub("_ipause.txt", "", cancers)

process_cancer <- function(cancer) {
    log_file <- paste0("./result/log/", cancer, "_log.txt")
    cat(sprintf("Start processing %s at %s\n", cancer, Sys.time()), file = log_file, append = TRUE)
    tryCatch({
        # Load data
        iqtm <- fread(paste0("./data/", cancer, "_trans_ipaqtm.txt"), header = TRUE, sep = "\t")
        eqtm <- fread(paste0("./data/", cancer, "_cis_eqtm.txt"), header = TRUE, sep = "\t")
        ipause <- as.matrix(fread(paste0("./data/", cancer, "_ipause.txt"), header = TRUE, sep = "\t")[, -1, with = FALSE])
        mdnause <- as.matrix(fread(paste0("./data/", cancer, "_mdnause.txt"), header = TRUE, sep = "\t")[, -1, with = FALSE])
        tpmuse <- as.matrix(fread(paste0("./data/", cancer, "_tpmuse.txt"), header = TRUE, sep = "\t")[, -1, with = FALSE])
        # Filter eQTM data (|r| >= 0.15)
        eqtm <- eqtm[abs(r) >= 0.15]
        cocg <- intersect(unique(iqtm$snps), unique(eqtm$snps))
        # Build bridge data (CpG-gene-IPA triplets)
        bridge <- merge(
            eqtm[snps %in% cocg, .(snps, gene)],
            iqtm[snps %in% cocg, .(snps, gene)],
            by = "snps",
            allow.cartesian=TRUE
        )
        bridge <- unique(bridge)
        setnames(bridge, c("cg", "gene", "ipa"))
        # Split into 12 chunks for parallel processing
        bridge_chunks <- split(bridge, cut(1:nrow(bridge), 12, labels = FALSE))
        # Inner parallelization (12 cores)
        results_list <- mclapply(bridge_chunks, function(chunk) {
            res <- vector("list", length = nrow(chunk))
            for (i in seq_len(nrow(chunk))) {
                cg <- as.character(chunk[i, 1])
                gene <- as.character(chunk[i, 2])
                ipa <- as.character(chunk[i, 3])
                # Build mediation model
                tmp <- data.frame(
                    Predictor = mdnause[, cg],
                    Outcome = ipause[, ipa],
                    Mediator = tpmuse[, gene]
                )
                #---------------- check before mediation
                fit1 <- lm(Outcome ~ Predictor, data = tmp) # cg to ipa
                fit1_coef <- coef(summary(fit1))["Predictor", "Estimate"]
                fit1_p <- coef(summary(fit1))["Predictor", "Pr(>|t|)"]                
                fit2 <- lm(Mediator ~ Predictor, data = tmp) # cg to gene
                fit2_coef <- coef(summary(fit2))["Predictor", "Estimate"]
                fit2_p <- coef(summary(fit2))["Predictor", "Pr(>|t|)"]                
                fit3 <- lm(Outcome ~ Mediator, data = tmp) # gene to ipa        
                fit3_coef <- coef(summary(fit3))["Mediator", "Estimate"]
                fit3_p <- coef(summary(fit3))["Mediator", "Pr(>|t|)"]
                # cg to ipa direction should keep same with gene involved direction
                direction_ok <- sign(fit1_coef) == sign(fit2_coef * fit3_coef)
                # pvalue should all within 0.05
                pvalue_ok <- fit1_p <= 0.05 && fit2_p <= 0.05 && fit3_p <= 0.05
                #---------------- start mediation
                if(pvalue_ok && direction_ok){
                    model_mediator <- lm(Mediator ~ Predictor, data = tmp)
                    model_outcome <- lm(Outcome ~ Predictor + Mediator, data =tmp)                
                    set.seed(20030622 + i)  # Ensure reproducibility
                    mediation_result <- mediate(
                        model_mediator, model_outcome,
                        treat = "Predictor", mediator = "Mediator",
                        bootstrap = 1000
                    )
                    rr1 <- summary(mediation_result)
                    # Log progress every 1000 iterations
                    if (i %% 1000 == 0) {cat(sprintf("Done %d/%d in %s\n", i, nrow(chunk), cancer), file = log_file, append = TRUE)}
                    # Store results
                    res[[i]] <- data.table(
                        cpg = cg, ipa = ipa, gene = gene,
                        ACME_Estimate = rr1$d0, ACME_pval = rr1$d0.p,
                        ACME_CI_lower = rr1$d0.ci[1], ACME_CI_Upper = rr1$d0.ci[2],
                        ADE_Estimate = rr1$z0, ADE_pval = rr1$z0.p,
                        ADE_CI_lowerr = rr1$z0.ci[1], ADE_CI_Upperr = rr1$z0.ci[2],
                        total_effect_Estimater = rr1$tau.coef, total_effect_pvalr = rr1$tau.p,
                        TE_CI_lower = rr1$tau.ci[1], TE_CI_Upper = rr1$tau.ci[2],
                        Prop.Mediated_Estimate = rr1$n.avg, Prop.Mediated_pval = rr1$n.avg.p,
                        Prop_CI_lower = rr1$n.avg.ci[1], Prop_CI_Upper = rr1$n.avg.ci[2],
                        cancer_type = cancer)
                }else{res[[i]] <- NULL}
            }
            rbindlist(res)  # Combine results for this chunk
        }, mc.cores = 12)  # Use 12 cores for inner parallelization
        # Combine all results and save
        rtmp <- rbindlist(results_list)
        save(rtmp, file = paste0("./result/", cancer, "_mediate_all.RData"))
        cat(sprintf("Finished %s at %s\n", cancer, Sys.time()), file = log_file, append = TRUE)
    }, error = function(e) {
        cat(sprintf("ERROR in %s: %s\n", cancer, e$message), file = log_file, append = TRUE)
    })
}

# Outer loop (serial processing)
lapply(cancers, process_cancer)