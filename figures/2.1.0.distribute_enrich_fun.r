# ------------------ build matrix -----------------------
    build_enrich_matrix <- function(cancer, type, sig_file, sig_samcol, bg , bg_samcol ,bg_regcol, region) { # cancer type, cis or trans, sig file pathway, sample column in sig files, background, sample column in background, feature column in background, the region should be analyze(all should include in bg$bg_samcol)
        stopifnot(file.exists(sig_file))
        stopifnot(type %in% c("cis", "trans"))

        sig_cg <- fread(sig_file)[[sig_samcol]] %>% unique() # sig
        bg$group <- ifelse(bg[[bg_samcol]] %in% sig_cg, "sig", "nosig")
        status <- c("sig_y", "sig_n", "nosig_y", "nosig_n") # colnames
        out <- matrix(0, nrow = length(region), ncol = length(status), dimnames = list(region, status)) #output matrix
        for (r in region) {
            in_region <- bg[[bg_regcol]] == r
            # sig 
            sig_in <- sum(bg$group == "sig" & in_region, na.rm = TRUE)
            sig_total <- sum(bg$group == "sig", na.rm = TRUE)
            sig_out <- sig_total - sig_in
            # nosig
            nosig_in <- sum(bg$group == "nosig" & in_region, na.rm = TRUE)
            nosig_total <- sum(bg$group == "nosig", na.rm = TRUE)
            nosig_out <- nosig_total - nosig_in
            # fill in matrix
            out[r, "sig_y"] <- sig_in
            out[r, "sig_n"] <- sig_out
            out[r, "nosig_y"] <- nosig_in
            out[r, "nosig_n"] <- nosig_out
        }
        out <- as.data.frame(out)
        out$cancer <- cancer
        out$type <- type
        return(out)
    }
# ------------------ calculate enrichment -----------------------
    compute_enrichment_results <- function(res_list, type = c("cis", "trans")) { #res_list usually refers to cis_df or trans_df after load(/data1/wuguojia/data/IPA_QTM_tcga/figures/*feature_data.RData)
        type <- match.arg(type)
        results <- list()
        for (cancer in names(res_list)) {
            df <- res_list[[cancer]]
            if (!all(c("sig_y", "sig_n", "nosig_y", "nosig_n") %in% colnames(df))) {
                warning(paste("Skipping", cancer, ": required columns missing"))
            next
            }
            cancer_results <- lapply(rownames(df), function(region) {
                sig_y <- df[region, "sig_y"]
                sig_n <- df[region, "sig_n"]
                nosig_y <- df[region, "nosig_y"]
                nosig_n <- df[region, "nosig_n"]
                # build matrix for calculation
                contingency <- matrix(c(sig_y, nosig_y, sig_n, nosig_n), nrow = 2, byrow = TRUE, dimnames = list(c("InRegion", "OutRegion"), c("Sig", "NoSig")))
                # Fisher test（one-sided greater）
                p_fisher <- tryCatch(fisher.test(contingency, alternative = "greater")$p.value, error = function(e) NA)
                # chip test
                p_chisq <- tryCatch(chisq.test(contingency)$p.value, error = function(e) NA)
                # G-test 
                p_gtest <- tryCatch({
                    observed <- as.vector(contingency)
                    expected <- outer(rowSums(contingency), colSums(contingency)) / sum(contingency)
                    G <- 2 * sum(ifelse(observed == 0, 0, observed * log(observed / expected)))
                    pchisq(G, df = 1, lower.tail = FALSE)
                }, error = function(e) NA)
                # OR calculation
                OR <- tryCatch((sig_y / sig_n) / (nosig_y / nosig_n), error = function(e) NA)
                data.frame(
                    Cancer = cancer,
                    Region = region,
                    P_value_F = p_fisher,
                    P_value_Chi = p_chisq,
                    P_value_G = p_gtest,
                    Odds_Ratio = OR,
                    Type = type,
                    stringsAsFactors = FALSE
                )
            })
            # save to result list
            results[[cancer]] <- do.call(rbind, cancer_results)
        }
        # combine all cancer results
        final_df <- do.call(rbind, results)
        # filter invalid OR
        final_df <- final_df[!is.na(final_df$Odds_Ratio) & is.finite(final_df$Odds_Ratio), ]
        rownames(final_df) <- NULL
        return(final_df)
    }
