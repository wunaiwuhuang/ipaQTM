library(jsonlite)
library(dplyr)
setwd("/data1/wuguojia/data/IPA_QTM_tcga/apa_reg/")
# deal with harmonizome data
    json_lines <- readLines("./harmonizome/all_polyadenylation.txt")
    attribute_gene_list <- list()
    for (line in json_lines) {
        if (trimws(line) == "") next
        tryCatch({
            # 解析 JSON
            data <- fromJSON(line)
            attr_name <- paste0(data$attribute$name,"-",data$dataset$name)
            gene_symbols <- data$associations$gene$symbol
            attribute_gene_list[[attr_name]] <- gene_symbols
        }, error = function(e) {
            message("Error parsing JSON: ", e$message)
            return(NULL)
        })
    }
    attribute_gene_list <- attribute_gene_list[!grepl("text-mining", names(attribute_gene_list), ignore.case = TRUE)] # remove text-mining
    harmonizome <- unlist(attribute_gene_list) %>% unique()
# deal with tranddb data
    files=list.files("./tranddb/14_APA_genes_whole_plus_dropouts",pattern = "*.txt",full.names = F)
    tranddb <- sub("^.*?_(.*?)\\.txt$", "\\1", files)

#save
    apa_reg <- unique(c(harmonizome,tranddb))
    save(apa_reg,harmonizome,tranddb,file = "./apa_regulators.RData")